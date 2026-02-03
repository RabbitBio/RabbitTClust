#include "dbscan.h"
#include "phmap.h"   // For phmap::flat_hash_map
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <omp.h>
#include <algorithm>
#include <cmath>

using std::cerr;
using std::endl;
using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::queue;

/**
 * @brief Calculate Jaccard similarity between two hash arrays
 * Implemented locally for DBSCAN independence
 */
static double jaccard(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry) {
    uint64_t common_elements = 0;
    uint64_t i = 0, j = 0;
    uint32_t sizeRef = hashesRef.size();
    uint32_t sizeQry = hashesQry.size();
    
    while (i < sizeRef && j < sizeQry) {
        if (hashesRef[i] < hashesQry[j]) {
            i++;
        } else if (hashesRef[i] > hashesQry[j]) {
            j++;
        } else {
            common_elements++;
            i++;
            j++;
        }
    }
    uint64_t union_elements = static_cast<uint64_t>(sizeRef) + sizeQry - common_elements;
    
    if (union_elements == 0) {
        return 0.0; 
    }
    
    return static_cast<double>(common_elements) / union_elements;
}

/**
 * @brief Calculate Mash distance from Jaccard similarity
 * Implemented locally for DBSCAN independence
 */
static double calculate_mash_distance(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry, double kmerSize) {
    double jaccard_ = jaccard(hashesRef, hashesQry);
    
    if (jaccard_ == 1.0) {
        return 0.0;
    }
    
    double dist = -log(2.0 * jaccard_ / (1.0 + jaccard_)) / kmerSize;
    
    if (dist > 1.0) {
        dist = 1.0;
    }
    
    return dist;
}

/**
 * @brief Build inverted index for all sketches (KSSD, 32-bit hashes)
 * Optimization 7: Skip empty sketches and ensure hash uniqueness
 */
phmap::flat_hash_map<uint32_t, vector<int>> buildInvertedIndex32(
    const vector<KssdSketchInfo>& sketches
) {
    phmap::flat_hash_map<uint32_t, vector<int>> index;
    index.reserve(10000000);  // Pre-allocate for large datasets
    
    for(size_t i = 0; i < sketches.size(); i++) {
        if(sketches[i].use64) continue;  // Skip 64-bit hashes for now
        
        // Optimization 7: Skip empty sketches
        if(sketches[i].hash32_arr.empty()) continue;
        
        // Optimization 7: Ensure hash uniqueness (should be done at sketch generation, but double-check)
        // Note: We assume hash32_arr is already sorted and unique from sketch generation
        // If not, we would need to sort and unique here, but that's expensive
        
        for(uint32_t hash : sketches[i].hash32_arr) {
            index[hash].push_back(i);
        }
    }
    
    return index;
}

/**
 * @brief Find all points within epsilon distance using inverted index (KSSD, 32-bit)
 * Optimized version using global stamp+count arrays and direct Jaccard calculation
 */
vector<int> findNeighborsKSSDWithIndex(
    const vector<KssdSketchInfo>& sketches,
    int point_idx,
    double jaccard_min,           // Pre-calculated jaccard_min (optimization 5)
    const phmap::flat_hash_map<uint32_t, vector<int>>& inverted_index,
    vector<uint32_t>& mark,      // Global mark array (epoch-based)
    vector<uint32_t>& cnt,        // Global count array
    vector<int>& touched,         // Global touched list
    uint32_t& cur_mark,          // Current epoch marker
    vector<vector<int>>& thread_neighbors,  // Optimization 2: Reuse thread-local buffers
    int threads
) {
    const auto& point = sketches[point_idx];
    
    // Optimization 6: Handle cur_mark overflow
    cur_mark++;
    if(cur_mark == 0) {
        std::fill(mark.begin(), mark.end(), 0);
        cur_mark = 1;
    }
    
    if(point.use64) {
        // Fallback to brute force for 64-bit hashes
        // Clear and reuse thread_neighbors
        for(auto& v : thread_neighbors) v.clear();
        
        #pragma omp parallel for schedule(dynamic, 100) num_threads(threads)
        for(size_t i = 0; i < sketches.size(); i++) {
            if(i == point_idx || !sketches[i].use64) continue;
            
            int tid = omp_get_thread_num();
            auto& buf = thread_neighbors[tid];
            
            // Compute distance for 64-bit hashes
            uint64_t common = 0;
            size_t i1 = 0, i2 = 0;
            size_t size1 = point.hash64_arr.size();
            size_t size2 = sketches[i].hash64_arr.size();
            
            while(i1 < size1 && i2 < size2) {
                if(point.hash64_arr[i1] < sketches[i].hash64_arr[i2]) {
                    i1++;
                } else if(point.hash64_arr[i1] > sketches[i].hash64_arr[i2]) {
                    i2++;
                } else {
                    common++;
                    i1++;
                    i2++;
                }
            }
            
            uint64_t union_size = size1 + size2 - common;
            double jaccard_ = (union_size == 0) ? 0.0 : (double)common / union_size;
            
            if(jaccard_ >= jaccard_min) {
                buf.push_back(i);
            }
        }
        
        // Merge thread-local results
        vector<int> neighbors;
        for(int t = 0; t < threads; t++) {
            neighbors.insert(neighbors.end(), thread_neighbors[t].begin(), thread_neighbors[t].end());
        }
        return neighbors;
    }
    
    // Use inverted index for 32-bit hashes with optimized stamp+count arrays
    touched.clear();
    
    int sizeRef = point.hash32_arr.size();
    if(sizeRef == 0) return vector<int>();  // Optimization 7: Skip empty sketches
    
    // Step 1: Count intersections using inverted index (no hash map allocation)
    for(uint32_t hash : point.hash32_arr) {
        auto it = inverted_index.find(hash);
        if(it != inverted_index.end()) {
            for(int candidate_id : it->second) {
                if(candidate_id == point_idx) continue;
                
                if(mark[candidate_id] != cur_mark) {
                    // First time seeing this candidate in current query
                    mark[candidate_id] = cur_mark;
                    cnt[candidate_id] = 1;
                    touched.push_back(candidate_id);
                } else {
                    // Already seen, increment count
                    cnt[candidate_id]++;
                }
            }
        }
    }
    
    // Optimization 3: Skip OpenMP for small touched sets
    if(touched.empty()) {
        return vector<int>();
    }
    
    // Step 2: Parallel evaluation using thread-local vectors (no critical section)
    // Optimization 2: Clear and reuse thread_neighbors
    for(auto& v : thread_neighbors) {
        v.clear();
        v.reserve(256);  // Pre-reserve
    }
    
    // Optimization 3: Threshold-based parallelization
    const size_t PARALLEL_THRESHOLD = 5000;
    if(touched.size() < PARALLEL_THRESHOLD || threads <= 1) {
        // Serial processing for small sets
        vector<int> neighbors;
        neighbors.reserve(touched.size());
        
        for(int candidate_id : touched) {
            const auto& candidate = sketches[candidate_id];
            if(candidate.use64) continue;
            
            int common = cnt[candidate_id];
            int sizeQry = candidate.hash32_arr.size();
            
            // Optimization 4: Use ceil for min_common_needed
            double v = jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min);
            int min_common_needed = (int)std::ceil(v - 1e-12);
            if(common < min_common_needed) {
                continue;
            }
            
            // Direct Jaccard calculation from common
            uint64_t union_size = sizeRef + sizeQry - common;
            double jaccard_ = (union_size == 0) ? 0.0 : (double)common / union_size;
            
            if(jaccard_ >= jaccard_min) {
                neighbors.push_back(candidate_id);
            }
        }
        return neighbors;
    }
    
    // Parallel processing for large sets
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        auto& buf = thread_neighbors[tid];
        
        #pragma omp for schedule(guided)
        for(size_t i = 0; i < touched.size(); i++) {
            int candidate_id = touched[i];
            
            const auto& candidate = sketches[candidate_id];
            if(candidate.use64) continue;
            
            int common = cnt[candidate_id];
            int sizeQry = candidate.hash32_arr.size();
            
            // Optimization 4: Use ceil for min_common_needed
            double v = jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min);
            int min_common_needed = (int)std::ceil(v - 1e-12);
            if(common < min_common_needed) {
                continue;
            }
            
            // Direct Jaccard calculation from common
            uint64_t union_size = sizeRef + sizeQry - common;
            double jaccard_ = (union_size == 0) ? 0.0 : (double)common / union_size;
            
            if(jaccard_ >= jaccard_min) {
                buf.push_back(candidate_id);
            }
        }
    }
    
    // Merge thread-local results
    vector<int> neighbors;
    neighbors.reserve(touched.size());
    for(int t = 0; t < threads; t++) {
        neighbors.insert(neighbors.end(), thread_neighbors[t].begin(), thread_neighbors[t].end());
    }
    
    return neighbors;
}

/**
 * @brief Find all points within epsilon distance of a given point (KSSD)
 * @deprecated Use findNeighborsKSSDWithIndex for better performance
 */
vector<int> findNeighborsKSSD(
    const vector<KssdSketchInfo>& sketches,
    int point_idx,
    double eps,
    int kmer_size,
    int threads
) {
    vector<int> neighbors;
    const auto& point = sketches[point_idx];
    
    #pragma omp parallel for schedule(dynamic, 100) num_threads(threads)
    for(size_t i = 0; i < sketches.size(); i++) {
        if(i == point_idx) continue;
        
        double dist;
        if(point.use64 && sketches[i].use64) {
            // Both use 64-bit hashes - need to compute jaccard manually
            uint64_t common = 0;
            size_t i1 = 0, i2 = 0;
            size_t size1 = point.hash64_arr.size();
            size_t size2 = sketches[i].hash64_arr.size();
            
            while(i1 < size1 && i2 < size2) {
                if(point.hash64_arr[i1] < sketches[i].hash64_arr[i2]) {
                    i1++;
                } else if(point.hash64_arr[i1] > sketches[i].hash64_arr[i2]) {
                    i2++;
                } else {
                    common++;
                    i1++;
                    i2++;
                }
            }
            
            uint64_t union_size = size1 + size2 - common;
            double jaccard_ = (union_size == 0) ? 0.0 : (double)common / union_size;
            
            if(jaccard_ == 1.0) {
                dist = 0.0;
            } else {
                dist = -log(2.0 * jaccard_ / (1.0 + jaccard_)) / kmer_size;
                if(dist > 1.0) dist = 1.0;
            }
        } else if(!point.use64 && !sketches[i].use64) {
            // Both use 32-bit hashes
            dist = calculate_mash_distance(point.hash32_arr, sketches[i].hash32_arr, kmer_size);
        } else {
            // Different hash sizes - cannot compare directly
            continue;
        }
        
        if(dist <= eps) {
            #pragma omp critical
            {
                neighbors.push_back(i);
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Find all points within epsilon distance of a given point (MinHash)
 */
vector<int> findNeighborsMinHash(
    const vector<SketchInfo>& sketches,
    int point_idx,
    double eps,
    int sketch_func_id,
    int threads
) {
    vector<int> neighbors;
    const auto& point = sketches[point_idx];
    
    #pragma omp parallel for schedule(dynamic, 100) num_threads(threads)
    for(size_t i = 0; i < sketches.size(); i++) {
        if(i == point_idx) continue;
        
        double dist;
        if(sketch_func_id == 0) {
            if(point.isContainment) {
                dist = point.minHash->containDistance(sketches[i].minHash);
            } else {
                dist = point.minHash->distance(sketches[i].minHash);
            }
        } else {
            // Unsupported sketch type
            continue;
        }
        
        if(dist <= eps) {
            #pragma omp critical
            {
                neighbors.push_back(i);
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief DBSCAN clustering for KSSD sketches
 */
DBSCANResult KssdDBSCAN(
    vector<KssdSketchInfo>& sketches,
    double eps,
    int minPts,
    int kmer_size,
    int threads
) {
    DBSCANResult result;
    int n = sketches.size();
    
    if(n == 0) {
        result.num_clusters = 0;
        result.num_noise = 0;
        return result;
    }
    
    // Initialize labels: -1 = unvisited, -2 = noise, >= 0 = cluster ID
    vector<int> labels(n, -1);
    int cluster_id = 0;
    
    cerr << "-----Running DBSCAN clustering (KSSD)..." << endl;
    cerr << "-----Parameters: eps=" << eps << ", minPts=" << minPts << endl;
    
    // Optimization 5: Pre-calculate jaccard_min once
    double x = std::exp(-eps * kmer_size);
    double jaccard_min = x / (2.0 - x);
    
    // Build inverted index for acceleration (only for 32-bit hashes)
    cerr << "-----Building inverted index for acceleration..." << endl;
    phmap::flat_hash_map<uint32_t, vector<int>> inverted_index = buildInvertedIndex32(sketches);
    cerr << "-----Inverted index built: " << inverted_index.size() << " unique hashes" << endl;
    
    // Initialize global stamp+count arrays (optimization 1: avoid per-query allocation)
    vector<uint32_t> mark(n, 0);      // Epoch-based marking
    vector<uint32_t> cnt(n, 0);       // Intersection counts
    vector<int> touched;               // Touched candidates list
    touched.reserve(10000);            // Pre-allocate
    uint32_t cur_mark = 0;            // Current epoch marker
    
    // Optimization 1: Use stamp array for queue tracking (instead of vector<bool>)
    vector<uint32_t> qmark(n, 0);      // Queue membership stamp
    uint32_t qepoch = 0;               // Queue epoch marker
    
    // Optimization 2: Pre-allocate thread-local buffers
    vector<vector<int>> thread_neighbors(threads);
    for(auto& v : thread_neighbors) {
        v.reserve(256);
    }
    
    // Process each point
    for(int i = 0; i < n; i++) {
        if(labels[i] != -1) continue;  // Already processed
        
        // Find neighbors using inverted index with optimized arrays
        vector<int> neighbors = findNeighborsKSSDWithIndex(
            sketches, i, jaccard_min, inverted_index, 
            mark, cnt, touched, cur_mark, thread_neighbors, threads
        );
        
        // Optimization 4: minPts includes the point itself
        if((int)neighbors.size() + 1 < minPts) {
            // Not a core point, mark as noise (may be reassigned later)
            labels[i] = -2;
            continue;
        }
        
        // Core point found, start a new cluster
        labels[i] = cluster_id;
        queue<int> seed_set;
        
        // Optimization 1: Use stamp array for queue tracking
        qepoch++;
        if(qepoch == 0) {
            std::fill(qmark.begin(), qmark.end(), 0);
            qepoch = 1;
        }
        
        auto inq = [&](int v) { return qmark[v] == qepoch; };
        auto set_inq = [&](int v) { qmark[v] = qepoch; };
        
        for(int neighbor : neighbors) {
            if(!inq(neighbor)) {
                seed_set.push(neighbor);
                set_inq(neighbor);
            }
        }
        
        // Expand cluster using seed set
        while(!seed_set.empty()) {
            int q = seed_set.front();
            seed_set.pop();
            
            if(labels[q] == -2) {
                // Reassign noise point to cluster
                labels[q] = cluster_id;
            } else if(labels[q] != -1) {
                // Already processed
                continue;
            }
            
            labels[q] = cluster_id;
            
            // Find neighbors of q using inverted index
            vector<int> q_neighbors = findNeighborsKSSDWithIndex(
                sketches, q, jaccard_min, inverted_index,
                mark, cnt, touched, cur_mark, thread_neighbors, threads
            );
            
            // Optimization 4: minPts includes the point itself
            if((int)q_neighbors.size() + 1 >= minPts) {
                // q is also a core point, add its neighbors to seed set (with deduplication)
                for(int neighbor : q_neighbors) {
                    if((labels[neighbor] == -1 || labels[neighbor] == -2) && !inq(neighbor)) {
                        seed_set.push(neighbor);
                        set_inq(neighbor);
                    }
                }
            }
        }
        
        cluster_id++;
        
        if((i + 1) % 1000 == 0) {
            cerr << "-----Processed " << (i + 1) << " / " << n << " points" << endl;
        }
    }
    
    // Organize results
    vector<vector<int>> clusters(cluster_id);
    vector<int> noise;
    
    for(int i = 0; i < n; i++) {
        if(labels[i] == -2) {
            noise.push_back(i);
        } else if(labels[i] >= 0) {
            clusters[labels[i]].push_back(i);
        }
    }
    
    // Remove empty clusters
    clusters.erase(
        std::remove_if(clusters.begin(), clusters.end(),
            [](const vector<int>& c) { return c.empty(); }),
        clusters.end()
    );
    
    result.clusters = clusters;
    result.noise = noise;
    result.num_clusters = clusters.size();
    result.num_noise = noise.size();
    
    cerr << "-----DBSCAN clustering complete!" << endl;
    cerr << "-----Found " << result.num_clusters << " clusters" << endl;
    cerr << "-----Found " << result.num_noise << " noise points (outliers)" << endl;
    
    return result;
}

/**
 * @brief DBSCAN clustering for MinHash sketches
 */
DBSCANResult MinHashDBSCAN(
    vector<SketchInfo>& sketches,
    double eps,
    int minPts,
    int sketch_func_id,
    int threads
) {
    DBSCANResult result;
    int n = sketches.size();
    
    if(n == 0) {
        result.num_clusters = 0;
        result.num_noise = 0;
        return result;
    }
    
    // Initialize labels: -1 = unvisited, -2 = noise, >= 0 = cluster ID
    vector<int> labels(n, -1);
    int cluster_id = 0;
    
    cerr << "-----Running DBSCAN clustering (MinHash)..." << endl;
    cerr << "-----Parameters: eps=" << eps << ", minPts=" << minPts << endl;
    
    // Process each point
    for(int i = 0; i < n; i++) {
        if(labels[i] != -1) continue;  // Already processed
        
        // Find neighbors
        vector<int> neighbors = findNeighborsMinHash(sketches, i, eps, sketch_func_id, threads);
        
        if(neighbors.size() < minPts) {
            // Not a core point, mark as noise (may be reassigned later)
            labels[i] = -2;
            continue;
        }
        
        // Core point found, start a new cluster
        labels[i] = cluster_id;
        queue<int> seed_set;
        for(int neighbor : neighbors) {
            seed_set.push(neighbor);
        }
        
        // Expand cluster using seed set
        while(!seed_set.empty()) {
            int q = seed_set.front();
            seed_set.pop();
            
            if(labels[q] == -2) {
                // Reassign noise point to cluster
                labels[q] = cluster_id;
            } else if(labels[q] != -1) {
                // Already processed
                continue;
            }
            
            labels[q] = cluster_id;
            
            // Find neighbors of q
            vector<int> q_neighbors = findNeighborsMinHash(sketches, q, eps, sketch_func_id, threads);
            
            if(q_neighbors.size() >= minPts) {
                // q is also a core point, add its neighbors to seed set
                for(int neighbor : q_neighbors) {
                    if(labels[neighbor] == -1 || labels[neighbor] == -2) {
                        seed_set.push(neighbor);
                    }
                }
            }
        }
        
        cluster_id++;
        
        if((i + 1) % 1000 == 0) {
            cerr << "-----Processed " << (i + 1) << " / " << n << " points" << endl;
        }
    }
    
    // Organize results
    vector<vector<int>> clusters(cluster_id);
    vector<int> noise;
    
    for(int i = 0; i < n; i++) {
        if(labels[i] == -2) {
            noise.push_back(i);
        } else if(labels[i] >= 0) {
            clusters[labels[i]].push_back(i);
        }
    }
    
    // Remove empty clusters
    clusters.erase(
        std::remove_if(clusters.begin(), clusters.end(),
            [](const vector<int>& c) { return c.empty(); }),
        clusters.end()
    );
    
    result.clusters = clusters;
    result.noise = noise;
    result.num_clusters = clusters.size();
    result.num_noise = noise.size();
    
    cerr << "-----DBSCAN clustering complete!" << endl;
    cerr << "-----Found " << result.num_clusters << " clusters" << endl;
    cerr << "-----Found " << result.num_noise << " noise points (outliers)" << endl;
    
    return result;
}

/**
 * @brief Print DBSCAN results for MinHash sketches
 */
void printDBSCANResult(
    const DBSCANResult& result,
    const vector<SketchInfo>& sketches,
    bool sketch_by_file,
    const string& output_file,
    double eps,
    int minPts
) {
    FILE *fp = fopen(output_file.c_str(), "w");
    if(!fp) {
        cerr << "Error in printDBSCANResult(), cannot open file: " << output_file << endl;
        exit(1);
    }
    
    fprintf(fp, "# DBSCAN Clustering Results\n");
    fprintf(fp, "# Parameters: eps=%.6f, minPts=%d\n", eps, minPts);
    fprintf(fp, "# Total clusters: %d\n", result.num_clusters);
    fprintf(fp, "# Total noise points (outliers): %d\n", result.num_noise);
    fprintf(fp, "#\n");
    
    // Print clusters
    for(size_t i = 0; i < result.clusters.size(); i++) {
        fprintf(fp, "# Cluster %zu (size: %zu)\n", i, result.clusters[i].size());
        for(int idx : result.clusters[i]) {
            if(sketch_by_file) {
                fprintf(fp, "%s\n", sketches[idx].fileName.c_str());
            } else {
                fprintf(fp, "%s\n", sketches[idx].seqInfo.name.c_str());
            }
        }
        fprintf(fp, "\n");
    }
    
    // Print noise points (outliers)
    if(!result.noise.empty()) {
        fprintf(fp, "# Noise points (outliers) - %d points\n", result.num_noise);
        for(int idx : result.noise) {
            if(sketch_by_file) {
                fprintf(fp, "%s\n", sketches[idx].fileName.c_str());
            } else {
                fprintf(fp, "%s\n", sketches[idx].seqInfo.name.c_str());
            }
        }
    }
    
    fclose(fp);
}

/**
 * @brief Print DBSCAN results for KSSD sketches
 */
void printKssdDBSCANResult(
    const DBSCANResult& result,
    const vector<KssdSketchInfo>& sketches,
    bool sketch_by_file,
    const string& output_file,
    double eps,
    int minPts
) {
    FILE *fp = fopen(output_file.c_str(), "w");
    if(!fp) {
        cerr << "Error in printKssdDBSCANResult(), cannot open file: " << output_file << endl;
        exit(1);
    }
    
    fprintf(fp, "# DBSCAN Clustering Results\n");
    fprintf(fp, "# Parameters: eps=%.6f, minPts=%d\n", eps, minPts);
    fprintf(fp, "# Total clusters: %d\n", result.num_clusters);
    fprintf(fp, "# Total noise points (outliers): %d\n", result.num_noise);
    fprintf(fp, "#\n");
    
    // Print clusters
    for(size_t i = 0; i < result.clusters.size(); i++) {
        fprintf(fp, "# Cluster %zu (size: %zu)\n", i, result.clusters[i].size());
        for(int idx : result.clusters[i]) {
            if(sketch_by_file) {
                fprintf(fp, "%s\n", sketches[idx].fileName.c_str());
            } else {
                fprintf(fp, "%s\n", sketches[idx].seqInfo.name.c_str());
            }
        }
        fprintf(fp, "\n");
    }
    
    // Print noise points (outliers)
    if(!result.noise.empty()) {
        fprintf(fp, "# Noise points (outliers) - %d points\n", result.num_noise);
        for(int idx : result.noise) {
            if(sketch_by_file) {
                fprintf(fp, "%s\n", sketches[idx].fileName.c_str());
            } else {
                fprintf(fp, "%s\n", sketches[idx].seqInfo.name.c_str());
            }
        }
    }
    
    fclose(fp);
}

