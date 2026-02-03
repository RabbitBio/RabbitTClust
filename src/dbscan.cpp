#include "dbscan.h"
#include "phmap.h"   // For phmap::flat_hash_map
#include <iostream>
#include <iomanip>   // For std::setprecision, std::fixed
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <functional>
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

struct InvertedIndexCSR {
    vector<uint32_t> keys;           // Unique hash keys
    vector<uint64_t> offsets;        // Offsets into postings (size = keys + 1)
    vector<int> postings;            // Flat postings array
    phmap::flat_hash_map<uint32_t, uint32_t> key_to_idx;
};

/**
 * @brief Build inverted index in CSR format for all sketches (KSSD, 32-bit hashes)
 * Optimization 7: Skip empty sketches and ensure hash uniqueness
 */
InvertedIndexCSR buildInvertedIndexCSR32(
    const vector<KssdSketchInfo>& sketches,
    int max_posting
) {
    InvertedIndexCSR index;
    phmap::flat_hash_map<uint32_t, uint32_t> counts;
    counts.reserve(10000000);  // Pre-allocate for large datasets
    
    // First pass: count postings per hash
    for(size_t i = 0; i < sketches.size(); i++) {
        if(sketches[i].use64) continue;
        if(sketches[i].hash32_arr.empty()) continue;
        for(uint32_t hash : sketches[i].hash32_arr) {
            auto it = counts.find(hash);
            if(it == counts.end()) {
                counts.emplace(hash, 1u);
            } else {
                it->second++;
            }
        }
    }
    
    index.keys.reserve(counts.size());
    index.offsets.reserve(counts.size() + 1);
    index.key_to_idx.reserve(counts.size());
    
    uint64_t total = 0;
    for(const auto& kv : counts) {
        if(max_posting > 0 && (int)kv.second > max_posting) {
            continue;
        }
        index.key_to_idx.emplace(kv.first, (uint32_t)index.keys.size());
        index.keys.push_back(kv.first);
        index.offsets.push_back(total);
        total += kv.second;
    }
    index.offsets.push_back(total);
    index.postings.resize((size_t)total);
    
    // Second pass: fill postings
    vector<uint64_t> cursor(index.keys.size(), 0);
    for(size_t i = 0; i < sketches.size(); i++) {
        if(sketches[i].use64) continue;
        if(sketches[i].hash32_arr.empty()) continue;
        for(uint32_t hash : sketches[i].hash32_arr) {
            auto it = index.key_to_idx.find(hash);
            if(it == index.key_to_idx.end()) continue;
            uint32_t idx = it->second;
            uint64_t pos = index.offsets[idx] + cursor[idx]++;
            index.postings[(size_t)pos] = (int)i;
        }
    }
    
    return index;
}

inline bool getPostingList(
    const InvertedIndexCSR& index,
    uint32_t hash,
    uint64_t& start,
    uint64_t& end
) {
    auto it = index.key_to_idx.find(hash);
    if(it == index.key_to_idx.end()) return false;
    uint32_t idx = it->second;
    start = index.offsets[idx];
    end = index.offsets[idx + 1];
    return start < end;
}

/**
 * @brief Build k-NN list for a single point (approximate DBSCAN accelerator)
 * Find and store the k nearest neighbors with Jaccard scores.
 *
 * Note: This is an approximate accelerator for DBSCAN. If eps-neighborhood
 * size exceeds k, DBSCAN results may be approximate.
 * 
 * @param sketches Input sketches
 * @param inverted_index Inverted index for fast candidate lookup
 * @param jaccard_min Minimum Jaccard threshold for filtering candidates
 * @param point_idx Point index to build k-NN for
 * @param k Number of nearest neighbors to keep
 * @param max_posting Max posting list size to consider (0 = disabled)
 * @param knn_graph Output k-NN graph: point_id -> vector of (neighbor_id, jaccard)
 * @param knn_built Per-point cache flag
 * @param knn_touched_sizes Stats: touched sizes
 * @param knn_skipped_postings Stats: skipped postings
 * @param mark Scratch: epoch marks
 * @param cnt Scratch: intersection counts
 * @param touched Scratch: touched candidates
 * @param cur_mark Scratch: current epoch
 */
void buildKNNForPoint(
    const vector<KssdSketchInfo>& sketches,
    const InvertedIndexCSR& inverted_index,
    double jaccard_min,
    int point_idx,
    int k,
    int max_posting,
    vector<vector<std::pair<int, float>>>& knn_graph,
    vector<char>& knn_built,
    vector<int>& knn_touched_sizes,
    uint64_t& knn_skipped_postings,
    vector<uint32_t>& mark,
    vector<uint32_t>& cnt,
    vector<int>& touched,
    uint32_t& cur_mark
) {
    int n = sketches.size();
    if(point_idx < 0 || point_idx >= n) return;
    if(knn_built[point_idx]) return;
    if(k <= 0) {
        knn_built[point_idx] = 1;
        return;
    }
    knn_built[point_idx] = 1;
    
    if(sketches[point_idx].use64 || sketches[point_idx].hash32_arr.empty()) {
        return;
    }
    
    // Find all candidates using inverted index
    touched.clear();
    cur_mark++;
    if(cur_mark == 0) {
        std::fill(mark.begin(), mark.end(), 0);
        cur_mark = 1;
    }
    
    int sizeRef = sketches[point_idx].hash32_arr.size();
    
    // Count intersections
    for(uint32_t hash : sketches[point_idx].hash32_arr) {
        uint64_t start = 0, end = 0;
        if(!getPostingList(inverted_index, hash, start, end)) continue;
        uint64_t len = end - start;
        if(max_posting > 0 && (int)len > max_posting) {
            knn_skipped_postings += 1;
            continue;
        }
        for(uint64_t p = start; p < end; p++) {
            int candidate_id = inverted_index.postings[(size_t)p];
            if(candidate_id == point_idx) continue;
            if(mark[candidate_id] != cur_mark) {
                mark[candidate_id] = cur_mark;
                cnt[candidate_id] = 1;
                touched.push_back(candidate_id);
            } else {
                cnt[candidate_id]++;
            }
        }
    }
    
    knn_touched_sizes.push_back((int)touched.size());
    if(touched.empty()) {
        knn_built[point_idx] = 1;
        return;
    }
    
    // Use min-heap to maintain top-k highest Jaccard
    using ScorePair = std::pair<float, int>;
    std::priority_queue<ScorePair, std::vector<ScorePair>, std::greater<ScorePair>> topk;  // min-heap
    
    // Serial processing with strong pruning
    for(int candidate_id : touched) {
        const auto& candidate = sketches[candidate_id];
        if(candidate.use64) continue;
        
        int common = cnt[candidate_id];
        int sizeQry = candidate.hash32_arr.size();
        
        // Filter by common
        double v = jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min);
        int min_common_needed = (int)std::ceil(v - 1e-12);
        if(common < min_common_needed) continue;
        
        // Calculate Jaccard
        uint64_t union_size = sizeRef + sizeQry - common;
        double jaccard_ = (union_size == 0) ? 0.0 : (double)common / union_size;
        if(jaccard_ < jaccard_min) continue;
        
        float score = (float)jaccard_;
        if((int)topk.size() < k) {
            topk.push({score, candidate_id});
        } else if(score > topk.top().first) {
            topk.pop();
            topk.push({score, candidate_id});
        }
    }
    
    // Extract k-NN neighbors
    knn_graph[point_idx].reserve(k);
    while(!topk.empty()) {
        knn_graph[point_idx].push_back({topk.top().second, topk.top().first});
        topk.pop();
    }
    
    knn_built[point_idx] = 1;
}

/**
 * @brief Find all points within epsilon distance using inverted index (KSSD, 32-bit)
 * Optimized version using global stamp+count arrays and direct Jaccard calculation
 */
vector<int> findNeighborsKSSDWithIndex(
    const vector<KssdSketchInfo>& sketches,
    int point_idx,
    double jaccard_min,           // Pre-calculated jaccard_min (optimization 5)
    double eps,                   // Epsilon threshold for DBSCAN
    const InvertedIndexCSR& inverted_index,
    const vector<vector<std::pair<int, float>>>& knn_graph,  // k-NN graph with jaccard
    int max_posting,              // Max posting list size to consider (0 = disabled)
    vector<uint32_t>& mark,      // Global mark array (epoch-based)
    vector<uint32_t>& cnt,        // Global count array
    vector<int>& touched,         // Global touched list
    uint32_t& cur_mark,          // Current epoch marker
    vector<vector<int>>& thread_neighbors,  // Optimization 2: Reuse thread-local buffers
    int threads
) {
    const auto& point = sketches[point_idx];
    
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
    
    // If k-NN graph is available, use it directly (approximate DBSCAN)
    if(!knn_graph.empty()) {
        vector<int> neighbors;
        neighbors.reserve(knn_graph[point_idx].size());
        for(const auto& entry : knn_graph[point_idx]) {
            if(entry.second >= jaccard_min) {
                neighbors.push_back(entry.first);
            }
        }
        return neighbors;
    }
    
    // Use inverted index for 32-bit hashes with optimized stamp+count arrays
    int sizeRef = point.hash32_arr.size();
    if(sizeRef == 0) return vector<int>();  // Optimization 7: Skip empty sketches
    
    // Standard mode: use inverted index to find all candidates
    touched.clear();
    
    // Step 1: Count intersections using inverted index (no hash map allocation)
    cur_mark++;
    if(cur_mark == 0) {
        std::fill(mark.begin(), mark.end(), 0);
        cur_mark = 1;
    }
    
    for(uint32_t hash : point.hash32_arr) {
        uint64_t start = 0, end = 0;
        if(!getPostingList(inverted_index, hash, start, end)) continue;
        uint64_t len = end - start;
        if(max_posting > 0 && (int)len > max_posting) continue;
        for(uint64_t p = start; p < end; p++) {
            int candidate_id = inverted_index.postings[(size_t)p];
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
    int threads,
    int knn_k,
    int max_posting
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
    
    if(knn_k > 0 && knn_k < (minPts - 1)) {
        cerr << "-----WARNING: knn_k (" << knn_k << ") < minPts-1 (" << (minPts - 1)
             << "). Adjusting knn_k to " << (minPts - 1) << "." << endl;
        knn_k = minPts - 1;
    } else if(knn_k > 0 && knn_k < 5 * (minPts - 1)) {
        cerr << "-----WARNING: knn_k (" << knn_k << ") may be too small for stable DBSCAN."
             << " Consider knn_k >= " << (5 * (minPts - 1)) << "." << endl;
    }
    
    // Build inverted index for acceleration (only for 32-bit hashes)
    cerr << "-----Building inverted index for acceleration..." << endl;
    InvertedIndexCSR inverted_index = buildInvertedIndexCSR32(sketches, max_posting);
    cerr << "-----Inverted index built: " << inverted_index.keys.size() << " unique hashes" << endl;
    
    // Build k-NN graph lazily for approximate acceleration (if enabled)
    vector<vector<std::pair<int, float>>> knn_graph;
    vector<char> knn_built;
    vector<int> knn_touched_sizes;
    uint64_t knn_skipped_postings = 0;
    if(knn_k > 0) {
        cerr << "-----WARNING: k-NN acceleration is approximate for DBSCAN (may miss eps neighbors if k is small)." << endl;
        knn_graph.resize(sketches.size());
        knn_built.assign(sketches.size(), 0);
        knn_touched_sizes.reserve(sketches.size() / 4 + 1);
    }
    
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
    
    // Statistics for progress reporting
    int processed_count = 0;
    int noise_count = 0;
    int last_reported = 0;
    const int REPORT_INTERVAL = std::max(1000, n / 100);  // Report every 1% or 1000 points, whichever is larger
    
    // RegionQuery stats
    uint64_t total_region_queries = 0;
    uint64_t total_neighbors = 0;
    uint64_t core_points = 0;
    
    // Process each point
    for(int i = 0; i < n; i++) {
        if(labels[i] != -1) continue;  // Already processed
        
        processed_count++;
        
        // Find neighbors using inverted index with optimized arrays (and k-NN graph if available)
        if(knn_k > 0 && !knn_built[i]) {
            buildKNNForPoint(sketches, inverted_index, jaccard_min, i, knn_k, max_posting,
                             knn_graph, knn_built, knn_touched_sizes, knn_skipped_postings,
                             mark, cnt, touched, cur_mark);
        }
        vector<int> neighbors = findNeighborsKSSDWithIndex(
            sketches, i, jaccard_min, eps, inverted_index, knn_graph, max_posting,
            mark, cnt, touched, cur_mark, thread_neighbors, threads
        );
        total_region_queries++;
        total_neighbors += neighbors.size();
        
        // Optimization 4: minPts includes the point itself
        if((int)neighbors.size() + 1 < minPts) {
            // Not a core point, mark as noise (may be reassigned later)
            labels[i] = -2;
            noise_count++;
            continue;
        }
        core_points++;
        
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
        int cluster_size = 1;  // Count cluster size for progress
        while(!seed_set.empty()) {
            int q = seed_set.front();
            seed_set.pop();
            
            if(labels[q] == -2) {
                // Reassign noise point to cluster
                labels[q] = cluster_id;
                noise_count--;  // Adjust noise count
                cluster_size++;
            } else if(labels[q] != -1) {
                // Already processed
                continue;
            }
            
            labels[q] = cluster_id;
            cluster_size++;
            processed_count++;
            
            // Find neighbors of q using inverted index (and k-NN graph if available)
            if(knn_k > 0 && !knn_built[q]) {
                buildKNNForPoint(sketches, inverted_index, jaccard_min, q, knn_k, max_posting,
                                 knn_graph, knn_built, knn_touched_sizes, knn_skipped_postings,
                                 mark, cnt, touched, cur_mark);
            }
            vector<int> q_neighbors = findNeighborsKSSDWithIndex(
                sketches, q, jaccard_min, eps, inverted_index, knn_graph, max_posting,
                mark, cnt, touched, cur_mark, thread_neighbors, threads
            );
            total_region_queries++;
            total_neighbors += q_neighbors.size();
            
            // Optimization 4: minPts includes the point itself
            if((int)q_neighbors.size() + 1 >= minPts) {
                // q is also a core point, add its neighbors to seed set (with deduplication)
                core_points++;
                for(int neighbor : q_neighbors) {
                    if((labels[neighbor] == -1 || labels[neighbor] == -2) && !inq(neighbor)) {
                        seed_set.push(neighbor);
                        set_inq(neighbor);
                    }
                }
            }
        }
        
        cluster_id++;
        
        // Progress report (similar to greedy style)
        if(processed_count - last_reported >= REPORT_INTERVAL || processed_count == n) {
            double progress_pct = 100.0 * processed_count / n;
            double clustering_rate = (n > 0) ? 100.0 * (n - noise_count) / n : 0.0;
            
            cerr << "Progress: " << processed_count << "/" << n 
                 << " (" << std::fixed << std::setprecision(1) << progress_pct << "%)"
                 << " | Clusters: " << cluster_id
                 << " | Noise: " << noise_count
                 << " | Clustered: " << std::fixed << std::setprecision(1) << clustering_rate << "%"
                 << endl;
            
            last_reported = processed_count;
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
    if(total_region_queries > 0) {
        double avg_neighbors = (double)total_neighbors / (double)total_region_queries;
        cerr << "-----RegionQuery stats: queries=" << total_region_queries
             << ", avg_neighbors=" << avg_neighbors << endl;
    }
    if(core_points > 0) {
        cerr << "-----Core points: " << core_points << " (" 
             << (100.0 * core_points / n) << "%)" << endl;
    }
    if(knn_k > 0 && !knn_touched_sizes.empty()) {
        std::vector<int> tmp = knn_touched_sizes;
        std::sort(tmp.begin(), tmp.end());
        size_t idx95 = (size_t)std::ceil(0.95 * tmp.size()) - 1;
        if(idx95 >= tmp.size()) idx95 = tmp.size() - 1;
        double avg_touched = 0.0;
        for(int v : knn_touched_sizes) avg_touched += v;
        avg_touched /= knn_touched_sizes.size();
        cerr << "-----k-NN stats: built=" << knn_touched_sizes.size()
             << ", avg_touched=" << avg_touched
             << ", p95_touched=" << tmp[idx95]
             << ", skipped_postings=" << knn_skipped_postings << endl;
    }
    
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
 * Format matches printResult() from MST_IO.cpp
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
    
    // Write threshold information at the beginning (matching printResult format)
    fprintf(fp, "# DBSCAN clustering parameters: eps=%.6f, minPts=%d\n", eps, minPts);
    fprintf(fp, "# Total clusters: %d\n", result.num_clusters);
    if(result.num_noise > 0) {
        fprintf(fp, "# Total noise points (outliers): %d\n", result.num_noise);
    }
    fprintf(fp, "#\n");
    
    if(sketch_by_file) {
        // Print clusters
        for(int i = 0; i < (int)result.clusters.size(); i++) {
            fprintf(fp, "the cluster %d is: \n", i);
            for(int j = 0; j < (int)result.clusters[i].size(); j++) {
                int curId = result.clusters[i][j];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printDBSCANResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                const char* seqName = "N/A";
                const char* seqComment = "N/A";
                if(!sketches[curId].fileSeqs.empty()) {
                    seqName = sketches[curId].fileSeqs[0].name.c_str();
                    seqComment = sketches[curId].fileSeqs[0].comment.c_str();
                }
                fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n",
                    j, curId, sketches[curId].totalSeqLength,
                    sketches[curId].fileName.c_str(), seqName, seqComment);
            }
            fprintf(fp, "\n");
        }
        
        // Print noise points as separate clusters (if any)
        if(!result.noise.empty()) {
            for(int i = 0; i < (int)result.noise.size(); i++) {
                int curId = result.noise[i];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printDBSCANResult(), invalid noise index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "the cluster %d is: \n", (int)result.clusters.size() + i);
                const char* seqName = "N/A";
                const char* seqComment = "N/A";
                if(!sketches[curId].fileSeqs.empty()) {
                    seqName = sketches[curId].fileSeqs[0].name.c_str();
                    seqComment = sketches[curId].fileSeqs[0].comment.c_str();
                }
                fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n",
                    0, curId, sketches[curId].totalSeqLength,
                    sketches[curId].fileName.c_str(), seqName, seqComment);
                fprintf(fp, "\n");
            }
        }
    } else {
        // sketch by sequence
        // Print clusters
        for(int i = 0; i < (int)result.clusters.size(); i++) {
            fprintf(fp, "the cluster %d is: \n", i);
            for(int j = 0; j < (int)result.clusters[i].size(); j++) {
                int curId = result.clusters[i][j];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printDBSCANResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n",
                    j, curId, sketches[curId].seqInfo.length,
                    sketches[curId].seqInfo.name.c_str(),
                    sketches[curId].seqInfo.comment.c_str());
            }
            fprintf(fp, "\n");
        }
        
        // Print noise points as separate clusters (if any)
        if(!result.noise.empty()) {
            for(int i = 0; i < (int)result.noise.size(); i++) {
                int curId = result.noise[i];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printDBSCANResult(), invalid noise index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "the cluster %d is: \n", (int)result.clusters.size() + i);
                fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n",
                    0, curId, sketches[curId].seqInfo.length,
                    sketches[curId].seqInfo.name.c_str(),
                    sketches[curId].seqInfo.comment.c_str());
                fprintf(fp, "\n");
            }
        }
    }
    
    fclose(fp);
}

/**
 * @brief Print DBSCAN results for KSSD sketches
 * Format matches printKssdResult() from MST_IO.cpp
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
    
    // Write threshold information at the beginning (matching printKssdResult format)
    fprintf(fp, "# DBSCAN clustering parameters: eps=%.6f, minPts=%d\n", eps, minPts);
    fprintf(fp, "# Total clusters: %d\n", result.num_clusters);
    if(result.num_noise > 0) {
        fprintf(fp, "# Total noise points (outliers): %d\n", result.num_noise);
    }
    fprintf(fp, "#\n");
    
    if(sketch_by_file) {
        // Print clusters
        for(int i = 0; i < (int)result.clusters.size(); i++) {
            fprintf(fp, "the cluster %d is: \n", i);
            for(int j = 0; j < (int)result.clusters[i].size(); j++) {
                int curId = result.clusters[i][j];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printKssdDBSCANResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                const char* seqName = "N/A";
                const char* seqComment = "N/A";
                if(!sketches[curId].fileSeqs.empty()) {
                    seqName = sketches[curId].fileSeqs[0].name.c_str();
                    seqComment = sketches[curId].fileSeqs[0].comment.c_str();
                }
                fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n",
                    j, curId, sketches[curId].totalSeqLength,
                    sketches[curId].fileName.c_str(), seqName, seqComment);
            }
            fprintf(fp, "\n");
        }
        
        // Print noise points as separate clusters (if any)
        if(!result.noise.empty()) {
            for(int i = 0; i < (int)result.noise.size(); i++) {
                int curId = result.noise[i];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printKssdDBSCANResult(), invalid noise index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "the cluster %d is: \n", (int)result.clusters.size() + i);
                const char* seqName = "N/A";
                const char* seqComment = "N/A";
                if(!sketches[curId].fileSeqs.empty()) {
                    seqName = sketches[curId].fileSeqs[0].name.c_str();
                    seqComment = sketches[curId].fileSeqs[0].comment.c_str();
                }
                fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n",
                    0, curId, sketches[curId].totalSeqLength,
                    sketches[curId].fileName.c_str(), seqName, seqComment);
                fprintf(fp, "\n");
            }
        }
    } else {
        // sketch by sequence
        // Print clusters
        for(int i = 0; i < (int)result.clusters.size(); i++) {
            fprintf(fp, "the cluster %d is: \n", i);
            for(int j = 0; j < (int)result.clusters[i].size(); j++) {
                int curId = result.clusters[i][j];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printKssdDBSCANResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n",
                    j, curId, sketches[curId].seqInfo.length,
                    sketches[curId].seqInfo.name.c_str(),
                    sketches[curId].seqInfo.comment.c_str());
            }
            fprintf(fp, "\n");
        }
        
        // Print noise points as separate clusters (if any)
        if(!result.noise.empty()) {
            for(int i = 0; i < (int)result.noise.size(); i++) {
                int curId = result.noise[i];
                if(curId < 0 || curId >= (int)sketches.size()) {
                    cerr << "ERROR: printKssdDBSCANResult(), invalid noise index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
                    continue;
                }
                fprintf(fp, "the cluster %d is: \n", (int)result.clusters.size() + i);
                fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n",
                    0, curId, sketches[curId].seqInfo.length,
                    sketches[curId].seqInfo.name.c_str(),
                    sketches[curId].seqInfo.comment.c_str());
                fprintf(fp, "\n");
            }
        }
    }
    
    fclose(fp);
}

