#ifdef LEIDEN_CLUST

#include "leiden.h"
#include "SketchInfo.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;

// ==================== Dynamic Inverted Index (same as greedy for now) ====================
class DynamicInvertedIndex {
private:
    unordered_map<uint32_t, vector<uint32_t>> index_map;
    unordered_set<int> representative_set;

public:
    void add_representative(int rep_id, const vector<uint32_t>& hash_array) {
        representative_set.insert(rep_id);
        for (uint32_t hash : hash_array) {
            index_map[hash].push_back(rep_id);
        }
    }

    bool is_representative(int id) const {
        return representative_set.find(id) != representative_set.end();
    }

    void calculate_intersections_sparse(const vector<uint32_t>& hash_array, 
                                       unordered_map<int, int>& intersection_map) const {
        intersection_map.clear();
        for (uint32_t hash : hash_array) {
            auto it = index_map.find(hash);
            if (it != index_map.end()) {
                for (uint32_t rep_id : it->second) {
                    intersection_map[rep_id]++;
                }
            }
        }
    }

    size_t num_representatives() const {
        return representative_set.size();
    }

    void clear() {
        index_map.clear();
        representative_set.clear();
    }
};

// ==================== Distance Calculation ====================
static inline double calculate_mash_distance_fast(int common, int size1, int size2, int k) {
    if (common == 0) return 1.0;
    
    int union_size = size1 + size2 - common;
    if (union_size == 0) return 1.0;
    
    double jaccard = static_cast<double>(common) / union_size;
    if (jaccard <= 0.0) return 1.0;
    if (jaccard >= 1.0) return 0.0;
    
    double mash_dist = -1.0 / k * log(2.0 * jaccard / (1.0 + jaccard));
    return max(0.0, min(1.0, mash_dist));
}

// ==================== Leiden Clustering (Placeholder) ====================
// Currently implements the same algorithm as greedy with inverted index
// Will be replaced with true Leiden algorithm later
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    int threads,
    int kmer_size)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) {
        return vector<vector<int>>();
    }

    cerr << "-----Starting Leiden clustering (currently using greedy inverted index as placeholder)" << endl;
    cerr << "-----Number of genomes: " << numGenomes << endl;
    cerr << "-----Threads: " << threads << endl;
    cerr << "-----Note: True Leiden algorithm with modularity optimization will be implemented in future versions" << endl;

    // Initialize dynamic inverted index
    DynamicInvertedIndex dynamic_index;
    
    // Cluster assignment: -1 means not assigned, >=0 is cluster id
    vector<int> cluster_assignment(numGenomes, -1);
    
    // First genome is always the first cluster representative
    cluster_assignment[0] = 0;
    dynamic_index.add_representative(0, sketches[0].hash32_arr);
    
    // Temporary threshold for placeholder (will be removed in true Leiden)
    // For now, we use a heuristic: no threshold, just find closest representative
    
    // Process each genome sequentially
    for (int j = 1; j < numGenomes; j++) {
        unordered_map<int, int> intersection_map;
        
        // Query inverted index
        dynamic_index.calculate_intersections_sparse(sketches[j].hash32_arr, intersection_map);
        
        if (intersection_map.empty()) {
            // No intersection with any representative, create new cluster
            int new_cluster_id = dynamic_index.num_representatives();
            cluster_assignment[j] = new_cluster_id;
            dynamic_index.add_representative(j, sketches[j].hash32_arr);
            continue;
        }
        
        // Prepare candidates
        vector<pair<int, int>> candidates;
        candidates.reserve(intersection_map.size());
        for (const auto& entry : intersection_map) {
            candidates.push_back(entry);
        }
        
        // Parallel distance calculation
        vector<int> best_rep_per_thread(threads, -1);
        vector<double> best_dist_per_thread(threads, 1.0);
        
        int size_j = sketches[j].sketchsize;
        
        #pragma omp parallel num_threads(threads)
        {
            int tid = omp_get_thread_num();
            int local_best_rep = -1;
            double local_best_dist = 1.0;
            
            #pragma omp for schedule(dynamic)
            for (size_t k = 0; k < candidates.size(); k++) {
                int repId = candidates[k].first;
                int common = candidates[k].second;
                
                int size_rep = sketches[repId].sketchsize;
                double dist = calculate_mash_distance_fast(common, size_j, size_rep, kmer_size);
                
                if (dist < local_best_dist) {
                    local_best_dist = dist;
                    local_best_rep = repId;
                }
            }
            
            best_rep_per_thread[tid] = local_best_rep;
            best_dist_per_thread[tid] = local_best_dist;
        }
        
        // Merge results
        int best_rep = -1;
        double best_dist = 1.0;
        for (int t = 0; t < threads; t++) {
            if (best_dist_per_thread[t] < best_dist) {
                best_dist = best_dist_per_thread[t];
                best_rep = best_rep_per_thread[t];
            }
        }
        
        if (best_rep != -1) {
            // Assign to existing cluster
            cluster_assignment[j] = cluster_assignment[best_rep];
        } else {
            // Create new cluster
            int new_cluster_id = dynamic_index.num_representatives();
            cluster_assignment[j] = new_cluster_id;
            dynamic_index.add_representative(j, sketches[j].hash32_arr);
        }
        
        // Progress reporting
        if ((j + 1) % 5000 == 0 || j + 1 == numGenomes) {
            cerr << "Progress: " << (j + 1) << "/" << numGenomes 
                 << " | Clusters: " << dynamic_index.num_representatives() << endl;
        }
    }
    
    // Convert cluster_assignment to cluster format
    unordered_map<int, vector<int>> cluster_map;
    for (int i = 0; i < numGenomes; i++) {
        cluster_map[cluster_assignment[i]].push_back(i);
    }
    
    vector<vector<int>> result;
    result.reserve(cluster_map.size());
    for (const auto& entry : cluster_map) {
        result.push_back(entry.second);
    }
    
    cerr << "-----Leiden clustering complete" << endl;
    cerr << "-----Total clusters: " << result.size() << endl;
    
    return result;
}

#endif // LEIDEN_CLUST

