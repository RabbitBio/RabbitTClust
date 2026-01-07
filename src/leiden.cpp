#ifdef LEIDEN_CLUST

#include "leiden.h"
#include "SketchInfo.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <igraph.h>

using namespace std;

// ==================== Global Inverted Index for Graph Construction ====================
class GlobalInvertedIndex {
private:
    unordered_map<uint32_t, vector<int>> index_map;
    int num_sequences;

public:
    GlobalInvertedIndex() : num_sequences(0) {}
    
    void build(const vector<KssdSketchInfo>& sketches, int threads) {
        num_sequences = sketches.size();
        
        cerr << "-----Building global inverted index..." << endl;
        
        vector<unordered_map<uint32_t, vector<int>>> thread_local_indices(threads);
        
        #pragma omp parallel num_threads(threads)
        {
            int tid = omp_get_thread_num();
            auto& local_index = thread_local_indices[tid];
            
            #pragma omp for schedule(dynamic, 100)
            for (int i = 0; i < num_sequences; i++) {
                const auto& sketch = sketches[i];
                for (uint32_t hash : sketch.hash32_arr) {
                    local_index[hash].push_back(i);
                }
            }
        }
        
        cerr << "-----Merging thread-local indices..." << endl;
        for (int tid = 0; tid < threads; tid++) {
            for (auto& entry : thread_local_indices[tid]) {
                uint32_t hash = entry.first;
                auto& seq_ids = entry.second;
                index_map[hash].insert(index_map[hash].end(), seq_ids.begin(), seq_ids.end());
            }
        }
        
        cerr << "-----Inverted index built: " << index_map.size() << " unique hashes" << endl;
    }
    
    const vector<int>* get_sequences_with_hash(uint32_t hash) const {
        auto it = index_map.find(hash);
        if (it != index_map.end()) {
            return &(it->second);
        }
        return nullptr;
    }
    
    void calculate_all_intersections(int seq_id, const vector<uint32_t>& hash_array,
                                     unordered_map<int, int>& intersection_map) const {
        intersection_map.clear();
        for (uint32_t hash : hash_array) {
            auto* seq_list = get_sequences_with_hash(hash);
            if (seq_list) {
                for (int other_id : *seq_list) {
                    if (other_id != seq_id) {
                        intersection_map[other_id]++;
                    }
                }
            }
        }
    }
    
    void clear() {
        index_map.clear();
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

// ==================== Edge Structure ====================
struct Edge {
    int from;
    int to;
    double weight;
    
    Edge(int f, int t, double w) : from(f), to(t), weight(w) {}
};

// ==================== Leiden Clustering using igraph ====================
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution,
    bool use_modularity)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) {
        return vector<vector<int>>();
    }

    cerr << "=============================" << endl;
    cerr << "Graph-based Clustering (Louvain)" << endl;
    cerr << "=============================" << endl;
    cerr << "Genomes: " << numGenomes << endl;
    cerr << "Edge threshold: " << threshold << endl;
    cerr << "Resolution: " << resolution << endl;
    cerr << "Threads: " << threads << endl;
    cerr << "=============================" << endl;

    // Step 1: Build inverted index
    GlobalInvertedIndex global_index;
    global_index.build(sketches, threads);

    // Step 2: Build edge list using inverted index
    cerr << "-----Building similarity graph..." << endl;
    
    vector<vector<Edge>> thread_local_edges(threads);
    long long total_comparisons = 0;
    long long edges_created = 0;
    
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        auto& local_edges = thread_local_edges[tid];
        unordered_map<int, int> intersection_map;
        
        long long local_comparisons = 0;
        long long local_edge_count = 0;
        
        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < numGenomes; i++) {
            const auto& sketch_i = sketches[i];
            int size_i = sketch_i.sketchsize;
            
            global_index.calculate_all_intersections(i, sketch_i.hash32_arr, intersection_map);
            
            for (const auto& entry : intersection_map) {
                int j = entry.first;
                
                if (i >= j) continue;
                
                int common = entry.second;
                int size_j = sketches[j].sketchsize;
                
                double size_ratio = (size_i < size_j) ? 
                                   (double)size_i / size_j : 
                                   (double)size_j / size_i;
                if (size_ratio < 0.5) continue;
                
                double dist = calculate_mash_distance_fast(common, size_i, size_j, kmer_size);
                local_comparisons++;
                
                if (dist < threshold) {
                    double weight = 1.0 - dist;
                    local_edges.push_back(Edge(i, j, weight));
                    local_edge_count++;
                }
            }
            
            if ((i + 1) % 1000 == 0) {
                #pragma omp critical
                {
                    cerr << "Progress: " << (i + 1) << "/" << numGenomes << endl;
                }
            }
        }
        
        #pragma omp critical
        {
            total_comparisons += local_comparisons;
            edges_created += local_edge_count;
        }
    }
    
    cerr << "-----Comparisons: " << total_comparisons << endl;
    cerr << "-----Edges: " << edges_created << endl;
    
    // Merge edges
    vector<Edge> all_edges;
    for (const auto& local_edges : thread_local_edges) {
        all_edges.insert(all_edges.end(), local_edges.begin(), local_edges.end());
    }
    
    if (all_edges.empty()) {
        cerr << "-----Warning: No edges! Each sequence in its own cluster." << endl;
        vector<vector<int>> result(numGenomes);
        for (int i = 0; i < numGenomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    // Step 3: Build igraph
    cerr << "-----Building igraph (parallel)..." << endl;
    
    igraph_t graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    igraph_vector_int_init(&edges_vec, all_edges.size() * 2);
    igraph_vector_init(&weights_vec, all_edges.size());
    
    // Parallel edge vector construction
    size_t num_edges = all_edges.size();
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < num_edges; i++) {
        VECTOR(edges_vec)[2 * i] = all_edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = all_edges[i].to;
        VECTOR(weights_vec)[i] = all_edges[i].weight;
    }
    
    igraph_create(&graph, &edges_vec, numGenomes, IGRAPH_UNDIRECTED);
    
    // Step 4: Run community detection algorithm
    cerr << "-----Running community detection algorithm..." << endl;
    
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, numGenomes);
    igraph_vector_init(&modularity_vec, 0);
    
    // Use Multilevel (Louvain) algorithm - more stable than Leiden for weighted graphs
    // This algorithm is well-tested and widely used
    int result_code = igraph_community_multilevel(
        &graph,
        &weights_vec,     // edge weights
        resolution,       // resolution parameter  
        &membership,
        NULL,             // memberships (intermediate steps, not needed)
        &modularity_vec   // final modularity (vector)
    );
    
    // Get final modularity value
    double modularity = 0.0;
    if (igraph_vector_size(&modularity_vec) > 0) {
        modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
    }
    
    cerr << "-----Community detection complete!" << endl;
    cerr << "-----Final modularity: " << modularity << endl;
    
    // Count clusters
    int nb_clusters = 0;
    for (int i = 0; i < numGenomes; i++) {
        if (VECTOR(membership)[i] + 1 > nb_clusters) {
            nb_clusters = VECTOR(membership)[i] + 1;
        }
    }
    
    if (result_code != IGRAPH_SUCCESS) {
        cerr << "-----ERROR: Community detection failed with code " << result_code << endl;
        igraph_vector_destroy(&modularity_vec);
        igraph_vector_int_destroy(&membership);
        igraph_vector_destroy(&weights_vec);
        igraph_vector_int_destroy(&edges_vec);
        igraph_destroy(&graph);
        
        // Return each sequence in its own cluster as fallback
        vector<vector<int>> result(numGenomes);
        for (int i = 0; i < numGenomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    cerr << "-----Number of clusters: " << nb_clusters << endl;
    
    // Step 5: Extract clusters
    unordered_map<int, vector<int>> cluster_map;
    for (int i = 0; i < numGenomes; i++) {
        int community = VECTOR(membership)[i];
        cluster_map[community].push_back(i);
    }
    
    vector<vector<int>> result;
    result.reserve(cluster_map.size());
    for (const auto& entry : cluster_map) {
        result.push_back(entry.second);
    }
    
    sort(result.begin(), result.end(), 
         [](const vector<int>& a, const vector<int>& b) {
             return a.size() > b.size();
         });
    
    cerr << "=============================" << endl;
    cerr << "Clustering Complete" << endl;
    cerr << "Total clusters: " << result.size() << endl;
    cerr << "Largest cluster: " << (result.empty() ? 0 : result[0].size()) << endl;
    cerr << "Average size: " << (result.empty() ? 0.0 : (double)numGenomes / result.size()) << endl;
    cerr << "Final modularity: " << modularity << endl;
    cerr << "=============================" << endl;
    
    // Cleanup
    igraph_vector_destroy(&modularity_vec);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_int_destroy(&edges_vec);
    igraph_destroy(&graph);
    
    return result;
}

#endif // LEIDEN_CLUST
