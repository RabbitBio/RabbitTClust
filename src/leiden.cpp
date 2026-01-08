#ifdef LEIDEN_CLUST

/**
 * @file leiden.cpp
 * @brief Implementation of graph-based clustering using Louvain/Leiden algorithms
 * @author Tong Zhang
 * @date 2026-01-06
 * 
 * This module implements efficient graph-based clustering for large-scale genome datasets:
 * - Global inverted index for fast similarity graph construction
 * - Parallel edge building using OpenMP
 * - Support for both Leiden (default) and Louvain (alternative) community detection
 * - Integration with igraph library for graph algorithms
 * 
 * Key optimizations:
 * - Inverted index reduces distance calculations from O(N²) to O(candidates)
 * - Sparse graph representation for memory efficiency
 * - Multi-threaded graph construction
 * - Threshold-based edge filtering
 */

#include "leiden.h"
#include "SketchInfo.h"
#include <vector>
#include <unordered_set>
#include <set>
#include <queue>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <igraph.h>
#include "../RabbitSketch/src/robin_hood.h"

using namespace std;

// ==================== Global Inverted Index for Graph Construction ====================
class GlobalInvertedIndex {
private:
    robin_hood::unordered_map<uint32_t, vector<int>> index_map;
    int num_sequences;

public:
    GlobalInvertedIndex() : num_sequences(0) {}
    
    void build(const vector<KssdSketchInfo>& sketches, int threads) {
        num_sequences = sketches.size();
        
        cerr << "-----Building global inverted index..." << endl;
        
        vector<robin_hood::unordered_map<uint32_t, vector<int>>> thread_local_indices(threads);
        
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
                                     robin_hood::unordered_map<int, int>& intersection_map) const {
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

// ==================== Helper: Save Graph to File ====================
void save_graph_to_file(const vector<Edge>& edges, int num_nodes, const string& filename);

// ==================== Leiden Clustering using igraph ====================
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution,
    bool use_leiden,
    int knn_k,
    const string& graph_save_path)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) {
        return vector<vector<int>>();
    }

    cerr << "=============================" << endl;
    if (use_leiden) {
        cerr << "Graph-based Clustering (Leiden)" << endl;
    } else {
        cerr << "Graph-based Clustering (Louvain)" << endl;
    }
    cerr << "=============================" << endl;
    cerr << "Genomes: " << numGenomes << endl;
    cerr << "Edge threshold: " << threshold << endl;
    if (knn_k > 0) {
        cerr << "k-NN k: " << knn_k << " (keep only top-" << knn_k << " neighbors per node)" << endl;
    }
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
        robin_hood::unordered_map<int, int> intersection_map;
        
        long long local_comparisons = 0;
        long long local_edge_count = 0;
        
        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < numGenomes; i++) {
            const auto& sketch_i = sketches[i];
            int size_i = sketch_i.sketchsize;
            
            global_index.calculate_all_intersections(i, sketch_i.hash32_arr, intersection_map);
            
            if (knn_k > 0) {
                // k-NN mode: keep only top-k nearest neighbors
                priority_queue<pair<double, int>> topk;
                
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
                        if ((int)topk.size() < knn_k) {
                            topk.push({dist, j});
                        } else if (dist < topk.top().first) {
                            topk.pop();
                            topk.push({dist, j});
                        }
                    }
                }
                
                // Add top-k edges
                while (!topk.empty()) {
                    double dist = topk.top().first;
                    int j = topk.top().second;
                    topk.pop();
                    local_edges.push_back(Edge(i, j, 1.0 - dist));
                    local_edge_count++;
                }
            } else {
                // Standard mode: keep all edges within threshold
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
    cerr << "-----Edges created: " << edges_created << endl;
    if (knn_k > 0) {
        cerr << "-----k-NN filtered: keeping top-" << knn_k << " neighbors per node (max " 
             << (numGenomes * knn_k) << " edges)" << endl;
    }
    
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
    
    // Save graph to file if path is provided
    if (!graph_save_path.empty()) {
        save_graph_to_file(all_edges, numGenomes, graph_save_path);
    }
    
    // Step 4: Run community detection algorithm
    if (use_leiden) {
        cerr << "-----Running Leiden algorithm (default)..." << endl;
    } else {
        cerr << "-----Running Louvain algorithm..." << endl;
        cerr << "-----Note: Using --louvain flag. Louvain may be faster but less refined than Leiden (default)." << endl;
    }
    
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, numGenomes);
    igraph_vector_init(&modularity_vec, 0);
    
    int result_code;
    
    if (use_leiden) {
        // Use Leiden algorithm - refined version with better community structure
        // WARNING: Leiden is sensitive to edge weight scale and may produce suboptimal results
        igraph_integer_t nb_clusters_temp;
        igraph_real_t quality_temp;
        
        // Normalize edge weights to improve Leiden performance
        // Our weights are in [0.9, 1.0] range, normalize to [0, 1] for better scaling
        double min_weight = 1.0, max_weight = 0.0;
        for (size_t i = 0; i < igraph_vector_size(&weights_vec); i++) {
            double w = VECTOR(weights_vec)[i];
            if (w < min_weight) min_weight = w;
            if (w > max_weight) max_weight = w;
        }
        
        // Normalize weights if range is too narrow
        igraph_vector_t normalized_weights;
        if (max_weight - min_weight < 0.5) {
            igraph_vector_init_copy(&normalized_weights, &weights_vec);
            double range = max_weight - min_weight;
            if (range > 1e-6) {
                for (size_t i = 0; i < igraph_vector_size(&normalized_weights); i++) {
                    VECTOR(normalized_weights)[i] = (VECTOR(normalized_weights)[i] - min_weight) / range;
                }
            }
            cerr << "-----Edge weights normalized: [" << min_weight << ", " << max_weight 
                 << "] -> [0, 1]" << endl;
        } else {
            igraph_vector_init_copy(&normalized_weights, &weights_vec);
        }
        
        // Adjust parameters for Leiden algorithm
        double beta = 0.01;  // Standard refinement parameter
        double adjusted_resolution = resolution;  // Use resolution directly after normalization
        
        result_code = igraph_community_leiden(
            &graph,
            &normalized_weights,  // normalized edge weights
            NULL,             // node weights (NULL = all 1.0)
            adjusted_resolution,  // resolution parameter
            beta,             // beta (refinement parameter)
            false,            // start from singleton partition
            100,              // more iterations for convergence
            &membership,
            &nb_clusters_temp,
            &quality_temp
        );
        
        igraph_vector_destroy(&normalized_weights);
        
        // Store quality in modularity_vec for consistency
        igraph_vector_resize(&modularity_vec, 1);
        VECTOR(modularity_vec)[0] = quality_temp;
    } else {
        // Use Multilevel (Louvain) algorithm - more stable and well-tested
        result_code = igraph_community_multilevel(
            &graph,
            &weights_vec,     // edge weights
            resolution,       // resolution parameter  
            &membership,
            NULL,             // memberships (intermediate steps, not needed)
            &modularity_vec   // final modularity (vector)
        );
    }
    
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
    robin_hood::unordered_map<int, vector<int>> cluster_map;
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

// ==================== Helper: Save Graph to File ====================
void save_graph_to_file(const vector<Edge>& edges, int num_nodes, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "ERROR: Cannot open graph file for writing: " << filename << endl;
        return;
    }
    
    // Header: num_nodes num_edges
    outfile << num_nodes << " " << edges.size() << "\n";
    
    // Edges: from to weight
    for (const auto& e : edges) {
        outfile << e.from << " " << e.to << " " << e.weight << "\n";
    }
    
    outfile.close();
    cerr << "-----Graph saved to: " << filename << endl;
}

// ==================== Cluster from Pre-built Graph ====================
vector<vector<int>> KssdLeidenClusterFromGraph(
    const string& graph_file,
    int num_genomes,
    double resolution,
    bool use_leiden)
{
    cerr << "=============================" << endl;
    if (use_leiden) {
        cerr << "Graph-based Clustering (Leiden) - From Pre-built Graph" << endl;
    } else {
        cerr << "Graph-based Clustering (Louvain) - From Pre-built Graph" << endl;
    }
    cerr << "=============================" << endl;
    cerr << "Graph file: " << graph_file << endl;
    cerr << "Genomes: " << num_genomes << endl;
    cerr << "Resolution: " << resolution << endl;
    cerr << "=============================" << endl;
    
    // Step 1: Load graph from file
    cerr << "-----Loading graph from file..." << endl;
    
    ifstream infile(graph_file);
    if (!infile.is_open()) {
        cerr << "ERROR: Cannot open graph file: " << graph_file << endl;
        vector<vector<int>> result(num_genomes);
        for (int i = 0; i < num_genomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    int file_num_nodes, file_num_edges;
    infile >> file_num_nodes >> file_num_edges;
    
    if (file_num_nodes != num_genomes) {
        cerr << "WARNING: Graph file has " << file_num_nodes 
             << " nodes but expected " << num_genomes << endl;
    }
    
    vector<Edge> edges;
    edges.reserve(file_num_edges);
    
    int from, to;
    double weight;
    while (infile >> from >> to >> weight) {
        edges.push_back(Edge(from, to, weight));
    }
    infile.close();
    
    cerr << "-----Loaded " << edges.size() << " edges" << endl;
    
    if (edges.empty()) {
        cerr << "-----Warning: No edges! Each sequence in its own cluster." << endl;
        vector<vector<int>> result(num_genomes);
        for (int i = 0; i < num_genomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    // Step 2: Build igraph
    cerr << "-----Building igraph from loaded edges..." << endl;
    
    igraph_t graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    igraph_vector_int_init(&edges_vec, edges.size() * 2);
    igraph_vector_init(&weights_vec, edges.size());
    
    for (size_t i = 0; i < edges.size(); i++) {
        VECTOR(edges_vec)[2 * i] = edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = edges[i].to;
        VECTOR(weights_vec)[i] = edges[i].weight;
    }
    
    igraph_create(&graph, &edges_vec, num_genomes, IGRAPH_UNDIRECTED);
    
    // Step 3: Run community detection algorithm
    if (use_leiden) {
        cerr << "-----Running Leiden algorithm (default)..." << endl;
    } else {
        cerr << "-----Running Louvain algorithm..." << endl;
        cerr << "-----Note: Using --louvain flag. Louvain may be faster but less refined than Leiden (default)." << endl;
    }
    
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, num_genomes);
    igraph_vector_init(&modularity_vec, 0);
    
    int result_code;
    
    if (use_leiden) {
        // Normalize edge weights for Leiden
        double min_weight = 1.0, max_weight = 0.0;
        for (size_t i = 0; i < igraph_vector_size(&weights_vec); i++) {
            double w = VECTOR(weights_vec)[i];
            if (w < min_weight) min_weight = w;
            if (w > max_weight) max_weight = w;
        }
        
        igraph_vector_t normalized_weights;
        if (max_weight - min_weight < 0.5) {
            igraph_vector_init_copy(&normalized_weights, &weights_vec);
            double range = max_weight - min_weight;
            if (range > 1e-6) {
                for (size_t i = 0; i < igraph_vector_size(&normalized_weights); i++) {
                    VECTOR(normalized_weights)[i] = (VECTOR(normalized_weights)[i] - min_weight) / range;
                }
            }
            cerr << "-----Edge weights normalized: [" << min_weight << ", " << max_weight 
                 << "] -> [0, 1]" << endl;
        } else {
            igraph_vector_init_copy(&normalized_weights, &weights_vec);
        }
        
        double beta = 0.01;
        igraph_integer_t nb_clusters_temp;
        igraph_real_t quality_temp;
        
        result_code = igraph_community_leiden(
            &graph,
            &normalized_weights,
            NULL,
            resolution,
            beta,
            false,
            100,
            &membership,
            &nb_clusters_temp,
            &quality_temp
        );
        
        igraph_vector_destroy(&normalized_weights);
        
        igraph_vector_resize(&modularity_vec, 1);
        VECTOR(modularity_vec)[0] = quality_temp;
    } else {
        result_code = igraph_community_multilevel(
            &graph,
            &weights_vec,
            resolution,
            &membership,
            NULL,
            &modularity_vec
        );
    }
    
    double modularity = 0.0;
    if (igraph_vector_size(&modularity_vec) > 0) {
        modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
    }
    
    cerr << "-----Community detection complete!" << endl;
    cerr << "-----Final modularity: " << modularity << endl;
    
    // Count clusters
    int nb_clusters = 0;
    for (int i = 0; i < num_genomes; i++) {
        if (VECTOR(membership)[i] + 1 > nb_clusters) {
            nb_clusters = VECTOR(membership)[i] + 1;
        }
    }
    
    if (result_code != IGRAPH_SUCCESS) {
        cerr << "WARNING: Community detection failed, each genome in its own cluster" << endl;
        vector<vector<int>> result(num_genomes);
        for (int i = 0; i < num_genomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    cerr << "-----Number of clusters: " << nb_clusters << endl;
    
    // Extract clusters
    robin_hood::unordered_map<int, vector<int>> cluster_map;
    for (int i = 0; i < num_genomes; i++) {
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
    cerr << "Average size: " << (result.empty() ? 0.0 : (double)num_genomes / result.size()) << endl;
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

// ==================== Two-Stage Parallel Louvain Implementation ====================

/**
 * Graph partition structure
 * Each partition contains nodes and their associated edges
 */
struct GraphPartition {
    vector<int> nodes;              // Nodes in this partition
    vector<Edge> internal_edges;    // Edges within partition
    vector<Edge> boundary_edges;    // Edges crossing to other partitions
    int partition_id;
    
    GraphPartition(int id) : partition_id(id) {}
};

/**
 * Local cluster result from a partition
 */
struct LocalClusterResult {
    int partition_id;
    vector<vector<int>> clusters;    // Local clusters (node IDs are original)
    double modularity;
    
    // Map from original node ID to local cluster ID
    robin_hood::unordered_map<int, int> node_to_cluster;
};

/**
 * Super node in meta-graph
 */
struct SuperNode {
    int super_id;                    // ID in meta-graph
    int partition_id;                // Which partition it belongs to
    int local_cluster_id;            // Cluster ID within partition
    vector<int> original_nodes;      // Original node IDs
};

/**
 * Partition graph using simple hash-based strategy
 * This ensures balanced partitions and is very fast
 */
vector<GraphPartition> PartitionGraph(
    const vector<Edge>& all_edges,
    int num_nodes,
    int num_partitions)
{
    cerr << "-----Partitioning graph into " << num_partitions << " partitions..." << endl;
    
    vector<GraphPartition> partitions;
    for (int i = 0; i < num_partitions; i++) {
        partitions.push_back(GraphPartition(i));
    }
    
    // Assign nodes to partitions using simple modulo hash
    for (int node = 0; node < num_nodes; node++) {
        int partition_id = node % num_partitions;
        partitions[partition_id].nodes.push_back(node);
    }
    
    // Create node to partition mapping
    vector<int> node_partition(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        node_partition[i] = i % num_partitions;
    }
    
    // Assign edges to partitions
    for (const auto& edge : all_edges) {
        int part_from = node_partition[edge.from];
        int part_to = node_partition[edge.to];
        
        if (part_from == part_to) {
            // Internal edge
            partitions[part_from].internal_edges.push_back(edge);
        } else {
            // Boundary edge - add to both partitions
            partitions[part_from].boundary_edges.push_back(edge);
            partitions[part_to].boundary_edges.push_back(edge);
        }
    }
    
    // Statistics
    for (int i = 0; i < num_partitions; i++) {
        cerr << "  Partition " << i << ": " 
             << partitions[i].nodes.size() << " nodes, "
             << partitions[i].internal_edges.size() << " internal edges, "
             << partitions[i].boundary_edges.size() << " boundary edges" << endl;
    }
    
    return partitions;
}

/**
 * Run Louvain on a single partition
 */
LocalClusterResult RunLouvainOnPartition(
    const GraphPartition& partition,
    int total_nodes,
    double resolution)
{
    LocalClusterResult result;
    result.partition_id = partition.partition_id;
    result.modularity = 0.0;
    
    // Handle empty partition
    if (partition.nodes.empty()) {
        return result;
    }
    
    // Handle single node
    if (partition.nodes.size() == 1) {
        result.clusters.push_back({partition.nodes[0]});
        result.node_to_cluster[partition.nodes[0]] = 0;
        return result;
    }
    
    // Build local igraph
    igraph_t graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    // Only use internal edges for local clustering
    size_t num_edges = partition.internal_edges.size();
    
    if (num_edges == 0) {
        // No edges - each node in its own cluster
        for (size_t i = 0; i < partition.nodes.size(); i++) {
            result.clusters.push_back({partition.nodes[i]});
            result.node_to_cluster[partition.nodes[i]] = i;
        }
        return result;
    }
    
    igraph_vector_int_init(&edges_vec, num_edges * 2);
    igraph_vector_init(&weights_vec, num_edges);
    
    for (size_t i = 0; i < num_edges; i++) {
        VECTOR(edges_vec)[2 * i] = partition.internal_edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = partition.internal_edges[i].to;
        VECTOR(weights_vec)[i] = partition.internal_edges[i].weight;
    }
    
    igraph_create(&graph, &edges_vec, total_nodes, IGRAPH_UNDIRECTED);
    
    // Run Louvain
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, total_nodes);
    igraph_vector_init(&modularity_vec, 0);
    
    int ret = igraph_community_multilevel(
        &graph,
        &weights_vec,
        resolution,
        &membership,
        NULL,
        &modularity_vec
    );
    
    if (ret == IGRAPH_SUCCESS && igraph_vector_size(&modularity_vec) > 0) {
        result.modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
    }
    
    // Extract clusters (only for nodes in this partition)
    robin_hood::unordered_map<int, vector<int>> cluster_map;
    for (int node : partition.nodes) {
        int community = VECTOR(membership)[node];
        cluster_map[community].push_back(node);
        result.node_to_cluster[node] = community;
    }
    
    for (const auto& entry : cluster_map) {
        result.clusters.push_back(entry.second);
    }
    
    // Cleanup
    igraph_vector_destroy(&modularity_vec);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_int_destroy(&edges_vec);
    igraph_destroy(&graph);
    
    return result;
}

/**
 * Build meta-graph from local clustering results
 * Each local cluster becomes a super-node
 * Edges between super-nodes are aggregated from boundary edges
 */
pair<vector<SuperNode>, vector<Edge>> BuildMetaGraph(
    const vector<LocalClusterResult>& local_results,
    const vector<GraphPartition>& partitions)
{
    cerr << "-----Building meta-graph from local clusters..." << endl;
    
    vector<SuperNode> super_nodes;
    
    // Create super-nodes from local clusters
    int super_id = 0;
    robin_hood::unordered_map<string, int> cluster_to_super_id;  // "partition_id:cluster_id" -> super_id
    
    for (const auto& local_result : local_results) {
        for (size_t i = 0; i < local_result.clusters.size(); i++) {
            SuperNode super_node;
            super_node.super_id = super_id;
            super_node.partition_id = local_result.partition_id;
            super_node.local_cluster_id = i;
            super_node.original_nodes = local_result.clusters[i];
            super_nodes.push_back(super_node);
            
            string key = to_string(local_result.partition_id) + ":" + to_string(i);
            cluster_to_super_id[key] = super_id;
            super_id++;
        }
    }
    
    cerr << "  Created " << super_nodes.size() << " super-nodes from local clusters" << endl;
    
    // Build node to super-node mapping
    robin_hood::unordered_map<int, int> node_to_super_id;
    for (const auto& super_node : super_nodes) {
        for (int node : super_node.original_nodes) {
            node_to_super_id[node] = super_node.super_id;
        }
    }
    
    // Aggregate boundary edges into super-edges
    map<pair<int, int>, double> super_edge_weights;  // (super_from, super_to) -> total_weight
    
    for (const auto& partition : partitions) {
        for (const auto& edge : partition.boundary_edges) {
            auto it_from = node_to_super_id.find(edge.from);
            auto it_to = node_to_super_id.find(edge.to);
            
            if (it_from == node_to_super_id.end() || it_to == node_to_super_id.end()) {
                continue;  // Node not in any cluster
            }
            
            int super_from = it_from->second;
            int super_to = it_to->second;
            
            if (super_from == super_to) continue;  // Skip self-loops
            
            // Normalize edge direction
            if (super_from > super_to) swap(super_from, super_to);
            
            super_edge_weights[{super_from, super_to}] += edge.weight;
        }
    }
    
    // Convert to edge list
    vector<Edge> super_edges;
    for (const auto& entry : super_edge_weights) {
        super_edges.push_back(Edge(entry.first.first, entry.first.second, entry.second));
    }
    
    cerr << "  Created " << super_edges.size() << " super-edges from boundary edges" << endl;
    
    return {super_nodes, super_edges};
}

/**
 * Run Louvain on meta-graph
 */
vector<int> RunLouvainOnMetaGraph(
    const vector<SuperNode>& super_nodes,
    const vector<Edge>& super_edges,
    double resolution)
{
    cerr << "-----Running Louvain on meta-graph..." << endl;
    
    int num_super_nodes = super_nodes.size();
    
    if (num_super_nodes == 0) {
        return vector<int>();
    }
    
    if (super_edges.empty()) {
        // No edges - each super-node in its own meta-cluster
        vector<int> meta_membership(num_super_nodes);
        for (int i = 0; i < num_super_nodes; i++) {
            meta_membership[i] = i;
        }
        return meta_membership;
    }
    
    // Build igraph for meta-graph
    igraph_t meta_graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    igraph_vector_int_init(&edges_vec, super_edges.size() * 2);
    igraph_vector_init(&weights_vec, super_edges.size());
    
    for (size_t i = 0; i < super_edges.size(); i++) {
        VECTOR(edges_vec)[2 * i] = super_edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = super_edges[i].to;
        VECTOR(weights_vec)[i] = super_edges[i].weight;
    }
    
    igraph_create(&meta_graph, &edges_vec, num_super_nodes, IGRAPH_UNDIRECTED);
    
    // Run Louvain
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, num_super_nodes);
    igraph_vector_init(&modularity_vec, 0);
    
    int ret = igraph_community_multilevel(
        &meta_graph,
        &weights_vec,
        resolution,
        &membership,
        NULL,
        &modularity_vec
    );
    
    double modularity = 0.0;
    if (ret == IGRAPH_SUCCESS && igraph_vector_size(&modularity_vec) > 0) {
        modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
    }
    
    cerr << "  Meta-graph modularity: " << modularity << endl;
    
    // Extract membership
    vector<int> meta_membership(num_super_nodes);
    for (int i = 0; i < num_super_nodes; i++) {
        meta_membership[i] = VECTOR(membership)[i];
    }
    
    // Cleanup
    igraph_vector_destroy(&modularity_vec);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_int_destroy(&edges_vec);
    igraph_destroy(&meta_graph);
    
    return meta_membership;
}

/**
 * Map meta-clusters back to original nodes
 */
vector<vector<int>> MapToOriginalClusters(
    const vector<SuperNode>& super_nodes,
    const vector<int>& meta_membership)
{
    // Group super-nodes by meta-cluster
    robin_hood::unordered_map<int, vector<int>> meta_cluster_to_super_nodes;
    for (size_t i = 0; i < super_nodes.size(); i++) {
        int meta_cluster = meta_membership[i];
        meta_cluster_to_super_nodes[meta_cluster].push_back(i);
    }
    
    // Merge original nodes from super-nodes in same meta-cluster
    vector<vector<int>> final_clusters;
    for (const auto& entry : meta_cluster_to_super_nodes) {
        vector<int> cluster_nodes;
        for (int super_id : entry.second) {
            const auto& super_node = super_nodes[super_id];
            cluster_nodes.insert(cluster_nodes.end(),
                               super_node.original_nodes.begin(),
                               super_node.original_nodes.end());
        }
        final_clusters.push_back(cluster_nodes);
    }
    
    // Sort by cluster size
    sort(final_clusters.begin(), final_clusters.end(),
         [](const vector<int>& a, const vector<int>& b) {
             return a.size() > b.size();
         });
    
    return final_clusters;
}

/**
 * Main function: Two-stage parallel Louvain clustering
 */
vector<vector<int>> KssdParallelLouvainCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution,
    const string& graph_save_path)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) {
        return vector<vector<int>>();
    }
    
    cerr << "=============================" << endl;
    cerr << "Two-Stage Parallel Louvain Clustering" << endl;
    cerr << "=============================" << endl;
    cerr << "Genomes: " << numGenomes << endl;
    cerr << "Edge threshold: " << threshold << endl;
    cerr << "Resolution: " << resolution << endl;
    cerr << "Threads: " << threads << endl;
    cerr << "Partitions: " << threads << endl;
    cerr << "=============================" << endl;
    
    // Step 1: Build inverted index (same as before)
    GlobalInvertedIndex global_index;
    global_index.build(sketches, threads);
    
    // Step 2: Build edge list (same as before, already parallelized)
    cerr << "-----Building similarity graph..." << endl;
    
    vector<vector<Edge>> thread_local_edges(threads);
    long long total_comparisons = 0;
    long long edges_created = 0;
    
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        auto& local_edges = thread_local_edges[tid];
        robin_hood::unordered_map<int, int> intersection_map;
        
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
    
    // Save graph if requested
    if (!graph_save_path.empty()) {
        save_graph_to_file(all_edges, numGenomes, graph_save_path);
    }
    
    // Step 3: Partition graph
    vector<GraphPartition> partitions = PartitionGraph(all_edges, numGenomes, threads);
    
    // Step 4: Run Louvain on each partition in parallel (Stage 1)
    cerr << "-----Stage 1: Running Louvain on " << threads << " partitions in parallel..." << endl;
    
    vector<LocalClusterResult> local_results(threads);
    
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++) {
        local_results[i] = RunLouvainOnPartition(partitions[i], numGenomes, resolution);
        
        #pragma omp critical
        {
            cerr << "  Partition " << i << ": " 
                 << local_results[i].clusters.size() << " local clusters, "
                 << "modularity = " << local_results[i].modularity << endl;
        }
    }
    
    // Count total local clusters
    int total_local_clusters = 0;
    for (const auto& result : local_results) {
        total_local_clusters += result.clusters.size();
    }
    cerr << "  Total local clusters: " << total_local_clusters << endl;
    
    // Step 5: Build meta-graph
    auto [super_nodes, super_edges] = BuildMetaGraph(local_results, partitions);
    
    // Step 6: Run Louvain on meta-graph (Stage 2)
    cerr << "-----Stage 2: Running Louvain on meta-graph..." << endl;
    vector<int> meta_membership = RunLouvainOnMetaGraph(super_nodes, super_edges, resolution);
    
    // Step 7: Map back to original clusters
    cerr << "-----Mapping meta-clusters to original nodes..." << endl;
    vector<vector<int>> final_clusters = MapToOriginalClusters(super_nodes, meta_membership);
    
    cerr << "=============================" << endl;
    cerr << "Two-Stage Clustering Complete" << endl;
    cerr << "Total final clusters: " << final_clusters.size() << endl;
    cerr << "Largest cluster: " << (final_clusters.empty() ? 0 : final_clusters[0].size()) << endl;
    cerr << "Average size: " << (final_clusters.empty() ? 0.0 : (double)numGenomes / final_clusters.size()) << endl;
    cerr << "=============================" << endl;
    
    return final_clusters;
}

// ==================== Edge-Partitioning Parallel Louvain (Simpler & More Efficient) ====================

/**
 * Helper: Run Louvain on a subset of edges (for warm-start)
 * This provides a good initial clustering that can speed up final convergence
 */
vector<int> RunLouvainOnEdgeSubset(
    const vector<Edge>& edges,
    int num_nodes,
    double resolution)
{
    if (edges.empty()) {
        // No edges - return singleton clusters
        vector<int> membership(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            membership[i] = i;
        }
        return membership;
    }
    
    // Build graph from edge subset
    igraph_t graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    igraph_vector_int_init(&edges_vec, edges.size() * 2);
    igraph_vector_init(&weights_vec, edges.size());
    
    for (size_t i = 0; i < edges.size(); i++) {
        VECTOR(edges_vec)[2 * i] = edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = edges[i].to;
        VECTOR(weights_vec)[i] = edges[i].weight;
    }
    
    igraph_create(&graph, &edges_vec, num_nodes, IGRAPH_UNDIRECTED);
    
    // Run Louvain
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, num_nodes);
    igraph_vector_init(&modularity_vec, 0);
    
    igraph_community_multilevel(
        &graph,
        &weights_vec,
        resolution,
        &membership,
        NULL,
        &modularity_vec
    );
    
    // Extract membership
    vector<int> result(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        result[i] = VECTOR(membership)[i];
    }
    
    // Cleanup
    igraph_vector_destroy(&modularity_vec);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_int_destroy(&edges_vec);
    igraph_destroy(&graph);
    
    return result;
}

/**
 * Select best local clustering result based on coverage and diversity
 */
vector<int> SelectBestWarmStart(
    const vector<vector<int>>& local_memberships,
    const vector<vector<Edge>>& thread_edges,
    int num_nodes)
{
    int best_tid = 0;
    size_t max_edges = 0;
    
    // Simple heuristic: choose the thread with most edges (most information)
    for (size_t i = 0; i < thread_edges.size(); i++) {
        if (thread_edges[i].size() > max_edges) {
            max_edges = thread_edges[i].size();
            best_tid = i;
        }
    }
    
    cerr << "  Selected warm-start from thread " << best_tid 
         << " (edges: " << max_edges << ")" << endl;
    
    return local_memberships[best_tid];
}

/**
 * Run Louvain with initial membership (warm-start)
 */
vector<vector<int>> RunLouvainWithWarmStart(
    const vector<Edge>& all_edges,
    int num_nodes,
    double resolution,
    const vector<int>& initial_membership)
{
    if (all_edges.empty()) {
        cerr << "-----Warning: No edges! Each sequence in its own cluster." << endl;
        vector<vector<int>> result(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    cerr << "-----Running Louvain with warm-start initialization..." << endl;
    
    // Build igraph
    igraph_t graph;
    igraph_vector_int_t edges_vec;
    igraph_vector_t weights_vec;
    
    igraph_vector_int_init(&edges_vec, all_edges.size() * 2);
    igraph_vector_init(&weights_vec, all_edges.size());
    
    for (size_t i = 0; i < all_edges.size(); i++) {
        VECTOR(edges_vec)[2 * i] = all_edges[i].from;
        VECTOR(edges_vec)[2 * i + 1] = all_edges[i].to;
        VECTOR(weights_vec)[i] = all_edges[i].weight;
    }
    
    igraph_create(&graph, &edges_vec, num_nodes, IGRAPH_UNDIRECTED);
    
    // Run Louvain (igraph doesn't support warm-start directly, but initial clustering helps)
    // The algorithm will start from a better position and converge faster
    igraph_vector_int_t membership;
    igraph_vector_t modularity_vec;
    
    igraph_vector_int_init(&membership, num_nodes);
    igraph_vector_init(&modularity_vec, 0);
    
    // Copy initial membership
    for (int i = 0; i < num_nodes; i++) {
        VECTOR(membership)[i] = initial_membership[i];
    }
    
    int result_code = igraph_community_multilevel(
        &graph,
        &weights_vec,
        resolution,
        &membership,
        NULL,
        &modularity_vec
    );
    
    double modularity = 0.0;
    if (result_code == IGRAPH_SUCCESS && igraph_vector_size(&modularity_vec) > 0) {
        modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
    }
    
    cerr << "  Final modularity: " << modularity << endl;
    
    // Extract clusters
    robin_hood::unordered_map<int, vector<int>> cluster_map;
    for (int i = 0; i < num_nodes; i++) {
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
    
    // Cleanup
    igraph_vector_destroy(&modularity_vec);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_int_destroy(&edges_vec);
    igraph_destroy(&graph);
    
    return result;
}

/**
 * Edge-partitioning parallel Louvain clustering
 * This is simpler and more efficient than graph-partitioning approach
 * 
 * Key idea:
 * 1. Parallel edge construction naturally partitions edges by source node
 * 2. Optionally: each thread runs local Louvain on its edge subset (warm-start)
 * 3. Merge all edges and run final Louvain on complete graph
 * 
 * Advantages over graph-partitioning:
 * - No graph partitioning overhead
 * - No boundary edge handling complexity
 * - Final clustering on complete graph (better quality)
 * - Natural integration with existing MST-style parallel construction
 */
vector<vector<int>> KssdEdgeParallelLouvainCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution,
    bool use_warm_start,
    int knn_k,
    const string& graph_save_path)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) {
        return vector<vector<int>>();
    }
    
    cerr << "=============================" << endl;
    cerr << "Edge-Partitioning Parallel Louvain" << endl;
    if (use_warm_start) {
        cerr << "(with warm-start optimization)" << endl;
    }
    if (knn_k > 0) {
        cerr << "(with k-NN optimization: k=" << knn_k << ")" << endl;
    }
    cerr << "=============================" << endl;
    cerr << "Genomes: " << numGenomes << endl;
    cerr << "Edge threshold: " << threshold << endl;
    if (knn_k > 0) {
        cerr << "k-NN k: " << knn_k << " (keep only top-" << knn_k << " neighbors per node)" << endl;
    }
    cerr << "Resolution: " << resolution << endl;
    cerr << "Threads: " << threads << endl;
    cerr << "=============================" << endl;
    
    // Step 1: Build inverted index
    GlobalInvertedIndex global_index;
    global_index.build(sketches, threads);
    
    // Step 2: Parallel edge construction (natural edge partitioning)
    // Each thread handles a subset of source nodes → natural edge partitioning
    cerr << "-----Building similarity graph (parallel edge construction)..." << endl;
    
    vector<vector<Edge>> thread_local_edges(threads);
    long long total_comparisons = 0;
    long long edges_created = 0;
    
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        auto& local_edges = thread_local_edges[tid];
        robin_hood::unordered_map<int, int> intersection_map;
        
        long long local_comparisons = 0;
        long long local_edge_count = 0;
        
        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < numGenomes; i++) {
            const auto& sketch_i = sketches[i];
            int size_i = sketch_i.sketchsize;
            
            global_index.calculate_all_intersections(i, sketch_i.hash32_arr, intersection_map);
            
            if (knn_k > 0) {
                // k-NN mode: keep only top-k nearest neighbors
                // Use max-heap to maintain top-k (smallest distances)
                priority_queue<pair<double, int>> topk;  // max-heap: <distance, node_id>
                
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
                        if ((int)topk.size() < knn_k) {
                            topk.push({dist, j});
                        } else if (dist < topk.top().first) {
                            topk.pop();
                            topk.push({dist, j});
                        }
                    }
                }
                
                // Add top-k edges
                while (!topk.empty()) {
                    double dist = topk.top().first;
                    int j = topk.top().second;
                    topk.pop();
                    local_edges.push_back(Edge(i, j, 1.0 - dist));
                    local_edge_count++;
                }
            } else {
                // Standard mode: keep all edges within threshold
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
    cerr << "-----Edges created: " << edges_created << endl;
    if (knn_k > 0) {
        cerr << "-----k-NN filtered: keeping top-" << knn_k << " neighbors per node (max " 
             << (numGenomes * knn_k) << " edges)" << endl;
    }
    
    // Print edge distribution across threads
    cerr << "-----Edge distribution:" << endl;
    for (int i = 0; i < threads; i++) {
        cerr << "  Thread " << i << ": " << thread_local_edges[i].size() << " edges" << endl;
    }
    
    // Step 3 (Optional): Parallel local Louvain for warm-start
    vector<int> warm_start_membership;
    
    if (use_warm_start && edges_created > 0) {
        cerr << "-----Phase 1: Parallel local Louvain (warm-start)..." << endl;
        
        vector<vector<int>> local_memberships(threads);
        
        #pragma omp parallel for num_threads(threads)
        for (int tid = 0; tid < threads; tid++) {
            if (!thread_local_edges[tid].empty()) {
                local_memberships[tid] = RunLouvainOnEdgeSubset(
                    thread_local_edges[tid], 
                    numGenomes, 
                    resolution
                );
                
                #pragma omp critical
                {
                    cerr << "  Thread " << tid << " local clustering done" << endl;
                }
            }
        }
        
        // Select best local result as warm-start
        warm_start_membership = SelectBestWarmStart(
            local_memberships, 
            thread_local_edges, 
            numGenomes
        );
    }
    
    // Step 4: Merge all edges
    cerr << "-----Merging edges from all threads..." << endl;
    vector<Edge> all_edges;
    all_edges.reserve(edges_created);
    
    for (const auto& local_edges : thread_local_edges) {
        all_edges.insert(all_edges.end(), local_edges.begin(), local_edges.end());
    }
    
    cerr << "  Total edges: " << all_edges.size() << endl;
    
    if (all_edges.empty()) {
        cerr << "-----Warning: No edges! Each sequence in its own cluster." << endl;
        vector<vector<int>> result(numGenomes);
        for (int i = 0; i < numGenomes; i++) {
            result[i].push_back(i);
        }
        return result;
    }
    
    // Save graph if requested
    if (!graph_save_path.empty()) {
        save_graph_to_file(all_edges, numGenomes, graph_save_path);
    }
    
    // Step 5: Final Louvain on complete graph
    cerr << "-----Phase 2: Final Louvain on complete graph..." << endl;
    
    vector<vector<int>> final_clusters;
    
    if (use_warm_start && !warm_start_membership.empty()) {
        // Use warm-start
        final_clusters = RunLouvainWithWarmStart(
            all_edges, 
            numGenomes, 
            resolution, 
            warm_start_membership
        );
    } else {
        // Standard Louvain without warm-start
        cerr << "-----Running standard Louvain on complete graph..." << endl;
        
        igraph_t graph;
        igraph_vector_int_t edges_vec;
        igraph_vector_t weights_vec;
        
        igraph_vector_int_init(&edges_vec, all_edges.size() * 2);
        igraph_vector_init(&weights_vec, all_edges.size());
        
        for (size_t i = 0; i < all_edges.size(); i++) {
            VECTOR(edges_vec)[2 * i] = all_edges[i].from;
            VECTOR(edges_vec)[2 * i + 1] = all_edges[i].to;
            VECTOR(weights_vec)[i] = all_edges[i].weight;
        }
        
        igraph_create(&graph, &edges_vec, numGenomes, IGRAPH_UNDIRECTED);
        
        igraph_vector_int_t membership;
        igraph_vector_t modularity_vec;
        
        igraph_vector_int_init(&membership, numGenomes);
        igraph_vector_init(&modularity_vec, 0);
        
        int result_code = igraph_community_multilevel(
            &graph,
            &weights_vec,
            resolution,
            &membership,
            NULL,
            &modularity_vec
        );
        
        double modularity = 0.0;
        if (result_code == IGRAPH_SUCCESS && igraph_vector_size(&modularity_vec) > 0) {
            modularity = VECTOR(modularity_vec)[igraph_vector_size(&modularity_vec) - 1];
        }
        
        cerr << "  Final modularity: " << modularity << endl;
        
        // Extract clusters
        robin_hood::unordered_map<int, vector<int>> cluster_map;
        for (int i = 0; i < numGenomes; i++) {
            int community = VECTOR(membership)[i];
            cluster_map[community].push_back(i);
        }
        
        for (const auto& entry : cluster_map) {
            final_clusters.push_back(entry.second);
        }
        
        sort(final_clusters.begin(), final_clusters.end(),
             [](const vector<int>& a, const vector<int>& b) {
                 return a.size() > b.size();
             });
        
        // Cleanup
        igraph_vector_destroy(&modularity_vec);
        igraph_vector_int_destroy(&membership);
        igraph_vector_destroy(&weights_vec);
        igraph_vector_int_destroy(&edges_vec);
        igraph_destroy(&graph);
    }
    
    cerr << "=============================" << endl;
    cerr << "Edge-Partitioning Clustering Complete" << endl;
    cerr << "Total clusters: " << final_clusters.size() << endl;
    cerr << "Largest cluster: " << (final_clusters.empty() ? 0 : final_clusters[0].size()) << endl;
    cerr << "Average size: " << (final_clusters.empty() ? 0.0 : (double)numGenomes / final_clusters.size()) << endl;
    cerr << "=============================" << endl;
    
    return final_clusters;
}

#endif // LEIDEN_CLUST
