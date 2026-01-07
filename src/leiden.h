#ifdef LEIDEN_CLUST
#ifndef H_LEIDEN
#define H_LEIDEN

/**
 * @file leiden.h
 * @brief Graph-based clustering for genome sequences using Louvain/Leiden algorithms
 * @author Tong Zhang
 * @date 2026-01-06
 * 
 * This module provides graph-based clustering functionality using inverted index
 * for efficient similarity graph construction and igraph library for community detection.
 * Supports both Louvain (default, stable) and Leiden (experimental) algorithms.
 */

#include "SketchInfo.h"
#include <vector>

using std::vector;

// Graph-based clustering for KSSD sketches using inverted index
// Supports both Louvain (default, stable) and Leiden (optional, refined) algorithms
// Parameters:
//   - sketches: input sketches
//   - sketch_func_id: sketch function ID (currently not used)
//   - threshold: similarity threshold for edge creation (0-1, lower = more similar)
//   - threads: number of threads for parallel computation
//   - kmer_size: k-mer size for Mash distance calculation
//   - resolution: controls cluster granularity (higher = more clusters)
//   - use_leiden: if true, use Leiden algorithm; if false (default), use Louvain
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution = 1.0,
    bool use_leiden = false,
    const std::string& graph_save_path = "");

// Cluster from pre-built graph (fast resolution adjustment without rebuilding graph)
// This function loads a pre-saved graph and runs community detection with specified resolution
// Parameters:
//   - graph_file: path to saved graph file
//   - num_genomes: number of genomes (nodes in graph)
//   - resolution: controls cluster granularity (higher = more clusters)
//   - use_leiden: if true, use Leiden algorithm; if false (default), use Louvain
vector<vector<int>> KssdLeidenClusterFromGraph(
    const std::string& graph_file,
    int num_genomes,
    double resolution = 1.0,
    bool use_leiden = false);

#endif
#endif

