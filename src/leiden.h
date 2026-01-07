#ifdef LEIDEN_CLUST
#ifndef H_LEIDEN
#define H_LEIDEN

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
    bool use_leiden = false);

#endif
#endif

