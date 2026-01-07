#ifdef LEIDEN_CLUST
#ifndef H_LEIDEN
#define H_LEIDEN

#include "SketchInfo.h"
#include <vector>

using std::vector;

// Leiden-style clustering for KSSD sketches using inverted index + Label Propagation Algorithm
// No external dependencies - pure C++ implementation
// Parameters:
//   - sketches: input sketches
//   - sketch_func_id: sketch function ID (currently not used)
//   - threshold: similarity threshold for edge creation (0-1, lower = more similar)
//   - threads: number of threads for parallel computation
//   - kmer_size: k-mer size for Mash distance calculation
//   - resolution: controls cluster granularity (higher = more clusters)
//   - max_iterations: maximum iterations for label propagation (default 100)
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    double resolution = 1.0,
    bool use_modularity = true);

#endif
#endif

