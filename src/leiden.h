#ifdef LEIDEN_CLUST
#ifndef H_LEIDEN
#define H_LEIDEN

#include "SketchInfo.h"
#include <vector>

using std::vector;

// Leiden clustering for KSSD sketches
// Currently uses inverted index approach (placeholder, will be replaced with true Leiden algorithm)
vector<vector<int>> KssdLeidenCluster(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    int threads,
    int kmer_size);

// Future: true Leiden algorithm with modularity optimization
// vector<vector<int>> KssdLeidenClusterTrue(
//     vector<KssdSketchInfo>& sketches,
//     int sketch_func_id,
//     int threads,
//     int kmer_size,
//     double resolution = 1.0);

#endif
#endif

