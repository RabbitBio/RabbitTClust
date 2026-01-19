#ifdef GREEDY_CLUST
#ifndef H_GREEDY
#define H_GREEDY

#include "SketchInfo.h"
#include <vector>

using std::vector;

vector<vector<int>> greedyCluster(vector<SketchInfo>& sketches, int sketch_func_id, double threshold, int threads);

vector<vector<int>> KssdgreedyCluster(vector<KssdSketchInfo>& sketches, int sketch_func_id, double threshold, int threads);

vector<vector<int>> KssdGreedyClusterWithInvertedIndex(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size);

vector<vector<int>> KssdGreedyClusterWithInvertedIndexBatched(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    int batch_size);

vector<vector<int>> MinHashGreedyClusterWithInvertedIndex(
    vector<SketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size);

#include "phmap.h"

struct KssdClusterState {
    vector<int> representative_ids;
    vector<KssdSketchInfo> representatives;
    vector<vector<int>> clusters;
    vector<KssdSketchInfo> all_sketches;
    KssdParameters params;
    double threshold;
    int kmer_size;

    phmap::flat_hash_map<uint32_t, vector<int>> inverted_index;

    bool save(const string& filepath) const;
    bool load(const string& filepath);

    void build_inverted_index();
    void update_inverted_index(int rep_idx);
};

vector<vector<int>> KssdIncrementalCluster(
    KssdClusterState& state,
    vector<KssdSketchInfo>& new_sketches,
    int threads);

KssdClusterState KssdInitialClusterWithState(
    vector<KssdSketchInfo>& sketches,
    const KssdParameters& params,
    double threshold,
    int threads,
    int kmer_size);

#endif
#endif
