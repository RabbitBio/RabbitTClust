#ifdef GREEDY_CLUST
#ifndef H_GREEDY
#define H_GREEDY

#include "SketchInfo.h"
#include <vector>

using std::vector;

// 原始的MinHash/KSSD贪心聚类
vector<vector<int>> greedyCluster(vector<SketchInfo>& sketches, int sketch_func_id, double threshold, int threads);

// KSSD贪心聚类（当前版本）
vector<vector<int>> KssdgreedyCluster(vector<KssdSketchInfo>& sketches, int sketch_func_id, double threshold, int threads);

// 基于倒排索引的KSSD贪心增量聚类（推荐）
// 外层串行 + 内层并行策略，使用动态倒排索引加速交集计算
vector<vector<int>> KssdGreedyClusterWithInvertedIndex(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size);

// 批处理版本（实验性）
// 批内并行处理，可能改变结果但提供更好的并行度
vector<vector<int>> KssdGreedyClusterWithInvertedIndexBatched(
    vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size,
    int batch_size);

// 基于倒排索引的MinHash贪心增量聚类
// 外层串行 + 内层并行策略，使用动态倒排索引加速交集计算
vector<vector<int>> MinHashGreedyClusterWithInvertedIndex(
    vector<SketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size);

#endif
#endif

