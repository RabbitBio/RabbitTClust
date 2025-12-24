#ifndef CLUSTER_POSTPROCESS_H
#define CLUSTER_POSTPROCESS_H

#include <vector>
#include "MST.h"
#include "SketchInfo.h"

// Build per-cluster candidate lists by collapsing "near-duplicate" nodes connected by
// forest edges with dist <= dedup_dist.
//
// - If dedup_dist <= 0, this is a no-op: candidates == clusters and node_to_rep[i] = i.
// - node_to_rep maps every node id -> chosen representative id of its dedup group.
std::vector<std::vector<int>> build_dedup_candidates_per_cluster(
    const std::vector<std::vector<int>>& clusters,
    const std::vector<EdgeInfo>& forest,
    const std::vector<KssdSketchInfo>& sketches,
    bool sketchByFile,
    double dedup_dist,
    std::vector<int>& node_to_rep);

// Select up to k representative ids per cluster, aiming to cover diversity.
// Uses farthest-first traversal (k-center heuristic) on the tree metric induced by the forest.
//
// clusters_original: original connected-components on forest
// candidates_per_cluster: selectable ids (e.g., output of build_dedup_candidates_per_cluster)
// node_to_rep: mapping from any node in clusters_original to a valid candidate id
std::vector<std::vector<int>> select_k_reps_per_cluster_tree(
    const std::vector<std::vector<int>>& clusters_original,
    const std::vector<std::vector<int>>& candidates_per_cluster,
    const std::vector<EdgeInfo>& forest,
    int N,
    const std::vector<int>& node_to_rep,
    int k);

#endif


