#include "cluster_postprocess.h"

#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include "UnionFind.h"

using std::pair;
using std::vector;

static inline uint64_t get_seq_len(const KssdSketchInfo& s, bool sketchByFile) {
  return sketchByFile ? s.totalSeqLength : (uint64_t)s.seqInfo.length;
}

static vector<vector<pair<int, double>>> build_adj_from_forest(int N, const vector<EdgeInfo>& forest) {
  vector<vector<pair<int, double>>> adj(N);
  adj.reserve(N);
  for (const auto& e : forest) {
    if (e.preNode < 0 || e.preNode >= N || e.sufNode < 0 || e.sufNode >= N) continue;
    adj[e.preNode].push_back({e.sufNode, e.dist});
    adj[e.sufNode].push_back({e.preNode, e.dist});
  }
  return adj;
}

static void distances_from(int start,
                           const vector<vector<pair<int, double>>>& ladj,
                           vector<double>& dist) {
  const int m = (int)ladj.size();
  dist.assign(m, -1.0);
  vector<int> st;
  vector<int> parent(m, -1);
  st.push_back(start);
  dist[start] = 0.0;
  parent[start] = start;
  while (!st.empty()) {
    int u = st.back();
    st.pop_back();
    for (const auto& vw : ladj[u]) {
      int v = vw.first;
      double w = vw.second;
      if (v == parent[u]) continue;
      parent[v] = u;
      dist[v] = dist[u] + w;
      st.push_back(v);
    }
  }
}

std::vector<std::vector<int>> build_dedup_candidates_per_cluster(
    const std::vector<std::vector<int>>& clusters,
    const std::vector<EdgeInfo>& forest,
    const std::vector<KssdSketchInfo>& sketches,
    bool sketchByFile,
    double dedup_dist,
    std::vector<int>& node_to_rep) {
  const int N = (int)sketches.size();
  node_to_rep.assign(N, -1);

  // No-op: identity mapping, candidates = clusters
  if (dedup_dist <= 0) {
    for (int i = 0; i < N; i++) node_to_rep[i] = i;
    return clusters;
  }

  // Build dedup subgraph adjacency list (only edges with dist <= dedup_dist)
  vector<vector<pair<int, double>>> dedupAdj(N);
  UnionFind uf(N);
  for (const auto& e : forest) {
    if (e.dist <= dedup_dist) {
      uf.merge(e.preNode, e.sufNode);
      if (e.preNode >= 0 && e.preNode < N && e.sufNode >= 0 && e.sufNode < N) {
        dedupAdj[e.preNode].push_back({e.sufNode, e.dist});
        dedupAdj[e.sufNode].push_back({e.preNode, e.dist});
      }
    }
  }

  // Group members by UF root
  std::unordered_map<int, vector<int>> groups;
  for (int i = 0; i < N; i++) {
    int r = uf.find(i);
    groups[r].push_back(i);
  }

  vector<int> bestRep(N, -1);

  // For each UF group, select tree-medoid (node with minimum total tree distance to others in group)
  for (const auto& kv : groups) {
    const vector<int>& members = kv.second;
    if (members.size() == 1) {
      bestRep[kv.first] = members[0];
      continue;
    }

    int chosenRep = members[0];
    double minTotalDist = std::numeric_limits<double>::infinity();
    uint64_t chosenLen = 0;

    for (int cand : members) {
      // Compute tree distances from cand to all nodes in N
      vector<double> dist;
      distances_from(cand, dedupAdj, dist);

      // Sum distances to other members in this group
      double totalDist = 0.0;
      for (int m : members) {
        if (m != cand && dist[m] >= 0) {
          totalDist += dist[m];
        }
      }

      uint64_t candLen = get_seq_len(sketches[cand], sketchByFile);
      // Choose if: smaller total distance, or tie + longer sequence, or tie + smaller id
      if (totalDist < minTotalDist ||
          (totalDist == minTotalDist && (candLen > chosenLen || (candLen == chosenLen && cand < chosenRep)))) {
        minTotalDist = totalDist;
        chosenRep = cand;
        chosenLen = candLen;
      }
    }

    bestRep[kv.first] = chosenRep;
  }

  for (int i = 0; i < N; i++) {
    int r = uf.find(i);
    node_to_rep[i] = (bestRep[r] == -1) ? i : bestRep[r];
  }
  uf.clear();

  vector<vector<int>> candidates;
  candidates.reserve(clusters.size());
  for (const auto& cl : clusters) {
    std::unordered_set<int> seen;
    vector<int> cand;
    cand.reserve(cl.size());
    for (int node : cl) {
      if (node < 0 || node >= N) continue;
      int rep = node_to_rep[node];
      if (rep < 0) rep = node;
      if (seen.insert(rep).second) cand.push_back(rep);
    }
    std::sort(cand.begin(), cand.end());
    candidates.push_back(std::move(cand));
  }
  return candidates;
}

static int farthest_node(int start,
                         const vector<vector<pair<int, double>>>& ladj,
                         vector<double>& dist_out) {
  distances_from(start, ladj, dist_out);
  int far = start;
  double best = -1.0;
  for (int i = 0; i < (int)dist_out.size(); i++) {
    if (dist_out[i] > best) {
      best = dist_out[i];
      far = i;
    }
  }
  return far;
}

std::vector<std::vector<int>> select_k_reps_per_cluster_tree(
    const std::vector<std::vector<int>>& clusters_original,
    const std::vector<std::vector<int>>& candidates_per_cluster,
    const std::vector<EdgeInfo>& forest,
    int N,
    const std::vector<int>& node_to_rep,
    int k) {
  vector<vector<int>> reps;
  reps.reserve(clusters_original.size());
  if (k <= 0) {
    reps.resize(clusters_original.size());
    return reps;
  }

  const auto adj = build_adj_from_forest(N, forest);
  const double INF = std::numeric_limits<double>::infinity();

  for (size_t ci = 0; ci < clusters_original.size(); ci++) {
    const auto& compNodes = clusters_original[ci];
    const auto& candidates = candidates_per_cluster[ci];

    if (candidates.empty()) {
      reps.push_back({});
      continue;
    }
    if ((int)candidates.size() <= k) {
      reps.push_back(candidates);
      continue;
    }

    // Localize node ids to 0..m-1 for fast traversal
    const int m = (int)compNodes.size();
    std::unordered_map<int, int> idx;
    idx.reserve(m * 2);
    vector<int> nodes = compNodes;
    for (int i = 0; i < m; i++) idx[nodes[i]] = i;

    vector<vector<pair<int, double>>> ladj(m);
    for (int i = 0; i < m; i++) {
      int u = nodes[i];
      for (const auto& vw : adj[u]) {
        int v = vw.first;
        double w = vw.second;
        auto it = idx.find(v);
        if (it != idx.end()) {
          ladj[i].push_back({it->second, w});
        }
      }
    }

    // Compute diameter endpoints on the full component tree
    vector<double> distTmp;
    int u = farthest_node(0, ladj, distTmp);
    int v = farthest_node(u, ladj, distTmp);

    // Candidate membership set
    std::unordered_set<int> candSet;
    candSet.reserve(candidates.size() * 2);
    for (int c : candidates) candSet.insert(c);

    auto map_to_candidate = [&](int nodeId) -> int {
      int repId = (nodeId >= 0 && nodeId < (int)node_to_rep.size()) ? node_to_rep[nodeId] : nodeId;
      if (candSet.find(repId) != candSet.end()) return repId;
      if (candSet.find(nodeId) != candSet.end()) return nodeId;
      return candidates[0];
    };

    vector<int> chosen;
    chosen.reserve(k);
    std::unordered_set<int> chosenSet;
    chosenSet.reserve(k * 2);

    int rep1 = map_to_candidate(nodes[u]);
    if (chosenSet.insert(rep1).second) chosen.push_back(rep1);
    if ((int)chosen.size() < k) {
      int rep2 = map_to_candidate(nodes[v]);
      if (chosenSet.insert(rep2).second) chosen.push_back(rep2);
    }

    // Prepare min distance to chosen reps for each local node
    vector<double> minDist(m, INF);
    vector<double> dist;

    auto add_rep_and_update = [&](int repGlobalId) {
      auto it = idx.find(repGlobalId);
      if (it == idx.end()) return;
      int start = it->second;
      distances_from(start, ladj, dist);
      for (int i = 0; i < m; i++) {
        if (dist[i] >= 0.0 && dist[i] < minDist[i]) minDist[i] = dist[i];
      }
    };

    for (int r : chosen) add_rep_and_update(r);

    // Convert candidates to local indices once
    vector<int> candLocal;
    candLocal.reserve(candidates.size());
    for (int c : candidates) {
      auto it = idx.find(c);
      if (it != idx.end()) candLocal.push_back(it->second);
    }

    while ((int)chosen.size() < k) {
      int bestLocal = -1;
      double bestScore = -1.0;
      for (int li : candLocal) {
        int gid = nodes[li];
        int mapped = map_to_candidate(gid);
        if (chosenSet.find(mapped) != chosenSet.end()) continue;
        double score = minDist[li];
        if (score > bestScore) {
          bestScore = score;
          bestLocal = li;
        }
      }
      if (bestLocal < 0) break;
      int nextRep = map_to_candidate(nodes[bestLocal]);
      if (!chosenSet.insert(nextRep).second) break;
      chosen.push_back(nextRep);
      add_rep_and_update(nextRep);
    }

    std::sort(chosen.begin(), chosen.end());
    reps.push_back(std::move(chosen));
  }

  return reps;
}


