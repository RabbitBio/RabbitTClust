#ifndef H_MST_STATE
#define H_MST_STATE

#include "SketchInfo.h"
#include "MST.h"
#include "phmap.h"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

// Persistent clustering state for clust-mst "--save-rep" / "--append" workflow.
//
// On an initial clust-mst run with --save-rep the pipeline collapses every cluster
// produced by the MST-cut forest to a single tree-medoid representative (ties broken
// by longer genome length, matching build_dedup_candidates_per_cluster). The set of
// representatives plus a hash -> rep_idx inverted index is serialised into
// folder_path + "/mst_cluster_state.bin".
//
// On a subsequent clust-mst --append the state is loaded and each new sketch is
// classified by probing the inverted index, computing true distance to candidate
// reps, then: (a) assign to the sole matching cluster, (b) merge all matching
// clusters (MST-consistent: the new sketch transitively bridges them by edges
// <= threshold) when more than one cluster matches, or (c) create a new cluster
// whose representative is the new sketch when nothing matches.
//
// Representatives are self-contained: rep sketches are reconstructed from serialised
// hashes with Sketch::MinHash::loadMinHashes / KssdSketchInfo::hash*_arr, so append
// does not need the original pre-sketch folder at all.

struct MinHashMstState {
    // ---- Parameters (must match pre-sketch params of any appended genomes) ----
    double threshold = 0.0;
    int kmer_size = 0;
    int sketch_size = 0;
    int contain_compress = 0;     // meaningful iff is_containment
    bool is_containment = false;

    // ---- Book-keeping ----
    int N = 0;                    // total node count (original + appended so far)
    bool sketch_by_file = true;

    // ---- Representatives (one per live cluster) ----
    // rep_idx i <-> clusters[i].
    // Retired reps (merged away) keep their slot but have clusters[i].empty() and
    // are transparently redirected by UnionFind.find() during probing.
    std::vector<int> representative_ids;                 // rep idx -> original node id
    std::vector<std::vector<uint64_t>> rep_hashes;       // rep idx -> MinHash hash list
    std::vector<std::string> rep_file_names;             // rep idx -> file/sequence name
    std::vector<uint64_t> rep_total_lens;                // rep idx -> genome length

    // Reconstructed in-memory after load/init (not serialised).
    std::vector<SketchInfo> representatives;

    // ---- Cluster membership (original node ids) ----
    std::vector<std::vector<int>> clusters;

    // ---- Metadata for printing the .cluster file on append ----
    std::vector<std::string> member_names;               // size N
    std::vector<uint64_t> member_lens;                   // size N

    // ---- Inverted index: hash -> rep_idx ----
    phmap::flat_hash_map<uint64_t, std::vector<int>> inverted_index;

    bool save(const std::string& path) const;
    bool load(const std::string& path);

    // Rebuild in-memory Sketch::MinHash* objects inside `representatives` from
    // rep_hashes / rep_file_names. Safe to call after load().
    void rehydrate_representatives();

    // Rebuild inverted_index from current reps' stored hashes.
    void build_inverted_index_from_reps();
};

struct KssdMstState {
    // ---- Parameters ----
    double threshold = 0.0;
    int kmer_size = 0;
    int half_k = 0;
    int half_subk = 0;
    int drlevel = 0;
    bool use64 = false;

    int N = 0;
    bool sketch_by_file = true;

    // ---- Representatives ----
    std::vector<int> representative_ids;
    std::vector<std::vector<uint32_t>> rep_hash32;       // populated iff !use64
    std::vector<std::vector<uint64_t>> rep_hash64;       // populated iff use64
    std::vector<std::string> rep_file_names;
    std::vector<uint64_t> rep_total_lens;

    // Reconstructed after load/init.
    std::vector<KssdSketchInfo> representatives;

    std::vector<std::vector<int>> clusters;

    std::vector<std::string> member_names;
    std::vector<uint64_t> member_lens;

    phmap::flat_hash_map<uint32_t, std::vector<int>> inverted_index_32;
    phmap::flat_hash_map<uint64_t, std::vector<int>> inverted_index_64;

    bool save(const std::string& path) const;
    bool load(const std::string& path);

    // Build KssdSketchInfo mirrors of rep_hash{32,64}.
    void rehydrate_representatives();

    void build_inverted_index_from_reps();
};

// ---------------------------------------------------------------------------
// Initial-state builders. Called from the MST branch of clust-mst right after
// `tmpClust = generateClusterWithBfs(forest, N)` has been computed.
// Each cluster yields exactly one tree-medoid representative
// (via build_dedup_candidates_per_cluster with dedup_dist = +inf).
// ---------------------------------------------------------------------------
MinHashMstState MinHashInitialMstState(
    const std::vector<SketchInfo>& sketches,
    const std::vector<EdgeInfo>& forest,
    const std::vector<std::vector<int>>& clusters,
    double threshold,
    int kmer_size,
    int sketch_size,
    int contain_compress,
    bool is_containment,
    bool sketch_by_file);

KssdMstState KssdInitialMstState(
    const std::vector<KssdSketchInfo>& sketches,
    const KssdParameters& params,
    const std::vector<EdgeInfo>& forest,
    const std::vector<std::vector<int>>& clusters,
    double threshold,
    bool sketch_by_file);

// ---------------------------------------------------------------------------
// Append: classify + merge-on-multi-match. Mutates `state` in place. Returns
// live clusters (empty slots removed) suitable for printMstStateClusterResult.
// ---------------------------------------------------------------------------
std::vector<std::vector<int>> MinHashMstAppendCluster(
    MinHashMstState& state,
    std::vector<SketchInfo>& new_sketches,
    int threads);

std::vector<std::vector<int>> KssdMstAppendCluster(
    KssdMstState& state,
    std::vector<KssdSketchInfo>& new_sketches,
    int threads);

// ---------------------------------------------------------------------------
// Lightweight printer that uses only metadata persisted in the state. Matches
// the columnar format of printResult / printKssdResult (threshold header +
// "the cluster N is:" sections) so downstream scripts keep working.
// ---------------------------------------------------------------------------
void printMstStateClusterResult(
    const std::vector<std::vector<int>>& clusters,
    const std::vector<std::string>& member_names,
    const std::vector<uint64_t>& member_lens,
    bool sketch_by_file,
    const std::string& output_file,
    double threshold);

// ---------------------------------------------------------------------------
// MST RepDB read-only operations (query / assign / stats)
// ---------------------------------------------------------------------------
//
// These mirror the greedy RepDB API (repdb_query / repdb_assign / repdb_stats
// in sub_command.cpp) but operate on the MST-specific state structs.
//
// query (top-k): return the k representative with the smallest distance to the
// query sketch. No threshold filter — result is always sorted by distance.
//
// assign: return the single closest representative whose distance is <= the
// state's persisted threshold, or rep_idx = -1 if no match is within threshold.

struct MstQueryHit {
    int    rep_idx       = -1;    // index into state.representatives (UF-resolved root)
    int    cluster_id    = -1;    // compact live-cluster id (post UF compaction)
    double distance      = 0.0;
    std::string rep_name;         // genome/file name of the representative
    int    cluster_size  = 0;     // number of members in the representative's live cluster
};

using MstAssignResult = MstQueryHit; // identical layout; rep_idx == -1 means unassigned / novel

// Compute top-k nearest representatives for a single MinHash query.
// `state` must have representatives rehydrated (load() / rehydrate_representatives()).
std::vector<MstQueryHit> MinHashMstQueryTopK(
    const MinHashMstState& state,
    SketchInfo& query,
    int topk,
    int threads);

MstAssignResult MinHashMstAssign(
    const MinHashMstState& state,
    SketchInfo& query,
    int threads);

// Analogous KSSD variants.
std::vector<MstQueryHit> KssdMstQueryTopK(
    const KssdMstState& state,
    const KssdSketchInfo& query,
    int topk,
    int threads);

MstAssignResult KssdMstAssign(
    const KssdMstState& state,
    const KssdSketchInfo& query,
    int threads);

// Stats: format a human-readable report of the state (parameters, #reps, #live
// clusters, cluster-size histogram, inverted-index stats). Mirrors
// KssdClusterState::print_stats / MinHashClusterState::print_stats on the
// greedy side.
void MinHashMstPrintStats(const MinHashMstState& state, std::ostream& os);
void KssdMstPrintStats(const KssdMstState& state, std::ostream& os);

#endif // H_MST_STATE
