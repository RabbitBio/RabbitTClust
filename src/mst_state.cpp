#include "mst_state.h"

#include "cluster_postprocess.h"
#include "UnionFind.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <omp.h>

using std::cerr;
using std::endl;
using std::string;
using std::vector;

namespace {

// ---------------- MinHash distance (mirrors greedy's MinHash path) ----------------
inline double minhash_distance(SketchInfo& qry, SketchInfo& ref, bool is_containment) {
    if (is_containment) return qry.minHash->containDistance(ref.minHash);
    return qry.minHash->distance(ref.minHash);
}

// ---------------- KSSD distance ----------------
// Jaccard + Mash formula matches greedy::jaccard / calculate_mash_distance but is
// kept local so mst_state.cpp does not depend on the GREEDY_CLUST-gated greedy.cpp.
template <typename T>
double jaccard_generic(const std::vector<T>& a, const std::vector<T>& b) {
    uint64_t common = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) i++;
        else if (a[i] > b[j]) j++;
        else { common++; i++; j++; }
    }
    uint64_t uni = (uint64_t)a.size() + b.size() - common;
    return uni == 0 ? 0.0 : (double)common / uni;
}

template <typename T>
double mash_distance(const std::vector<T>& a, const std::vector<T>& b, double k) {
    double j = jaccard_generic(a, b);
    if (j == 1.0) return 0.0;
    if (j == 0.0) return 1.0;
    double d = -std::log(2.0 * j / (1.0 + j)) / k;
    return d > 1.0 ? 1.0 : d;
}

// I/O helpers ----------------------------------------------------------------
template <typename T>
inline void write_pod(std::ofstream& ofs, const T& v) {
    ofs.write(reinterpret_cast<const char*>(&v), sizeof(T));
}
template <typename T>
inline void read_pod(std::ifstream& ifs, T& v) {
    ifs.read(reinterpret_cast<char*>(&v), sizeof(T));
}

inline void write_string(std::ofstream& ofs, const string& s) {
    uint32_t n = (uint32_t)s.size();
    write_pod(ofs, n);
    if (n) ofs.write(s.data(), n);
}
inline void read_string(std::ifstream& ifs, string& s) {
    uint32_t n = 0;
    read_pod(ifs, n);
    s.resize(n);
    if (n) ifs.read(&s[0], n);
}

template <typename T>
inline void write_vec_pod(std::ofstream& ofs, const std::vector<T>& v) {
    uint64_t n = (uint64_t)v.size();
    write_pod(ofs, n);
    if (n) ofs.write(reinterpret_cast<const char*>(v.data()), sizeof(T) * n);
}
template <typename T>
inline void read_vec_pod(std::ifstream& ifs, std::vector<T>& v) {
    uint64_t n = 0;
    read_pod(ifs, n);
    v.resize(n);
    if (n) ifs.read(reinterpret_cast<char*>(v.data()), sizeof(T) * n);
}

constexpr const char MINHASH_MAGIC[10] = "MHMSTST01";
constexpr const char KSSD_MAGIC[10]    = "KSMSTST01";

}  // namespace

// ===========================================================================
// MinHashMstState
// ===========================================================================

void MinHashMstState::rehydrate_representatives() {
    representatives.clear();
    representatives.reserve(rep_hashes.size());
    for (size_t i = 0; i < rep_hashes.size(); i++) {
        SketchInfo s;
        s.id = (i < representative_ids.size()) ? representative_ids[i] : (int)i;
        s.fileName = (i < rep_file_names.size()) ? rep_file_names[i] : string();
        s.totalSeqLength = (i < rep_total_lens.size()) ? rep_total_lens[i] : 0;
        s.isContainment = is_containment;
        s.minHash = new Sketch::MinHash(kmer_size, sketch_size);
        s.minHash->loadMinHashes(rep_hashes[i]);
        s.KSSD = nullptr;
        s.WMinHash = nullptr;
        s.HLL = nullptr;
        s.OMH = nullptr;
        representatives.push_back(std::move(s));
    }
}

void MinHashMstState::build_inverted_index_from_reps() {
    inverted_index.clear();
    inverted_index.reserve(representatives.size() * 100);
    for (size_t i = 0; i < rep_hashes.size(); i++) {
        for (uint64_t h : rep_hashes[i]) {
            inverted_index[h].push_back((int)i);
        }
    }
}

bool MinHashMstState::save(const string& path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        cerr << "ERROR: MinHashMstState::save cannot open " << path << endl;
        return false;
    }
    ofs.write(MINHASH_MAGIC, 9);

    write_pod(ofs, threshold);
    write_pod(ofs, kmer_size);
    write_pod(ofs, sketch_size);
    write_pod(ofs, contain_compress);
    write_pod(ofs, is_containment);
    write_pod(ofs, sketch_by_file);
    write_pod(ofs, N);

    // Representatives
    uint64_t R = (uint64_t)rep_hashes.size();
    write_pod(ofs, R);
    for (uint64_t i = 0; i < R; i++) {
        int rid = (i < representative_ids.size()) ? representative_ids[i] : -1;
        write_pod(ofs, rid);
        uint64_t tlen = (i < rep_total_lens.size()) ? rep_total_lens[i] : 0;
        write_pod(ofs, tlen);
        write_string(ofs, (i < rep_file_names.size()) ? rep_file_names[i] : string());
        write_vec_pod(ofs, rep_hashes[i]);
    }

    // Clusters
    uint64_t C = (uint64_t)clusters.size();
    write_pod(ofs, C);
    for (const auto& cl : clusters) {
        write_vec_pod(ofs, cl);
    }

    // Members
    uint64_t M = (uint64_t)member_names.size();
    write_pod(ofs, M);
    for (uint64_t i = 0; i < M; i++) {
        write_string(ofs, member_names[i]);
    }
    write_vec_pod(ofs, member_lens);

    // Inverted index
    uint64_t H = (uint64_t)inverted_index.size();
    write_pod(ofs, H);
    for (const auto& kv : inverted_index) {
        write_pod(ofs, kv.first);
        write_vec_pod(ofs, kv.second);
    }

    ofs.close();
    cerr << "Saved MinHash MST state to " << path
         << " (reps=" << R << ", clusters=" << C << ", members=" << M
         << ", unique_hashes=" << H << ")" << endl;
    return true;
}

bool MinHashMstState::load(const string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        cerr << "ERROR: MinHashMstState::load cannot open " << path << endl;
        return false;
    }
    char magic[10] = {0};
    ifs.read(magic, 9);
    if (std::memcmp(magic, MINHASH_MAGIC, 9) != 0) {
        cerr << "ERROR: " << path << " is not a MinHash MST state (bad magic)" << endl;
        return false;
    }
    read_pod(ifs, threshold);
    read_pod(ifs, kmer_size);
    read_pod(ifs, sketch_size);
    read_pod(ifs, contain_compress);
    read_pod(ifs, is_containment);
    read_pod(ifs, sketch_by_file);
    read_pod(ifs, N);

    uint64_t R = 0;
    read_pod(ifs, R);
    representative_ids.resize(R);
    rep_total_lens.resize(R);
    rep_file_names.resize(R);
    rep_hashes.resize(R);
    for (uint64_t i = 0; i < R; i++) {
        read_pod(ifs, representative_ids[i]);
        read_pod(ifs, rep_total_lens[i]);
        read_string(ifs, rep_file_names[i]);
        read_vec_pod(ifs, rep_hashes[i]);
    }

    uint64_t C = 0;
    read_pod(ifs, C);
    clusters.resize(C);
    for (uint64_t i = 0; i < C; i++) read_vec_pod(ifs, clusters[i]);

    uint64_t M = 0;
    read_pod(ifs, M);
    member_names.resize(M);
    for (uint64_t i = 0; i < M; i++) read_string(ifs, member_names[i]);
    read_vec_pod(ifs, member_lens);

    uint64_t H = 0;
    read_pod(ifs, H);
    inverted_index.clear();
    inverted_index.reserve(H);
    for (uint64_t i = 0; i < H; i++) {
        uint64_t h;
        read_pod(ifs, h);
        vector<int> lst;
        read_vec_pod(ifs, lst);
        inverted_index.emplace(h, std::move(lst));
    }

    ifs.close();
    cerr << "Loaded MinHash MST state from " << path
         << " (reps=" << R << ", clusters=" << C << ", members=" << M
         << ", unique_hashes=" << H << ")" << endl;

    rehydrate_representatives();
    return true;
}

// ===========================================================================
// KssdMstState
// ===========================================================================

void KssdMstState::rehydrate_representatives() {
    size_t R = use64 ? rep_hash64.size() : rep_hash32.size();
    representatives.clear();
    representatives.reserve(R);
    for (size_t i = 0; i < R; i++) {
        KssdSketchInfo s;
        s.id = (i < representative_ids.size()) ? representative_ids[i] : (int)i;
        s.fileName = (i < rep_file_names.size()) ? rep_file_names[i] : string();
        s.totalSeqLength = (i < rep_total_lens.size()) ? rep_total_lens[i] : 0;
        s.use64 = use64;
        if (use64) {
            s.hash64_arr = rep_hash64[i];
            s.sketchsize = (uint32_t)s.hash64_arr.size();
        } else {
            s.hash32_arr = rep_hash32[i];
            s.sketchsize = (uint32_t)s.hash32_arr.size();
        }
        representatives.push_back(std::move(s));
    }
}

void KssdMstState::build_inverted_index_from_reps() {
    inverted_index_32.clear();
    inverted_index_64.clear();
    if (use64) {
        inverted_index_64.reserve(rep_hash64.size() * 100);
        for (size_t i = 0; i < rep_hash64.size(); i++) {
            for (uint64_t h : rep_hash64[i]) inverted_index_64[h].push_back((int)i);
        }
    } else {
        inverted_index_32.reserve(rep_hash32.size() * 100);
        for (size_t i = 0; i < rep_hash32.size(); i++) {
            for (uint32_t h : rep_hash32[i]) inverted_index_32[h].push_back((int)i);
        }
    }
}

bool KssdMstState::save(const string& path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        cerr << "ERROR: KssdMstState::save cannot open " << path << endl;
        return false;
    }
    ofs.write(KSSD_MAGIC, 9);

    write_pod(ofs, threshold);
    write_pod(ofs, kmer_size);
    write_pod(ofs, half_k);
    write_pod(ofs, half_subk);
    write_pod(ofs, drlevel);
    write_pod(ofs, use64);
    write_pod(ofs, sketch_by_file);
    write_pod(ofs, N);

    uint64_t R = use64 ? rep_hash64.size() : rep_hash32.size();
    write_pod(ofs, R);
    for (uint64_t i = 0; i < R; i++) {
        int rid = (i < representative_ids.size()) ? representative_ids[i] : -1;
        write_pod(ofs, rid);
        uint64_t tlen = (i < rep_total_lens.size()) ? rep_total_lens[i] : 0;
        write_pod(ofs, tlen);
        write_string(ofs, (i < rep_file_names.size()) ? rep_file_names[i] : string());
        if (use64) write_vec_pod(ofs, rep_hash64[i]);
        else       write_vec_pod(ofs, rep_hash32[i]);
    }

    uint64_t C = clusters.size();
    write_pod(ofs, C);
    for (const auto& cl : clusters) write_vec_pod(ofs, cl);

    uint64_t M = member_names.size();
    write_pod(ofs, M);
    for (uint64_t i = 0; i < M; i++) write_string(ofs, member_names[i]);
    write_vec_pod(ofs, member_lens);

    if (use64) {
        uint64_t H = inverted_index_64.size();
        write_pod(ofs, H);
        for (const auto& kv : inverted_index_64) {
            write_pod(ofs, kv.first);
            write_vec_pod(ofs, kv.second);
        }
    } else {
        uint64_t H = inverted_index_32.size();
        write_pod(ofs, H);
        for (const auto& kv : inverted_index_32) {
            write_pod(ofs, kv.first);
            write_vec_pod(ofs, kv.second);
        }
    }

    ofs.close();
    cerr << "Saved KSSD MST state to " << path
         << " (reps=" << R << ", clusters=" << C << ", members=" << M << ")" << endl;
    return true;
}

bool KssdMstState::load(const string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        cerr << "ERROR: KssdMstState::load cannot open " << path << endl;
        return false;
    }
    char magic[10] = {0};
    ifs.read(magic, 9);
    if (std::memcmp(magic, KSSD_MAGIC, 9) != 0) {
        cerr << "ERROR: " << path << " is not a KSSD MST state (bad magic)" << endl;
        return false;
    }
    read_pod(ifs, threshold);
    read_pod(ifs, kmer_size);
    read_pod(ifs, half_k);
    read_pod(ifs, half_subk);
    read_pod(ifs, drlevel);
    read_pod(ifs, use64);
    read_pod(ifs, sketch_by_file);
    read_pod(ifs, N);

    uint64_t R = 0;
    read_pod(ifs, R);
    representative_ids.resize(R);
    rep_total_lens.resize(R);
    rep_file_names.resize(R);
    if (use64) { rep_hash64.resize(R); rep_hash32.clear(); }
    else       { rep_hash32.resize(R); rep_hash64.clear(); }
    for (uint64_t i = 0; i < R; i++) {
        read_pod(ifs, representative_ids[i]);
        read_pod(ifs, rep_total_lens[i]);
        read_string(ifs, rep_file_names[i]);
        if (use64) read_vec_pod(ifs, rep_hash64[i]);
        else       read_vec_pod(ifs, rep_hash32[i]);
    }

    uint64_t C = 0;
    read_pod(ifs, C);
    clusters.resize(C);
    for (uint64_t i = 0; i < C; i++) read_vec_pod(ifs, clusters[i]);

    uint64_t M = 0;
    read_pod(ifs, M);
    member_names.resize(M);
    for (uint64_t i = 0; i < M; i++) read_string(ifs, member_names[i]);
    read_vec_pod(ifs, member_lens);

    uint64_t H = 0;
    read_pod(ifs, H);
    inverted_index_32.clear();
    inverted_index_64.clear();
    if (use64) {
        inverted_index_64.reserve(H);
        for (uint64_t i = 0; i < H; i++) {
            uint64_t h;
            read_pod(ifs, h);
            vector<int> lst;
            read_vec_pod(ifs, lst);
            inverted_index_64.emplace(h, std::move(lst));
        }
    } else {
        inverted_index_32.reserve(H);
        for (uint64_t i = 0; i < H; i++) {
            uint32_t h;
            read_pod(ifs, h);
            vector<int> lst;
            read_vec_pod(ifs, lst);
            inverted_index_32.emplace(h, std::move(lst));
        }
    }

    ifs.close();
    cerr << "Loaded KSSD MST state from " << path
         << " (reps=" << R << ", clusters=" << C << ", members=" << M << ")" << endl;

    rehydrate_representatives();
    return true;
}

// ===========================================================================
// Init builders
// ===========================================================================

MinHashMstState MinHashInitialMstState(
    const vector<SketchInfo>& sketches,
    const vector<EdgeInfo>& forest,
    const vector<vector<int>>& clusters,
    double threshold,
    int kmer_size,
    int sketch_size,
    int contain_compress,
    bool is_containment,
    bool sketch_by_file) {

    MinHashMstState st;
    st.threshold = threshold;
    st.kmer_size = kmer_size;
    st.sketch_size = sketch_size;
    st.contain_compress = contain_compress;
    st.is_containment = is_containment;
    st.sketch_by_file = sketch_by_file;
    st.N = (int)sketches.size();

    // Member metadata for later printing.
    st.member_names.reserve(st.N);
    st.member_lens.reserve(st.N);
    for (const auto& s : sketches) {
        st.member_names.push_back(s.fileName);
        st.member_lens.push_back(sketch_by_file ? s.totalSeqLength : (uint64_t)s.seqInfo.length);
    }

    // Collapse every cluster to a single tree-medoid (dedup_dist = +inf).
    vector<int> node_to_rep;
    auto cands = build_dedup_candidates_per_cluster(
        clusters, forest, sketches, sketch_by_file,
        std::numeric_limits<double>::max(), node_to_rep);

    st.representative_ids.reserve(clusters.size());
    st.rep_hashes.reserve(clusters.size());
    st.rep_file_names.reserve(clusters.size());
    st.rep_total_lens.reserve(clusters.size());
    st.clusters.reserve(clusters.size());

    for (size_t i = 0; i < clusters.size(); i++) {
        const auto& cl = clusters[i];
        if (cl.empty()) continue;
        int rep_id = cands[i].empty() ? cl[0] : cands[i][0];
        if (rep_id < 0 || rep_id >= st.N) rep_id = cl[0];
        const SketchInfo& rep_sk = sketches[rep_id];

        st.representative_ids.push_back(rep_id);
        st.rep_hashes.push_back(rep_sk.minHash ? rep_sk.minHash->storeMinHashes() : vector<uint64_t>());
        st.rep_file_names.push_back(rep_sk.fileName);
        st.rep_total_lens.push_back(sketch_by_file ? rep_sk.totalSeqLength : (uint64_t)rep_sk.seqInfo.length);
        st.clusters.push_back(cl);
    }

    st.rehydrate_representatives();
    st.build_inverted_index_from_reps();
    return st;
}

KssdMstState KssdInitialMstState(
    const vector<KssdSketchInfo>& sketches,
    const KssdParameters& params,
    const vector<EdgeInfo>& forest,
    const vector<vector<int>>& clusters,
    double threshold,
    bool sketch_by_file) {

    KssdMstState st;
    st.threshold = threshold;
    st.half_k = params.half_k;
    st.half_subk = params.half_subk;
    st.drlevel = params.drlevel;
    st.kmer_size = params.half_k * 2;
    st.sketch_by_file = sketch_by_file;
    st.N = (int)sketches.size();
    st.use64 = sketches.empty() ? false : sketches.front().use64;

    st.member_names.reserve(st.N);
    st.member_lens.reserve(st.N);
    for (const auto& s : sketches) {
        st.member_names.push_back(s.fileName);
        st.member_lens.push_back(sketch_by_file ? s.totalSeqLength : (uint64_t)s.seqInfo.length);
    }

    vector<int> node_to_rep;
    auto cands = build_dedup_candidates_per_cluster(
        clusters, forest, sketches, sketch_by_file,
        std::numeric_limits<double>::max(), node_to_rep);

    st.representative_ids.reserve(clusters.size());
    st.rep_file_names.reserve(clusters.size());
    st.rep_total_lens.reserve(clusters.size());
    if (st.use64) st.rep_hash64.reserve(clusters.size());
    else          st.rep_hash32.reserve(clusters.size());
    st.clusters.reserve(clusters.size());

    for (size_t i = 0; i < clusters.size(); i++) {
        const auto& cl = clusters[i];
        if (cl.empty()) continue;
        int rep_id = cands[i].empty() ? cl[0] : cands[i][0];
        if (rep_id < 0 || rep_id >= st.N) rep_id = cl[0];
        const KssdSketchInfo& rep_sk = sketches[rep_id];

        st.representative_ids.push_back(rep_id);
        st.rep_file_names.push_back(rep_sk.fileName);
        st.rep_total_lens.push_back(sketch_by_file ? rep_sk.totalSeqLength : (uint64_t)rep_sk.seqInfo.length);
        if (st.use64) st.rep_hash64.push_back(rep_sk.hash64_arr);
        else          st.rep_hash32.push_back(rep_sk.hash32_arr);
        st.clusters.push_back(cl);
    }

    st.rehydrate_representatives();
    st.build_inverted_index_from_reps();
    return st;
}

// ===========================================================================
// Append: merge-on-multi-match, UF over reps, lazy inverted-index compaction.
// ===========================================================================

namespace {

// Assigns a new sketch based on its matches to existing reps.
// matches: list of (rep_root_idx, dist) with dist <= threshold.
// survivor_out: rep_root_idx of the cluster that gets the new sketch.
// merged_out:  other rep roots (to be merged into survivor).
void decide_assignment(
    const vector<std::pair<int, double>>& matches,
    int& survivor_out,
    vector<int>& merged_out) {
    survivor_out = -1;
    merged_out.clear();
    if (matches.empty()) return;

    // Pick closest rep as survivor; all others (if any) will be merged into it.
    int best_idx = 0;
    for (size_t i = 1; i < matches.size(); i++) {
        if (matches[i].second < matches[best_idx].second) best_idx = (int)i;
    }
    survivor_out = matches[best_idx].first;
    for (size_t i = 0; i < matches.size(); i++) {
        if ((int)i == best_idx) continue;
        // deduplicate (same root can appear twice if two rep_idx share root already)
        if (std::find(merged_out.begin(), merged_out.end(), matches[i].first) == merged_out.end()) {
            merged_out.push_back(matches[i].first);
        }
    }
}

vector<vector<int>> collect_live_clusters(
    const vector<vector<int>>& clusters,
    UnionFind& uf) {
    vector<vector<int>> out;
    out.reserve(clusters.size());
    for (size_t i = 0; i < clusters.size(); i++) {
        if (clusters[i].empty()) continue;
        // Only keep roots; retired reps have clusters[i] emptied already, but
        // guard anyway so the caller does not see duplicates.
        if (uf.find((int)i) == (int)i) out.push_back(clusters[i]);
    }
    return out;
}

// Compact state by dropping retired reps (empty clusters) and rebuilding the
// inverted index. Ensures a clean starting point for the next append cycle.
void compact_minhash_state(MinHashMstState& state, UnionFind& uf) {
    const int R = (int)state.representatives.size();
    std::vector<int> old_to_new(R, -1);
    int new_r = 0;
    for (int i = 0; i < R; i++) {
        if (!state.clusters[i].empty() && uf.find(i) == i) {
            old_to_new[i] = new_r++;
        }
    }
    if (new_r == R) return;  // nothing retired

    std::vector<int> new_rep_ids;          new_rep_ids.reserve(new_r);
    std::vector<std::vector<uint64_t>> new_rep_hashes; new_rep_hashes.reserve(new_r);
    std::vector<std::string> new_rep_names; new_rep_names.reserve(new_r);
    std::vector<uint64_t> new_rep_lens;     new_rep_lens.reserve(new_r);
    std::vector<SketchInfo> new_reps;       new_reps.reserve(new_r);
    std::vector<std::vector<int>> new_clusters; new_clusters.reserve(new_r);

    for (int i = 0; i < R; i++) {
        if (old_to_new[i] < 0) {
            if (state.representatives[i].minHash) delete state.representatives[i].minHash;
            continue;
        }
        new_rep_ids.push_back(state.representative_ids[i]);
        new_rep_hashes.push_back(std::move(state.rep_hashes[i]));
        new_rep_names.push_back(std::move(state.rep_file_names[i]));
        new_rep_lens.push_back(state.rep_total_lens[i]);
        new_reps.push_back(std::move(state.representatives[i]));
        new_clusters.push_back(std::move(state.clusters[i]));
    }
    state.representative_ids = std::move(new_rep_ids);
    state.rep_hashes         = std::move(new_rep_hashes);
    state.rep_file_names     = std::move(new_rep_names);
    state.rep_total_lens     = std::move(new_rep_lens);
    state.representatives    = std::move(new_reps);
    state.clusters           = std::move(new_clusters);
    state.build_inverted_index_from_reps();
}

void compact_kssd_state(KssdMstState& state, UnionFind& uf) {
    const int R = (int)state.representatives.size();
    std::vector<int> old_to_new(R, -1);
    int new_r = 0;
    for (int i = 0; i < R; i++) {
        if (!state.clusters[i].empty() && uf.find(i) == i) {
            old_to_new[i] = new_r++;
        }
    }
    if (new_r == R) return;

    std::vector<int> new_rep_ids; new_rep_ids.reserve(new_r);
    std::vector<std::vector<uint32_t>> new_h32; new_h32.reserve(new_r);
    std::vector<std::vector<uint64_t>> new_h64; new_h64.reserve(new_r);
    std::vector<std::string> new_rep_names; new_rep_names.reserve(new_r);
    std::vector<uint64_t> new_rep_lens; new_rep_lens.reserve(new_r);
    std::vector<KssdSketchInfo> new_reps; new_reps.reserve(new_r);
    std::vector<std::vector<int>> new_clusters; new_clusters.reserve(new_r);

    for (int i = 0; i < R; i++) {
        if (old_to_new[i] < 0) continue;
        new_rep_ids.push_back(state.representative_ids[i]);
        if (state.use64) new_h64.push_back(std::move(state.rep_hash64[i]));
        else             new_h32.push_back(std::move(state.rep_hash32[i]));
        new_rep_names.push_back(std::move(state.rep_file_names[i]));
        new_rep_lens.push_back(state.rep_total_lens[i]);
        new_reps.push_back(std::move(state.representatives[i]));
        new_clusters.push_back(std::move(state.clusters[i]));
    }
    state.representative_ids = std::move(new_rep_ids);
    state.rep_hash32         = std::move(new_h32);
    state.rep_hash64         = std::move(new_h64);
    state.rep_file_names     = std::move(new_rep_names);
    state.rep_total_lens     = std::move(new_rep_lens);
    state.representatives    = std::move(new_reps);
    state.clusters           = std::move(new_clusters);
    state.build_inverted_index_from_reps();
}

}  // namespace

vector<vector<int>> MinHashMstAppendCluster(
    MinHashMstState& state,
    vector<SketchInfo>& new_sketches,
    int threads) {

    cerr << "\n[clust-mst append] using inverted-index state" << endl;
    cerr << "  existing members: " << state.N
         << ", existing reps: " << state.representatives.size()
         << ", new sketches: " << new_sketches.size() << endl;

    UnionFind uf((int)state.representatives.size());
    if (new_sketches.empty()) return collect_live_clusters(state.clusters, uf);

    // Greedy-style jaccard-at-threshold lower bound. Reused for every new sketch.
    // jaccard_min = e^{-d*k} / (2 - e^{-d*k}); derived from mash distance formula.
    const double exp_dk = std::exp(-state.threshold * (double)state.kmer_size);
    const double jaccard_min = exp_dk / (2.0 - exp_dk);

    uint64_t total_candidates = 0;
    uint64_t total_candidates_after_filter = 0;
    uint64_t total_distance_calcs = 0;
    int assigned_existing = 0;
    int merged_clusters = 0;
    int new_clusters = 0;

    for (size_t ni = 0; ni < new_sketches.size(); ni++) {
        SketchInfo& ns = new_sketches[ni];
        if (!ns.minHash) {
            cerr << "  [warn] new sketch " << ni << " has null minHash; skipping" << endl;
            continue;
        }
        vector<uint64_t> nhashes = ns.minHash->storeMinHashes();
        const int sizeQry = (int)nhashes.size();

        // Probe inverted index -> per-rep intersection counts. Using per-rep
        // (not per-root) counts lets the `min_common_needed` filter below be
        // exact relative to the survivor's sketch (representatives[root]),
        // mirroring the greedy path.
        phmap::flat_hash_map<int, int> hits_per_rep;
        hits_per_rep.reserve(64);
        for (uint64_t h : nhashes) {
            auto it = state.inverted_index.find(h);
            if (it == state.inverted_index.end()) continue;
            for (int rep_idx : it->second) {
                if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
                hits_per_rep[rep_idx]++;
            }
        }

        // Deduplicate by root: for each live cluster we only care about the
        // survivor's own hashes (representatives[root]). hits_per_rep[root]
        // is the exact intersection count with the survivor.
        vector<int> cand_roots;
        cand_roots.reserve(hits_per_rep.size());
        {
            phmap::flat_hash_set<int> seen_roots;
            seen_roots.reserve(hits_per_rep.size());
            for (const auto& kv : hits_per_rep) {
                int r = uf.find(kv.first);
                if (seen_roots.insert(r).second) cand_roots.push_back(r);
            }
        }
        total_candidates += cand_roots.size();

        // Apply greedy-style pre-filter using the exact common(query, survivor).
        vector<std::pair<int, int>> filtered; // (root, common)
        filtered.reserve(cand_roots.size());
        for (int r : cand_roots) {
            auto it = hits_per_rep.find(r);
            if (it == hits_per_rep.end()) continue; // survivor contributes nothing -> distance ~= 1.0
            const int common = it->second;
            const int sizeRef = (int)state.rep_hashes[r].size();
            if (sizeRef == 0) continue;

            int min_common_needed;
            if (state.is_containment) {
                const int minSz = std::min(sizeQry, sizeRef);
                min_common_needed = (int)(jaccard_min * minSz);
            } else {
                min_common_needed = (int)(jaccard_min * (sizeQry + sizeRef) / (1.0 + jaccard_min));
            }
            if (common < min_common_needed) continue;
            filtered.push_back({r, common});
        }
        total_candidates_after_filter += filtered.size();

        // Direct jaccard-from-count avoids a second sort-merge pass: `common` is
        // already the exact intersection size with representatives[root].
        vector<std::pair<int, double>> matches;
        matches.reserve(filtered.size());
        const double inv_k = 1.0 / (double)state.kmer_size;

        const int F = (int)filtered.size();
        vector<std::vector<std::pair<int, double>>> local_matches(threads);
        #pragma omp parallel num_threads(threads)
        {
            int tid = omp_get_thread_num();
            auto& lm = local_matches[tid];
            #pragma omp for schedule(dynamic, 16) reduction(+:total_distance_calcs)
            for (int ci = 0; ci < F; ci++) {
                const int r = filtered[ci].first;
                const int common = filtered[ci].second;
                const int sizeRef = (int)state.rep_hashes[r].size();
                total_distance_calcs++;

                double jaccard;
                if (state.is_containment) {
                    jaccard = (double)common / (double)std::min(sizeQry, sizeRef);
                } else {
                    const int denom = sizeQry + sizeRef - common;
                    if (denom <= 0) continue;
                    jaccard = (double)common / (double)denom;
                }
                double d;
                if (jaccard >= 1.0) d = 0.0;
                else if (jaccard <= 0.0) d = 1.0;
                else {
                    d = -std::log(2.0 * jaccard / (1.0 + jaccard)) * inv_k;
                    if (d > 1.0) d = 1.0;
                }
                if (d <= state.threshold && !std::isnan(d) && !std::isinf(d)) {
                    lm.push_back({r, d});
                }
            }
        }
        for (auto& lm : local_matches) {
            matches.insert(matches.end(), lm.begin(), lm.end());
        }

        int survivor = -1;
        vector<int> merged;
        decide_assignment(matches, survivor, merged);

        int new_node_id = state.N++;
        state.member_names.push_back(ns.fileName);
        state.member_lens.push_back(state.sketch_by_file ? ns.totalSeqLength
                                                         : (uint64_t)ns.seqInfo.length);

        if (survivor == -1) {
            // New cluster: this sketch becomes its own rep.
            int new_rep_idx = (int)state.representatives.size();
            state.representative_ids.push_back(new_node_id);
            state.rep_file_names.push_back(ns.fileName);
            state.rep_total_lens.push_back(state.member_lens.back());
            state.rep_hashes.push_back(nhashes);

            // Transfer ownership into representatives (create a standalone SketchInfo).
            SketchInfo dup;
            dup.id = new_node_id;
            dup.fileName = ns.fileName;
            dup.totalSeqLength = ns.totalSeqLength;
            dup.seqInfo = ns.seqInfo;
            dup.isContainment = state.is_containment;
            dup.minHash = new Sketch::MinHash(state.kmer_size, state.sketch_size);
            dup.minHash->loadMinHashes(nhashes);
            dup.KSSD = nullptr;
            dup.WMinHash = nullptr;
            dup.HLL = nullptr;
            dup.OMH = nullptr;
            state.representatives.push_back(std::move(dup));

            state.clusters.push_back({new_node_id});
            uf.extend(1);

            for (uint64_t h : nhashes) state.inverted_index[h].push_back(new_rep_idx);
            new_clusters++;
        } else {
            // Merge all other matched clusters into survivor, then add new node.
            // NOTE: UnionFind::merge picks the new root by rank, so we resolve
            // the post-merge root first and copy cluster data into *that* slot.
            // inverted_index entries for retired reps stay put; uf.find() at the
            // next probe transparently collapses them to the current root.
            for (int other : merged) {
                int other_root = uf.find(other);
                int surv_root = uf.find(survivor);
                if (other_root == surv_root) continue;
                uf.merge(surv_root, other_root);
                int new_root = uf.find(surv_root);
                int loser_root = (new_root == surv_root) ? other_root : surv_root;
                auto& src = state.clusters[loser_root];
                auto& dst = state.clusters[new_root];
                dst.insert(dst.end(), src.begin(), src.end());
                src.clear();
                merged_clusters++;
            }
            int final_root = uf.find(survivor);
            state.clusters[final_root].push_back(new_node_id);
            assigned_existing++;
        }

        if ((ni + 1) % 10000 == 0) {
            cerr << "  [append] processed " << (ni + 1) << " new sketches" << endl;
        }
    }

    cerr << "\n[clust-mst append summary]" << endl
         << "  assigned to existing : " << assigned_existing << endl
         << "  new clusters         : " << new_clusters << endl
         << "  cluster merges       : " << merged_clusters << endl
         << "  avg candidates/query : "
         << (new_sketches.empty() ? 0.0 : (double)total_candidates / new_sketches.size()) << endl
         << "  avg kept after filter: "
         << (new_sketches.empty() ? 0.0 : (double)total_candidates_after_filter / new_sketches.size()) << endl
         << "  distance calcs       : " << total_distance_calcs << endl;

    auto live = collect_live_clusters(state.clusters, uf);
    compact_minhash_state(state, uf);
    return live;
}

vector<vector<int>> KssdMstAppendCluster(
    KssdMstState& state,
    vector<KssdSketchInfo>& new_sketches,
    int threads) {

    cerr << "\n[clust-mst --fast append] using inverted-index state" << endl;
    cerr << "  existing members: " << state.N
         << ", existing reps: " << state.representatives.size()
         << ", new sketches: " << new_sketches.size() << endl;

    UnionFind uf((int)state.representatives.size());

    // Greedy-style pre-filter constants. `radio` is the symmetric size-ratio
    // cap beyond which the mash distance cannot be <= threshold. `jaccard_min`
    // is the minimum jaccard that can achieve the threshold distance.
    const double exp_dk = std::exp(-state.threshold * (double)state.kmer_size);
    const double jaccard_min = exp_dk / (2.0 - exp_dk);
    const double radio = std::pow(exp_dk, -1.0); // same shape as greedy's calculateMaxSizeRatio

    uint64_t total_candidates = 0;
    uint64_t total_candidates_after_filter = 0;
    uint64_t total_distance_calcs = 0;
    int assigned_existing = 0;
    int merged_clusters = 0;
    int new_clusters = 0;

    for (size_t ni = 0; ni < new_sketches.size(); ni++) {
        KssdSketchInfo& ns = new_sketches[ni];
        if (ns.use64 != state.use64) {
            cerr << "  [warn] new sketch " << ni
                 << " has use64=" << ns.use64 << " but state.use64=" << state.use64
                 << "; skipping" << endl;
            continue;
        }

        const int sizeQry = state.use64 ? (int)ns.hash64_arr.size() : (int)ns.hash32_arr.size();

        // Per-rep intersection counts (not per-root). Lets the min_common_needed
        // pre-filter be exact against the survivor sketch representatives[root].
        phmap::flat_hash_map<int, int> hits_per_rep;
        hits_per_rep.reserve(64);
        if (state.use64) {
            for (uint64_t h : ns.hash64_arr) {
                auto it = state.inverted_index_64.find(h);
                if (it == state.inverted_index_64.end()) continue;
                for (int rep_idx : it->second) {
                    if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
                    hits_per_rep[rep_idx]++;
                }
            }
        } else {
            for (uint32_t h : ns.hash32_arr) {
                auto it = state.inverted_index_32.find(h);
                if (it == state.inverted_index_32.end()) continue;
                for (int rep_idx : it->second) {
                    if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
                    hits_per_rep[rep_idx]++;
                }
            }
        }

        // Dedup to live roots. Only the survivor's own hits count for the
        // distance-to-survivor semantics.
        vector<int> cand_roots;
        cand_roots.reserve(hits_per_rep.size());
        {
            phmap::flat_hash_set<int> seen_roots;
            seen_roots.reserve(hits_per_rep.size());
            for (const auto& kv : hits_per_rep) {
                int r = uf.find(kv.first);
                if (seen_roots.insert(r).second) cand_roots.push_back(r);
            }
        }
        total_candidates += cand_roots.size();

        // Size-ratio + min-common pre-filter, mirroring the greedy KSSD path.
        vector<std::pair<int, int>> filtered; // (root, common)
        filtered.reserve(cand_roots.size());
        for (int r : cand_roots) {
            auto it = hits_per_rep.find(r);
            if (it == hits_per_rep.end()) continue;
            const int common = it->second;
            const int sizeRef = state.use64
                ? (int)state.rep_hash64[r].size()
                : (int)state.rep_hash32[r].size();
            if (sizeRef == 0) continue;

            const double ratio = (double)sizeQry / (double)sizeRef;
            if (ratio > radio || ratio < 1.0 / radio) continue;

            const int min_common_needed = (int)(jaccard_min * (sizeQry + sizeRef) / (1.0 + jaccard_min));
            if (common < min_common_needed) continue;
            filtered.push_back({r, common});
        }
        total_candidates_after_filter += filtered.size();

        // Direct jaccard-from-count: `common` is already the exact intersection
        // with representatives[root] (since hits_per_rep[r] was incremented
        // only when a query hash matched r's own sketch).
        vector<std::pair<int, double>> matches;
        matches.reserve(filtered.size());
        const double inv_k = 1.0 / (double)state.kmer_size;
        const int F = (int)filtered.size();

        vector<std::vector<std::pair<int, double>>> local_matches(threads);
        #pragma omp parallel num_threads(threads)
        {
            int tid = omp_get_thread_num();
            auto& lm = local_matches[tid];
            #pragma omp for schedule(dynamic, 16) reduction(+:total_distance_calcs)
            for (int ci = 0; ci < F; ci++) {
                const int r = filtered[ci].first;
                const int common = filtered[ci].second;
                const int sizeRef = state.use64
                    ? (int)state.rep_hash64[r].size()
                    : (int)state.rep_hash32[r].size();
                total_distance_calcs++;

                const int denom = sizeQry + sizeRef - common;
                if (denom <= 0) continue;
                const double jaccard = (double)common / (double)denom;

                double d;
                if (jaccard >= 1.0) d = 0.0;
                else if (jaccard <= 0.0) d = 1.0;
                else {
                    d = -std::log(2.0 * jaccard / (1.0 + jaccard)) * inv_k;
                    if (d > 1.0) d = 1.0;
                }
                if (d <= state.threshold && !std::isnan(d) && !std::isinf(d)) {
                    lm.push_back({r, d});
                }
            }
        }
        for (auto& lm : local_matches) {
            matches.insert(matches.end(), lm.begin(), lm.end());
        }

        int survivor = -1;
        vector<int> merged;
        decide_assignment(matches, survivor, merged);

        int new_node_id = state.N++;
        state.member_names.push_back(ns.fileName);
        state.member_lens.push_back(state.sketch_by_file ? ns.totalSeqLength
                                                         : (uint64_t)ns.seqInfo.length);

        if (survivor == -1) {
            int new_rep_idx = (int)state.representatives.size();
            state.representative_ids.push_back(new_node_id);
            state.rep_file_names.push_back(ns.fileName);
            state.rep_total_lens.push_back(state.member_lens.back());

            KssdSketchInfo dup;
            dup.id = new_node_id;
            dup.fileName = ns.fileName;
            dup.totalSeqLength = ns.totalSeqLength;
            dup.seqInfo = ns.seqInfo;
            dup.use64 = state.use64;
            if (state.use64) {
                dup.hash64_arr = ns.hash64_arr;
                dup.sketchsize = (uint32_t)dup.hash64_arr.size();
                state.rep_hash64.push_back(ns.hash64_arr);
                for (uint64_t h : ns.hash64_arr) state.inverted_index_64[h].push_back(new_rep_idx);
            } else {
                dup.hash32_arr = ns.hash32_arr;
                dup.sketchsize = (uint32_t)dup.hash32_arr.size();
                state.rep_hash32.push_back(ns.hash32_arr);
                for (uint32_t h : ns.hash32_arr) state.inverted_index_32[h].push_back(new_rep_idx);
            }
            state.representatives.push_back(std::move(dup));
            state.clusters.push_back({new_node_id});
            uf.extend(1);
            new_clusters++;
        } else {
            for (int other : merged) {
                int other_root = uf.find(other);
                int surv_root = uf.find(survivor);
                if (other_root == surv_root) continue;
                uf.merge(surv_root, other_root);
                int new_root = uf.find(surv_root);
                int loser_root = (new_root == surv_root) ? other_root : surv_root;
                auto& src = state.clusters[loser_root];
                auto& dst = state.clusters[new_root];
                dst.insert(dst.end(), src.begin(), src.end());
                src.clear();
                merged_clusters++;
            }
            int final_root = uf.find(survivor);
            state.clusters[final_root].push_back(new_node_id);
            assigned_existing++;
        }

        if ((ni + 1) % 10000 == 0) {
            cerr << "  [append] processed " << (ni + 1) << " new sketches" << endl;
        }
    }

    cerr << "\n[clust-mst --fast append summary]" << endl
         << "  assigned to existing : " << assigned_existing << endl
         << "  new clusters         : " << new_clusters << endl
         << "  cluster merges       : " << merged_clusters << endl
         << "  avg candidates/query : "
         << (new_sketches.empty() ? 0.0 : (double)total_candidates / new_sketches.size()) << endl
         << "  avg kept after filter: "
         << (new_sketches.empty() ? 0.0 : (double)total_candidates_after_filter / new_sketches.size()) << endl
         << "  distance calcs       : " << total_distance_calcs << endl;

    auto live = collect_live_clusters(state.clusters, uf);
    compact_kssd_state(state, uf);
    return live;
}

// ===========================================================================
// Output
// ===========================================================================

void printMstStateClusterResult(
    const vector<vector<int>>& clusters,
    const vector<string>& member_names,
    const vector<uint64_t>& member_lens,
    bool sketch_by_file,
    const string& output_file,
    double threshold) {
    FILE* fp = std::fopen(output_file.c_str(), "w");
    if (!fp) {
        cerr << "ERROR: printMstStateClusterResult cannot open " << output_file << endl;
        return;
    }
    if (threshold >= 0.0) {
        std::fprintf(fp, "# Clustering threshold: %.6f\n", threshold);
        std::fprintf(fp, "# Total clusters: %zu\n", clusters.size());
        std::fprintf(fp, "#\n");
    }
    for (size_t i = 0; i < clusters.size(); i++) {
        std::fprintf(fp, "the cluster %zu is: \n", i);
        for (size_t j = 0; j < clusters[i].size(); j++) {
            int id = clusters[i][j];
            const char* name = "N/A";
            uint64_t len = 0;
            if (id >= 0 && (size_t)id < member_names.size()) {
                name = member_names[id].c_str();
                len = (size_t)id < member_lens.size() ? member_lens[id] : 0;
            }
            if (sketch_by_file) {
                std::fprintf(fp, "\t%5zu\t%6d\t%12lunt\t%20s\n",
                             j, id, (unsigned long)len, name);
            } else {
                std::fprintf(fp, "\t%6zu\t%6d\t%12lunt\t%20s\n",
                             j, id, (unsigned long)len, name);
            }
        }
        std::fprintf(fp, "\n");
    }
    std::fclose(fp);
}

// ---------------------------------------------------------------------------
// MST RepDB read-only operations: query / assign / stats
// ---------------------------------------------------------------------------

namespace {

// Live-cluster-id map: rep_idx -> compact id among non-empty clusters. Clusters
// whose rep has been merged away (clusters[i].empty()) are skipped. Retired reps
// resolve to their survivor via the inverted-index (loaded state) or the
// caller-provided UnionFind during append; for load()-only state we treat every
// non-empty cluster as a root.
struct LiveClusterIndex {
    std::vector<int> rep_to_live; // rep_idx -> live id, -1 if retired
    int live_count = 0;
    int size_of(int rep_idx) const {
        if (rep_idx < 0) return 0;
        return -1; // filled by build_from(); use size_at(live_id) instead
    }
};

template <typename StateT>
LiveClusterIndex build_live_index(const StateT& state) {
    LiveClusterIndex idx;
    idx.rep_to_live.assign(state.clusters.size(), -1);
    int live = 0;
    for (size_t i = 0; i < state.clusters.size(); i++) {
        if (!state.clusters[i].empty()) idx.rep_to_live[i] = live++;
    }
    idx.live_count = live;
    return idx;
}

// Top-k helper: maintain a max-heap of size k by distance; smaller is better.
struct HitLess {
    bool operator()(const MstQueryHit& a, const MstQueryHit& b) const {
        return a.distance < b.distance; // max-heap -> keep smallest k
    }
};

template <typename StateT>
std::vector<MstQueryHit> sort_topk(std::vector<std::pair<int, double>>& cand_dists,
                                   const StateT& state,
                                   const LiveClusterIndex& live,
                                   int topk) {
    std::sort(cand_dists.begin(), cand_dists.end(),
              [](const std::pair<int,double>& a, const std::pair<int,double>& b){
                  return a.second < b.second;
              });
    std::vector<MstQueryHit> hits;
    int k = topk > 0 ? topk : (int)cand_dists.size();
    if ((int)cand_dists.size() < k) k = (int)cand_dists.size();
    hits.reserve(k);
    for (int i = 0; i < k; i++) {
        int r = cand_dists[i].first;
        MstQueryHit h;
        h.rep_idx   = r;
        h.cluster_id = (r >= 0 && r < (int)live.rep_to_live.size()) ? live.rep_to_live[r] : -1;
        h.distance  = cand_dists[i].second;
        h.rep_name  = (r >= 0 && r < (int)state.rep_file_names.size())
                          ? state.rep_file_names[r] : std::string("rep_") + std::to_string(r);
        h.cluster_size = (r >= 0 && r < (int)state.clusters.size())
                             ? (int)state.clusters[r].size() : 0;
        hits.push_back(std::move(h));
    }
    return hits;
}

} // namespace

// ---------------- MinHash query / assign --------------------------------------

std::vector<MstQueryHit> MinHashMstQueryTopK(
    const MinHashMstState& state, SketchInfo& query, int topk, int threads) {
    std::vector<MstQueryHit> empty;
    if (state.representatives.empty() || !query.minHash) return empty;

    vector<uint64_t> qhashes = query.minHash->storeMinHashes();

    // Collect candidate rep indices via the inverted index.
    phmap::flat_hash_map<int, int> hits;
    hits.reserve(64);
    for (uint64_t h : qhashes) {
        auto it = state.inverted_index.find(h);
        if (it == state.inverted_index.end()) continue;
        for (int rep_idx : it->second) {
            if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
            // For load()-only states there's no live UF: skip retired reps
            // (rep_idx whose cluster is empty).
            if (state.clusters[rep_idx].empty()) continue;
            hits[rep_idx]++;
        }
    }

    vector<int> cand;
    cand.reserve(hits.size());
    for (auto& kv : hits) cand.push_back(kv.first);

    vector<std::pair<int, double>> cand_dists(cand.size());
    #pragma omp parallel for num_threads(threads > 0 ? threads : 1) schedule(dynamic)
    for (size_t i = 0; i < cand.size(); i++) {
        int r = cand[i];
        SketchInfo& rep = const_cast<SketchInfo&>(state.representatives[r]);
        double d = std::numeric_limits<double>::infinity();
        if (rep.minHash) d = minhash_distance(query, rep, state.is_containment);
        if (std::isnan(d)) d = std::numeric_limits<double>::infinity();
        cand_dists[i] = {r, d};
    }

    LiveClusterIndex live = build_live_index(state);
    return sort_topk(cand_dists, state, live, topk);
}

MstAssignResult MinHashMstAssign(
    const MinHashMstState& state, SketchInfo& query, int threads) {
    MstAssignResult out;
    auto hits = MinHashMstQueryTopK(state, query, 1, threads);
    if (hits.empty()) return out;                                   // rep_idx=-1 (novel)
    if (hits[0].distance > state.threshold) return out;             // beyond threshold -> novel
    return hits[0];
}

// ---------------- KSSD query / assign -----------------------------------------

std::vector<MstQueryHit> KssdMstQueryTopK(
    const KssdMstState& state, const KssdSketchInfo& query, int topk, int threads) {
    std::vector<MstQueryHit> empty;
    if (state.representatives.empty()) return empty;

    phmap::flat_hash_map<int, int> hits;
    hits.reserve(64);
    if (state.use64) {
        for (uint64_t h : query.hash64_arr) {
            auto it = state.inverted_index_64.find(h);
            if (it == state.inverted_index_64.end()) continue;
            for (int rep_idx : it->second) {
                if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
                if (state.clusters[rep_idx].empty()) continue;
                hits[rep_idx]++;
            }
        }
    } else {
        for (uint32_t h : query.hash32_arr) {
            auto it = state.inverted_index_32.find(h);
            if (it == state.inverted_index_32.end()) continue;
            for (int rep_idx : it->second) {
                if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;
                if (state.clusters[rep_idx].empty()) continue;
                hits[rep_idx]++;
            }
        }
    }

    vector<int> cand;
    cand.reserve(hits.size());
    for (auto& kv : hits) cand.push_back(kv.first);

    vector<std::pair<int, double>> cand_dists(cand.size());
    #pragma omp parallel for num_threads(threads > 0 ? threads : 1) schedule(dynamic)
    for (size_t i = 0; i < cand.size(); i++) {
        int r = cand[i];
        double d;
        if (state.use64) {
            d = mash_distance<uint64_t>(state.representatives[r].hash64_arr,
                                        query.hash64_arr, (double)state.kmer_size);
        } else {
            d = mash_distance<uint32_t>(state.representatives[r].hash32_arr,
                                        query.hash32_arr, (double)state.kmer_size);
        }
        if (std::isnan(d)) d = std::numeric_limits<double>::infinity();
        cand_dists[i] = {r, d};
    }

    LiveClusterIndex live = build_live_index(state);
    return sort_topk(cand_dists, state, live, topk);
}

MstAssignResult KssdMstAssign(
    const KssdMstState& state, const KssdSketchInfo& query, int threads) {
    MstAssignResult out;
    auto hits = KssdMstQueryTopK(state, query, 1, threads);
    if (hits.empty()) return out;
    if (hits[0].distance > state.threshold) return out;
    return hits[0];
}

// ---------------- Stats -------------------------------------------------------

namespace {

template <typename StateT>
void print_cluster_size_histogram(const StateT& state, std::ostream& os) {
    // Histogram buckets: 1, 2, 3-5, 6-10, 11-100, 101-1000, >1000.
    int b1 = 0, b2 = 0, b3_5 = 0, b6_10 = 0, b11_100 = 0, b101_1000 = 0, bgt = 0;
    int live = 0;
    int total_members = 0;
    int max_size = 0;
    int min_size = std::numeric_limits<int>::max();
    for (size_t i = 0; i < state.clusters.size(); i++) {
        int sz = (int)state.clusters[i].size();
        if (sz == 0) continue;
        live++;
        total_members += sz;
        if (sz > max_size) max_size = sz;
        if (sz < min_size) min_size = sz;
        if (sz == 1) b1++;
        else if (sz == 2) b2++;
        else if (sz <= 5) b3_5++;
        else if (sz <= 10) b6_10++;
        else if (sz <= 100) b11_100++;
        else if (sz <= 1000) b101_1000++;
        else bgt++;
    }
    if (live == 0) min_size = 0;

    os << "  Live clusters:    " << live << "\n";
    os << "  Total members:    " << total_members << "\n";
    os << "  Cluster size:     min=" << min_size
       << " max=" << max_size
       << " avg=" << std::fixed << std::setprecision(2)
       << (live ? (double)total_members / live : 0.0) << "\n";
    os << "  Size histogram:\n";
    os << "    size=1         : " << b1 << "\n";
    os << "    size=2         : " << b2 << "\n";
    os << "    size=3-5       : " << b3_5 << "\n";
    os << "    size=6-10      : " << b6_10 << "\n";
    os << "    size=11-100    : " << b11_100 << "\n";
    os << "    size=101-1000  : " << b101_1000 << "\n";
    os << "    size>1000      : " << bgt << "\n";
}

} // namespace

void MinHashMstPrintStats(const MinHashMstState& state, std::ostream& os) {
    os << "========== MinHash MST RepDB stats ==========\n";
    os << "  Kmer size:        " << state.kmer_size << "\n";
    os << "  Sketch size:      " << state.sketch_size << "\n";
    os << "  Containment:      " << (state.is_containment ? "yes" : "no") << "\n";
    if (state.is_containment) {
        os << "  Contain compress: " << state.contain_compress << "\n";
    }
    os << "  Threshold:        " << std::fixed << std::setprecision(6) << state.threshold << "\n";
    os << "  Total reps slots: " << state.representatives.size() << "\n";
    os << "  sketch_by_file:   " << (state.sketch_by_file ? "yes" : "no") << "\n";
    os << "  Total members N:  " << state.N << "\n";
    os << "  Inverted index:   " << state.inverted_index.size() << " unique hashes\n";
    print_cluster_size_histogram(state, os);
    os << "==============================================\n";
}

void KssdMstPrintStats(const KssdMstState& state, std::ostream& os) {
    os << "========== KSSD MST RepDB stats ==========\n";
    os << "  Kmer size:        " << state.kmer_size << "\n";
    os << "  half_k:           " << state.half_k << "\n";
    os << "  half_subk:        " << state.half_subk << "\n";
    os << "  drlevel:          " << state.drlevel << "\n";
    os << "  use64:            " << (state.use64 ? "yes" : "no") << "\n";
    os << "  Threshold:        " << std::fixed << std::setprecision(6) << state.threshold << "\n";
    os << "  Total reps slots: " << state.representatives.size() << "\n";
    os << "  sketch_by_file:   " << (state.sketch_by_file ? "yes" : "no") << "\n";
    os << "  Total members N:  " << state.N << "\n";
    if (state.use64) {
        os << "  Inverted index:   " << state.inverted_index_64.size() << " unique hashes (64-bit)\n";
    } else {
        os << "  Inverted index:   " << state.inverted_index_32.size() << " unique hashes (32-bit)\n";
    }
    print_cluster_size_histogram(state, os);
    os << "==========================================\n";
}
