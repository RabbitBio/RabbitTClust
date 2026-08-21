#ifndef RABBITTCLUST_CONTAINMENT_SCC_H
#define RABBITTCLUST_CONTAINMENT_SCC_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rabbittclust {

struct ContainmentSccOptions {
    std::string input_file;
    std::string output_file;
    std::string load_sketch_file;
    std::string save_sketch_file;
    bool sketch_by_file = false;
    bool no_save = false;
    bool write_edges = false;
    int threads = 1;
    uint64_t min_length = 10000;
    int kmer_size = 19;
    uint32_t scaled = 10000;
    uint64_t seed = 42;
    double distance_threshold = 0.02;
    double containment_threshold = -1.0;
    uint32_t candidate_hashes = 128;
    uint32_t max_candidates = 256;
    uint32_t max_neighbors = 64;
    uint32_t max_posting = 100000;
    uint32_t posting_scan_limit = 1024;
    uint32_t min_candidate_hits = 1;
    uint32_t min_shared_hashes = 3;
};

struct DirectedContainmentEdge {
    uint32_t source = 0;
    uint32_t target = 0;
    uint32_t common_hashes = 0;
    float containment = 0.0f;
    float reciprocal_containment = 0.0f;
};

struct SccPartition {
    std::vector<uint32_t> component_of;
    std::vector<std::vector<uint32_t>> components;
};

struct CondensationEdge {
    uint32_t source_component = 0;
    uint32_t target_component = 0;
    uint32_t source_node = 0;
    uint32_t target_node = 0;
    uint32_t common_hashes = 0;
    float containment = 0.0f;
    float reciprocal_containment = 0.0f;
};

SccPartition computeStronglyConnectedComponents(
    size_t node_count,
    const std::vector<DirectedContainmentEdge>& edges);

std::vector<CondensationEdge> buildCondensationDag(
    const SccPartition& partition,
    const std::vector<DirectedContainmentEdge>& edges);

int runContainmentScc(const ContainmentSccOptions& options);

} // namespace rabbittclust

#endif // RABBITTCLUST_CONTAINMENT_SCC_H
