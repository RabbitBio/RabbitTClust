#include "containment_scc.h"

#include "fracminhash.h"
#include "kseq.h"
#include "phmap.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

namespace rabbittclust {
namespace {

using Clock = std::chrono::steady_clock;

struct FracNode {
    uint32_t id = 0;
    std::string file_name;
    std::string sequence_name;
    std::string comment;
    uint64_t length = 0;
    std::vector<uint64_t> hashes;
};

struct PostingIndex {
    phmap::flat_hash_map<uint64_t, uint32_t> hash_to_slot;
    std::vector<uint32_t> frequencies;
    std::vector<uint64_t> offsets;
    std::vector<uint32_t> postings;
};

struct GraphBuildStats {
    uint64_t retained_hashes = 0;
    uint64_t unique_hashes = 0;
    uint64_t posting_visits = 0;
    uint64_t exact_comparisons = 0;
    uint64_t passing_pairs = 0;
};

struct PairHit {
    uint32_t target = 0;
    uint32_t common = 0;
    float source_containment = 0.0f;
    float target_containment = 0.0f;
    bool source_to_target = false;
    bool target_to_source = false;
};

struct HashSlot {
    uint32_t frequency = 0;
    uint32_t slot = 0;
    uint64_t hash = 0;
};

struct CandidateVote {
    uint32_t node = 0;
    uint16_t hits = 0;
};

struct OutputStats {
    double sketch_seconds = 0.0;
    double graph_seconds = 0.0;
    double scc_seconds = 0.0;
    double output_seconds = 0.0;
    double total_seconds = 0.0;
    GraphBuildStats graph;
};

struct SinkProjection {
    SccPartition partition;
    // For each SCC: >=0 is its unique reachable sink SCC; -1 is ambiguous.
    std::vector<int64_t> sink_of_component;
};

double elapsedSeconds(const Clock::time_point& begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}

uint64_t mix64(uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

std::string cleanField(std::string value) {
    for (char& ch : value) {
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = ' ';
    }
    return value.empty() ? "N/A" : value;
}

bool readGenomeList(const std::string& path,
                    std::vector<std::string>& files,
                    std::string& error) {
    std::ifstream input(path);
    if (!input) {
        error = "cannot open genome list: " + path;
        return false;
    }
    std::string line;
    while (std::getline(input, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        const size_t last = line.find_last_not_of(" \t\r\n");
        files.push_back(line.substr(first, last - first + 1));
    }
    if (files.empty()) {
        error = "genome list is empty: " + path;
        return false;
    }
    return true;
}

bool sketchGenomeFile(const std::string& path,
                      uint32_t scaled,
                      int kmer_size,
                      uint64_t seed,
                      uint64_t min_length,
                      std::unique_ptr<FracNode>& output,
                      std::string& error) {
    gzFile file = gzopen(path.c_str(), "r");
    if (file == nullptr) {
        error = "cannot open FASTA: " + path;
        return false;
    }
    kseq_t* record = kseq_init(file);
    Sketch::FracMinHash sketch(scaled, kmer_size, seed);
    uint64_t total_length = 0;
    std::string first_name;
    std::string first_comment;
    int status = 0;
    while ((status = kseq_read(record)) >= 0) {
        if (first_name.empty()) {
            first_name.assign(record->name.s, record->name.l);
            if (record->comment.l > 0)
                first_comment.assign(record->comment.s, record->comment.l);
        }
        total_length += static_cast<uint64_t>(record->seq.l);
        sketch.update(record->seq.s, static_cast<uint64_t>(record->seq.l));
    }
    kseq_destroy(record);
    gzclose(file);
    if (status < -1) {
        error = "malformed FASTA/FASTQ record in: " + path;
        return false;
    }
    if (total_length < min_length) return true;

    sketch.finalize();
    output.reset(new FracNode());
    output->file_name = path;
    output->sequence_name = first_name.empty() ? path : first_name;
    output->comment = first_comment;
    output->length = total_length;
    output->hashes = sketch.hashes();
    return true;
}

bool sketchFileList(const ContainmentSccOptions& options,
                    std::vector<FracNode>& nodes,
                    std::string& error) {
    std::vector<std::string> files;
    if (!readGenomeList(options.input_file, files, error)) return false;

    std::vector<std::unique_ptr<FracNode>> results(files.size());
    std::vector<std::string> errors(files.size());
    std::atomic<uint64_t> completed(0);
    std::atomic<uint64_t> filtered(0);

    #pragma omp parallel for num_threads(options.threads) schedule(dynamic, 1)
    for (size_t i = 0; i < files.size(); ++i) {
        if (!sketchGenomeFile(files[i], options.scaled, options.kmer_size,
                              options.seed, options.min_length,
                              results[i], errors[i])) {
            // The per-index error is reported after the parallel region.
        } else if (!results[i]) {
            filtered.fetch_add(1, std::memory_order_relaxed);
        }
        const uint64_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;
        if (done % 10000 == 0 || done == files.size()) {
            #pragma omp critical(containment_progress)
            std::cerr << "-----FracMinHash sketched " << done << "/"
                      << files.size() << " genome files\n";
        }
    }

    for (size_t i = 0; i < errors.size(); ++i) {
        if (!errors[i].empty()) {
            error = errors[i];
            return false;
        }
    }
    nodes.reserve(files.size() - static_cast<size_t>(filtered.load()));
    for (std::unique_ptr<FracNode>& item : results) {
        if (!item) continue;
        item->id = static_cast<uint32_t>(nodes.size());
        nodes.push_back(std::move(*item));
    }
    std::cerr << "-----filtered " << filtered.load()
              << " genome files below --min-length\n";
    return true;
}

bool sketchSequences(const ContainmentSccOptions& options,
                     std::vector<FracNode>& nodes,
                     std::string& error) {
    gzFile file = gzopen(options.input_file.c_str(), "r");
    if (file == nullptr) {
        error = "cannot open FASTA: " + options.input_file;
        return false;
    }
    kseq_t* record = kseq_init(file);
    uint64_t filtered = 0;
    int status = 0;
    while ((status = kseq_read(record)) >= 0) {
        const uint64_t length = static_cast<uint64_t>(record->seq.l);
        if (length < options.min_length) {
            ++filtered;
            continue;
        }
        Sketch::FracMinHash sketch(options.scaled, options.kmer_size,
                                   options.seed);
        sketch.update(record->seq.s, length);
        sketch.finalize();
        FracNode node;
        node.id = static_cast<uint32_t>(nodes.size());
        node.file_name = options.input_file;
        node.sequence_name.assign(record->name.s, record->name.l);
        if (record->comment.l > 0)
            node.comment.assign(record->comment.s, record->comment.l);
        node.length = length;
        node.hashes = sketch.hashes();
        nodes.push_back(std::move(node));
    }
    kseq_destroy(record);
    gzclose(file);
    if (status < -1) {
        error = "malformed FASTA/FASTQ record in: " + options.input_file;
        return false;
    }
    std::cerr << "-----filtered " << filtered
              << " sequences below --min-length\n";
    return true;
}

template <typename T>
bool writeValue(std::ofstream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return static_cast<bool>(output);
}

template <typename T>
bool readValue(std::ifstream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    return static_cast<bool>(input);
}

bool writeString(std::ofstream& output, const std::string& value) {
    if (value.size() > std::numeric_limits<uint32_t>::max()) return false;
    const uint32_t size = static_cast<uint32_t>(value.size());
    if (!writeValue(output, size)) return false;
    output.write(value.data(), size);
    return static_cast<bool>(output);
}

bool readString(std::ifstream& input, std::string& value) {
    uint32_t size = 0;
    if (!readValue(input, size)) return false;
    value.resize(size);
    if (size > 0) input.read(&value[0], size);
    return static_cast<bool>(input);
}

bool saveSketchCache(const std::string& path,
                     const std::vector<FracNode>& nodes,
                     int kmer_size,
                     uint32_t scaled,
                     uint64_t seed,
                     bool sketch_by_file,
                     std::string& error) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        error = "cannot create sketch cache: " + path;
        return false;
    }
    const char magic[8] = {'R','T','C','F','M','H','0','1'};
    output.write(magic, sizeof(magic));
    const uint32_t version = 1;
    const uint32_t k = static_cast<uint32_t>(kmer_size);
    const uint8_t by_file = sketch_by_file ? 1 : 0;
    const uint64_t count = static_cast<uint64_t>(nodes.size());
    if (!writeValue(output, version) || !writeValue(output, k) ||
        !writeValue(output, scaled) || !writeValue(output, seed) ||
        !writeValue(output, by_file) || !writeValue(output, count)) {
        error = "failed writing sketch-cache header: " + path;
        return false;
    }
    for (const FracNode& node : nodes) {
        const uint64_t hash_count = static_cast<uint64_t>(node.hashes.size());
        if (!writeValue(output, node.length) ||
            !writeString(output, node.file_name) ||
            !writeString(output, node.sequence_name) ||
            !writeString(output, node.comment) ||
            !writeValue(output, hash_count)) {
            error = "failed writing sketch-cache metadata: " + path;
            return false;
        }
        if (hash_count > 0) {
            output.write(reinterpret_cast<const char*>(node.hashes.data()),
                         static_cast<std::streamsize>(hash_count * sizeof(uint64_t)));
        }
        if (!output) {
            error = "failed writing sketch-cache hashes: " + path;
            return false;
        }
    }
    return true;
}

bool loadSketchCache(const std::string& path,
                     std::vector<FracNode>& nodes,
                     int& kmer_size,
                     uint32_t& scaled,
                     uint64_t& seed,
                     bool& sketch_by_file,
                     std::string& error) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        error = "cannot open sketch cache: " + path;
        return false;
    }
    char magic[8] = {};
    input.read(magic, sizeof(magic));
    const char expected[8] = {'R','T','C','F','M','H','0','1'};
    if (!input || std::memcmp(magic, expected, sizeof(magic)) != 0) {
        error = "invalid FracMinHash sketch-cache magic: " + path;
        return false;
    }
    uint32_t version = 0;
    uint32_t k = 0;
    uint8_t by_file = 0;
    uint64_t count = 0;
    if (!readValue(input, version) || !readValue(input, k) ||
        !readValue(input, scaled) || !readValue(input, seed) ||
        !readValue(input, by_file) || !readValue(input, count) ||
        version != 1 || k < 1 || k > 32 || scaled == 0 ||
        count > std::numeric_limits<uint32_t>::max()) {
        error = "invalid FracMinHash sketch-cache header: " + path;
        return false;
    }
    kmer_size = static_cast<int>(k);
    sketch_by_file = by_file != 0;
    nodes.clear();
    nodes.resize(static_cast<size_t>(count));
    const uint64_t threshold = scaled == 1
        ? std::numeric_limits<uint64_t>::max()
        : std::numeric_limits<uint64_t>::max() / scaled;
    for (size_t i = 0; i < nodes.size(); ++i) {
        FracNode& node = nodes[i];
        node.id = static_cast<uint32_t>(i);
        uint64_t hash_count = 0;
        if (!readValue(input, node.length) ||
            !readString(input, node.file_name) ||
            !readString(input, node.sequence_name) ||
            !readString(input, node.comment) ||
            !readValue(input, hash_count) ||
            hash_count > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
            error = "truncated FracMinHash sketch cache: " + path;
            return false;
        }
        node.hashes.resize(static_cast<size_t>(hash_count));
        if (hash_count > 0) {
            input.read(reinterpret_cast<char*>(node.hashes.data()),
                       static_cast<std::streamsize>(hash_count * sizeof(uint64_t)));
        }
        if (!input || !std::is_sorted(node.hashes.begin(), node.hashes.end()) ||
            std::adjacent_find(node.hashes.begin(), node.hashes.end()) != node.hashes.end() ||
            (!node.hashes.empty() && node.hashes.back() > threshold)) {
            error = "invalid hash vector in FracMinHash sketch cache: " + path;
            return false;
        }
    }
    return true;
}

PostingIndex buildPostingIndex(const std::vector<FracNode>& nodes,
                               GraphBuildStats& stats) {
    PostingIndex index;
    uint64_t total_hashes = 0;
    for (const FracNode& node : nodes) total_hashes += node.hashes.size();
    stats.retained_hashes = total_hashes;

    const uint64_t reserve64 = std::min<uint64_t>(
        total_hashes / 2 + 1, UINT64_C(50000000));
    index.hash_to_slot.reserve(static_cast<size_t>(reserve64));
    index.frequencies.reserve(static_cast<size_t>(reserve64));

    uint64_t visited = 0;
    for (const FracNode& node : nodes) {
        for (uint64_t hash : node.hashes) {
            auto result = index.hash_to_slot.try_emplace(
                hash, static_cast<uint32_t>(index.frequencies.size()));
            if (result.second) index.frequencies.push_back(0);
            ++index.frequencies[result.first->second];
            ++visited;
        }
        if (visited / 50000000 != (visited - node.hashes.size()) / 50000000)
            std::cerr << "-----indexed " << visited << " FracMinHash postings\n";
    }

    stats.unique_hashes = index.frequencies.size();
    index.offsets.resize(index.frequencies.size() + 1, 0);
    for (size_t i = 0; i < index.frequencies.size(); ++i)
        index.offsets[i + 1] = index.offsets[i] + index.frequencies[i];
    index.postings.resize(static_cast<size_t>(index.offsets.back()));
    std::vector<uint64_t> cursor(index.offsets.begin(), index.offsets.end() - 1);
    for (const FracNode& node : nodes) {
        for (uint64_t hash : node.hashes) {
            const uint32_t slot = index.hash_to_slot.find(hash)->second;
            index.postings[static_cast<size_t>(cursor[slot]++)] = node.id;
        }
    }
    return index;
}

uint32_t intersectionSize(const std::vector<uint64_t>& left,
                          const std::vector<uint64_t>& right) {
    size_t i = 0;
    size_t j = 0;
    uint64_t common = 0;
    while (i < left.size() && j < right.size()) {
        if (left[i] < right[j]) ++i;
        else if (right[j] < left[i]) ++j;
        else {
            ++common;
            ++i;
            ++j;
        }
    }
    return static_cast<uint32_t>(std::min<uint64_t>(
        common, std::numeric_limits<uint32_t>::max()));
}

bool betterCandidate(const CandidateVote& left, const CandidateVote& right) {
    if (left.hits != right.hits) return left.hits > right.hits;
    return left.node < right.node;
}

bool betterPairHit(const PairHit& left, const PairHit& right) {
    const float left_primary = std::max(left.source_containment,
                                        left.target_containment);
    const float right_primary = std::max(right.source_containment,
                                         right.target_containment);
    if (left_primary != right_primary) return left_primary > right_primary;
    const float left_secondary = std::min(left.source_containment,
                                          left.target_containment);
    const float right_secondary = std::min(right.source_containment,
                                           right.target_containment);
    if (left_secondary != right_secondary)
        return left_secondary > right_secondary;
    if (left.common != right.common) return left.common > right.common;
    return left.target < right.target;
}

void appendDiversePairHits(std::vector<PairHit>& group, size_t quota,
                           uint32_t source,
                           std::vector<PairHit>& selected) {
    if (quota == 0 || group.empty()) return;
    std::sort(group.begin(), group.end(), betterPairHit);
    if (group.size() <= quota) {
        selected.insert(selected.end(), group.begin(), group.end());
        return;
    }

    // The highest-containment neighbors are useful representatives, but a
    // list containing only those neighbors can stay trapped inside a dense
    // lineage. Keep a deterministic exploratory quarter from the remaining
    // passing pairs so real cross-lineage links survive sparsification.
    const size_t exploratory = quota >= 4 ? std::max<size_t>(1, quota / 4) : 0;
    const size_t best = quota - exploratory;
    selected.insert(selected.end(), group.begin(), group.begin() + best);
    if (exploratory == 0) return;

    const auto diverse_order = [source](const PairHit& left,
                                        const PairHit& right) {
        const uint64_t left_key = mix64(
            (static_cast<uint64_t>(source) << 32) ^ left.target);
        const uint64_t right_key = mix64(
            (static_cast<uint64_t>(source) << 32) ^ right.target);
        if (left_key != right_key) return left_key < right_key;
        return left.target < right.target;
    };
    auto remainder_begin = group.begin() + best;
    std::nth_element(remainder_begin, remainder_begin + exploratory,
                     group.end(), diverse_order);
    selected.insert(selected.end(), remainder_begin,
                    remainder_begin + exploratory);
}

void trimPairHits(std::vector<PairHit>& hits, uint32_t max_neighbors,
                  uint32_t source) {
    if (max_neighbors == 0 || hits.size() <= max_neighbors) {
        std::sort(hits.begin(), hits.end(), betterPairHit);
        return;
    }
    std::vector<PairHit> mutual;
    std::vector<PairHit> one_way;
    mutual.reserve(hits.size());
    one_way.reserve(hits.size());
    for (PairHit& hit : hits) {
        if (hit.source_to_target && hit.target_to_source)
            mutual.push_back(std::move(hit));
        else
            one_way.push_back(std::move(hit));
    }
    size_t one_quota = 0;
    if (!one_way.empty()) {
        one_quota = mutual.empty()
            ? std::min<size_t>(max_neighbors, one_way.size())
            : std::min<size_t>(std::max<uint32_t>(1, max_neighbors / 4),
                               one_way.size());
    }
    size_t mutual_quota = std::min<size_t>(
        mutual.size(), static_cast<size_t>(max_neighbors) - one_quota);
    size_t remaining = static_cast<size_t>(max_neighbors) -
                       mutual_quota - one_quota;
    const size_t extra_mutual = std::min(remaining,
                                         mutual.size() - mutual_quota);
    mutual_quota += extra_mutual;
    remaining -= extra_mutual;
    one_quota += std::min(remaining, one_way.size() - one_quota);

    hits.clear();
    hits.reserve(mutual_quota + one_quota);
    appendDiversePairHits(mutual, mutual_quota, source, hits);
    appendDiversePairHits(one_way, one_quota, source, hits);
    std::sort(hits.begin(), hits.end(), betterPairHit);
}

std::vector<DirectedContainmentEdge> buildContainmentGraph(
    const std::vector<FracNode>& nodes,
    const ContainmentSccOptions& options,
    double containment_threshold,
    GraphBuildStats& stats) {
    PostingIndex index = buildPostingIndex(nodes, stats);
    std::cerr << "-----FracMinHash inverted index: " << stats.unique_hashes
              << " unique hashes, " << stats.retained_hashes
              << " postings\n";

    const int thread_count = std::max(1, options.threads);
    std::vector<std::vector<DirectedContainmentEdge>> thread_edges(
        static_cast<size_t>(thread_count));
    std::atomic<uint64_t> completed(0);
    uint64_t posting_visits = 0;
    uint64_t exact_comparisons = 0;
    uint64_t passing_pairs = 0;

    #pragma omp parallel num_threads(thread_count) reduction(+:posting_visits,exact_comparisons,passing_pairs)
    {
        const int tid = omp_get_thread_num();
        std::vector<uint32_t> marks(nodes.size(),
                                    std::numeric_limits<uint32_t>::max());
        std::vector<uint16_t> hit_counts(nodes.size(), 0);
        std::vector<uint32_t> touched;
        std::vector<HashSlot> selected;
        std::vector<CandidateVote> candidates;
        std::vector<CandidateVote> random_candidates;
        std::vector<CandidateVote> ranked_candidates;
        std::vector<uint32_t> random_candidate_ids;
        std::vector<PairHit> pair_hits;
        if (options.max_neighbors > 0) {
            thread_edges[static_cast<size_t>(tid)].reserve(
                nodes.size() * static_cast<size_t>(options.max_neighbors) /
                static_cast<size_t>(thread_count) + 1024);
        }

        #pragma omp for schedule(dynamic, 8)
        for (size_t source = 0; source < nodes.size(); ++source) {
            const FracNode& source_node = nodes[source];
            if (source_node.hashes.empty()) continue;
            selected.clear();
            selected.reserve(source_node.hashes.size());
            for (uint64_t hash : source_node.hashes) {
                const auto found = index.hash_to_slot.find(hash);
                if (found == index.hash_to_slot.end()) continue;
                const uint32_t slot = found->second;
                const uint32_t frequency = index.frequencies[slot];
                if (frequency <= 1) continue;
                if (options.max_posting > 0 && frequency > options.max_posting)
                    continue;
                selected.push_back(HashSlot{frequency, slot, hash});
            }
            const auto rarer = [](const HashSlot& left, const HashSlot& right) {
                if (left.frequency != right.frequency)
                    return left.frequency < right.frequency;
                return left.hash < right.hash;
            };
            if (options.candidate_hashes > 0 &&
                selected.size() > options.candidate_hashes) {
                // Pure rare-hash selection over-specializes candidate recall
                // in very large clonal populations: accessory/lineage hashes
                // are rare, while the core hashes needed to connect lineages
                // have longer postings. Keep half of the budget rare and
                // stratify the other half across the full frequency range.
                std::sort(selected.begin(), selected.end(), rarer);
                const size_t budget = options.candidate_hashes;
                const size_t rare_count = (budget + 1) / 2;
                const size_t stratified_count = budget - rare_count;
                std::vector<HashSlot> sampled;
                sampled.reserve(budget);
                sampled.insert(sampled.end(), selected.begin(),
                               selected.begin() + rare_count);
                const size_t span = selected.size() - rare_count;
                for (size_t sample = 0; sample < stratified_count; ++sample) {
                    // Ceiling-based quantiles include the most frequent
                    // eligible hash in the last stratum.
                    const size_t rank = rare_count +
                        ((sample + 1) * span + stratified_count - 1) /
                            stratified_count - 1;
                    sampled.push_back(selected[rank]);
                }
                selected.swap(sampled);
            }
            std::sort(selected.begin(), selected.end(), rarer);

            touched.clear();
            const uint32_t stamp = static_cast<uint32_t>(source);
            for (const HashSlot& chosen : selected) {
                const uint64_t begin = index.offsets[chosen.slot];
                const uint64_t end = index.offsets[chosen.slot + 1];
                const uint64_t frequency = end - begin;
                uint64_t step = 1;
                uint64_t offset = 0;
                if (options.posting_scan_limit > 0 &&
                    frequency > options.posting_scan_limit) {
                    step = (frequency + options.posting_scan_limit - 1) /
                           options.posting_scan_limit;
                    offset = mix64(chosen.hash ^
                                   (static_cast<uint64_t>(source) << 32)) % step;
                }
                for (uint64_t position = begin + offset;
                     position < end; position += step) {
                    ++posting_visits;
                    const uint32_t target = index.postings[position];
                    if (target == source) continue;
                    if (marks[target] != stamp) {
                        marks[target] = stamp;
                        hit_counts[target] = 1;
                        touched.push_back(target);
                    } else if (hit_counts[target] !=
                               std::numeric_limits<uint16_t>::max()) {
                        ++hit_counts[target];
                    }
                }
            }

            candidates.clear();
            candidates.reserve(touched.size());
            for (uint32_t target : touched) {
                if (hit_counts[target] >= options.min_candidate_hits)
                    candidates.push_back(CandidateVote{target, hit_counts[target]});
            }
            if (options.max_candidates > 0 &&
                candidates.size() > options.max_candidates) {
                // High-vote candidates preserve local accuracy but, in a
                // population with thousands of genomes, can all come from one
                // lineage. Reserve one quarter of the exact-comparison budget
                // for a deterministic uniform sample of recalled nodes so the
                // sparse graph can retain real cross-lineage bridges.
                const size_t budget = options.max_candidates;
                const size_t random_count = std::max<size_t>(1, budget / 4);
                const auto random_order = [source](
                    const CandidateVote& left, const CandidateVote& right) {
                    const uint64_t left_key = mix64(
                        (static_cast<uint64_t>(source) << 32) ^ left.node);
                    const uint64_t right_key = mix64(
                        (static_cast<uint64_t>(source) << 32) ^ right.node);
                    if (left_key != right_key) return left_key < right_key;
                    return left.node < right.node;
                };
                random_candidates = candidates;
                std::nth_element(random_candidates.begin(),
                                 random_candidates.begin() + random_count,
                                 random_candidates.end(), random_order);
                random_candidates.resize(random_count);
                random_candidate_ids.clear();
                random_candidate_ids.reserve(random_count);
                for (const CandidateVote& candidate : random_candidates)
                    random_candidate_ids.push_back(candidate.node);
                std::sort(random_candidate_ids.begin(),
                          random_candidate_ids.end());

                ranked_candidates = candidates;
                const size_t ranked_pool = std::min(
                    ranked_candidates.size(), budget + random_count);
                if (ranked_pool < ranked_candidates.size()) {
                    std::nth_element(ranked_candidates.begin(),
                                     ranked_candidates.begin() + ranked_pool,
                                     ranked_candidates.end(), betterCandidate);
                    ranked_candidates.resize(ranked_pool);
                }
                std::sort(ranked_candidates.begin(), ranked_candidates.end(),
                          betterCandidate);

                candidates = random_candidates;
                for (const CandidateVote& candidate : ranked_candidates) {
                    if (candidates.size() == budget) break;
                    if (!std::binary_search(random_candidate_ids.begin(),
                                            random_candidate_ids.end(),
                                            candidate.node))
                        candidates.push_back(candidate);
                }
            }
            std::sort(candidates.begin(), candidates.end(), betterCandidate);

            pair_hits.clear();
            pair_hits.reserve(candidates.size());
            for (const CandidateVote& candidate : candidates) {
                ++exact_comparisons;
                const FracNode& target_node = nodes[candidate.node];
                if (target_node.hashes.empty()) continue;
                const uint32_t common = intersectionSize(source_node.hashes,
                                                         target_node.hashes);
                if (common < options.min_shared_hashes) continue;
                const float source_containment = static_cast<float>(
                    static_cast<double>(common) / source_node.hashes.size());
                const float target_containment = static_cast<float>(
                    static_cast<double>(common) / target_node.hashes.size());
                const bool source_to_target =
                    source_containment + 1e-12 >= containment_threshold;
                const bool target_to_source =
                    target_containment + 1e-12 >= containment_threshold;
                // Each query owns only its outgoing neighborhood. Although an
                // exact comparison also reveals target->source, inserting that
                // reverse edge here would bypass the target's candidate and
                // max-neighbor budgets and turn one sampled bridge into an
                // unconditional two-way shortcut. The reverse direction is
                // retained only when the target independently recalls it.
                if (!source_to_target) continue;
                ++passing_pairs;
                pair_hits.push_back(PairHit{
                    candidate.node, common, source_containment,
                    target_containment, source_to_target, target_to_source});
            }
            trimPairHits(pair_hits, options.max_neighbors,
                         static_cast<uint32_t>(source));
            std::vector<DirectedContainmentEdge>& local_edges =
                thread_edges[static_cast<size_t>(tid)];
            for (const PairHit& hit : pair_hits) {
                if (hit.source_to_target) {
                    local_edges.push_back(DirectedContainmentEdge{
                        static_cast<uint32_t>(source), hit.target, hit.common,
                        hit.source_containment, hit.target_containment});
                }
            }

            const uint64_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            if (done % 10000 == 0 || done == nodes.size()) {
                #pragma omp critical(containment_progress)
                std::cerr << "-----containment candidates processed " << done
                          << "/" << nodes.size() << " nodes\n";
            }
        }
    }

    stats.posting_visits = posting_visits;
    stats.exact_comparisons = exact_comparisons;
    stats.passing_pairs = passing_pairs;

    size_t edge_count = 0;
    for (const auto& edges : thread_edges) edge_count += edges.size();
    std::vector<DirectedContainmentEdge> edges;
    edges.reserve(edge_count);
    for (std::vector<DirectedContainmentEdge>& local : thread_edges) {
        edges.insert(edges.end(),
                     std::make_move_iterator(local.begin()),
                     std::make_move_iterator(local.end()));
    }
    const auto edge_order = [](const DirectedContainmentEdge& left,
                               const DirectedContainmentEdge& right) {
        if (left.source != right.source) return left.source < right.source;
        if (left.target != right.target) return left.target < right.target;
        if (left.containment != right.containment)
            return left.containment > right.containment;
        return left.common_hashes > right.common_hashes;
    };
    std::sort(edges.begin(), edges.end(), edge_order);
    size_t write = 0;
    for (size_t read = 0; read < edges.size(); ++read) {
        if (write > 0 && edges[read].source == edges[write - 1].source &&
            edges[read].target == edges[write - 1].target) {
            continue;
        }
        if (write != read) edges[write] = edges[read];
        ++write;
    }
    edges.resize(write);
    edges.shrink_to_fit();
    return edges;
}

void buildCsr(size_t node_count,
              const std::vector<DirectedContainmentEdge>& edges,
              bool reverse,
              std::vector<uint64_t>& offsets,
              std::vector<uint32_t>& adjacency) {
    offsets.assign(node_count + 1, 0);
    for (const DirectedContainmentEdge& edge : edges) {
        const uint32_t source = reverse ? edge.target : edge.source;
        if (source < node_count) ++offsets[source + 1];
    }
    for (size_t i = 0; i < node_count; ++i) offsets[i + 1] += offsets[i];
    adjacency.resize(edges.size());
    std::vector<uint64_t> cursor(offsets.begin(), offsets.end() - 1);
    for (const DirectedContainmentEdge& edge : edges) {
        const uint32_t source = reverse ? edge.target : edge.source;
        const uint32_t target = reverse ? edge.source : edge.target;
        if (source < node_count && target < node_count)
            adjacency[static_cast<size_t>(cursor[source]++)] = target;
    }
}

SccPartition normalizeComponents(std::vector<std::vector<uint32_t>> raw,
                                 size_t node_count) {
    SccPartition result;
    result.component_of.assign(node_count,
                               std::numeric_limits<uint32_t>::max());
    std::vector<uint32_t> order(raw.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](uint32_t left, uint32_t right) {
        if (raw[left].size() != raw[right].size())
            return raw[left].size() > raw[right].size();
        const uint32_t left_min = *std::min_element(
            raw[left].begin(), raw[left].end());
        const uint32_t right_min = *std::min_element(
            raw[right].begin(), raw[right].end());
        return left_min < right_min;
    });
    result.components.resize(raw.size());
    for (size_t new_id = 0; new_id < order.size(); ++new_id) {
        result.components[new_id] = std::move(raw[order[new_id]]);
        std::sort(result.components[new_id].begin(),
                  result.components[new_id].end());
        for (uint32_t node : result.components[new_id])
            result.component_of[node] = static_cast<uint32_t>(new_id);
    }
    return result;
}

SccPartition computeWeakComponents(
    size_t node_count,
    const std::vector<DirectedContainmentEdge>& edges) {
    std::vector<uint32_t> parent(node_count);
    std::vector<uint8_t> rank(node_count, 0);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find_root = [&](uint32_t node, auto&& self) -> uint32_t {
        if (parent[node] != node) parent[node] = self(parent[node], self);
        return parent[node];
    };
    for (const DirectedContainmentEdge& edge : edges) {
        uint32_t left = find_root(edge.source, find_root);
        uint32_t right = find_root(edge.target, find_root);
        if (left == right) continue;
        if (rank[left] < rank[right]) std::swap(left, right);
        parent[right] = left;
        if (rank[left] == rank[right]) ++rank[left];
    }
    phmap::flat_hash_map<uint32_t, uint32_t> root_to_component;
    root_to_component.reserve(node_count);
    std::vector<std::vector<uint32_t>> raw;
    for (uint32_t node = 0; node < node_count; ++node) {
        const uint32_t root = find_root(node, find_root);
        auto inserted = root_to_component.try_emplace(
            root, static_cast<uint32_t>(raw.size()));
        if (inserted.second) raw.emplace_back();
        raw[inserted.first->second].push_back(node);
    }
    return normalizeComponents(std::move(raw), node_count);
}

SccPartition computeReciprocalComponents(
    size_t node_count,
    const std::vector<DirectedContainmentEdge>& edges,
    double containment_threshold) {
    std::vector<uint32_t> parent(node_count);
    std::vector<uint8_t> rank(node_count, 0);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find_root = [&](uint32_t node, auto&& self) -> uint32_t {
        if (parent[node] != node) parent[node] = self(parent[node], self);
        return parent[node];
    };
    for (const DirectedContainmentEdge& edge : edges) {
        // A retained directed edge already passes in the forward direction.
        // Its stored reciprocal value lets us isolate pairs that pass in both
        // directions without rebuilding or searching the edge list.
        if (static_cast<double>(edge.reciprocal_containment) + 1e-7 <
            containment_threshold) continue;
        uint32_t left = find_root(edge.source, find_root);
        uint32_t right = find_root(edge.target, find_root);
        if (left == right) continue;
        if (rank[left] < rank[right]) std::swap(left, right);
        parent[right] = left;
        if (rank[left] == rank[right]) ++rank[left];
    }
    phmap::flat_hash_map<uint32_t, uint32_t> root_to_component;
    root_to_component.reserve(node_count);
    std::vector<std::vector<uint32_t>> raw;
    for (uint32_t node = 0; node < node_count; ++node) {
        const uint32_t root = find_root(node, find_root);
        auto inserted = root_to_component.try_emplace(
            root, static_cast<uint32_t>(raw.size()));
        if (inserted.second) raw.emplace_back();
        raw[inserted.first->second].push_back(node);
    }
    return normalizeComponents(std::move(raw), node_count);
}

SccPartition computeMutualKnnComponents(
    size_t node_count,
    const std::vector<DirectedContainmentEdge>& edges) {
    std::vector<uint64_t> offsets;
    std::vector<uint32_t> adjacency;
    buildCsr(node_count, edges, false, offsets, adjacency);

    std::vector<uint32_t> parent(node_count);
    std::vector<uint8_t> rank(node_count, 0);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find_root = [&](uint32_t node, auto&& self) -> uint32_t {
        if (parent[node] != node) parent[node] = self(parent[node], self);
        return parent[node];
    };
    for (const DirectedContainmentEdge& edge : edges) {
        if (edge.source >= edge.target) continue;
        const auto reverse_begin = adjacency.begin() + offsets[edge.target];
        const auto reverse_end = adjacency.begin() + offsets[edge.target + 1];
        if (!std::binary_search(reverse_begin, reverse_end, edge.source))
            continue;
        uint32_t left = find_root(edge.source, find_root);
        uint32_t right = find_root(edge.target, find_root);
        if (left == right) continue;
        if (rank[left] < rank[right]) std::swap(left, right);
        parent[right] = left;
        if (rank[left] == rank[right]) ++rank[left];
    }

    phmap::flat_hash_map<uint32_t, uint32_t> root_to_component;
    root_to_component.reserve(node_count);
    std::vector<std::vector<uint32_t>> raw;
    for (uint32_t node = 0; node < node_count; ++node) {
        const uint32_t root = find_root(node, find_root);
        auto inserted = root_to_component.try_emplace(
            root, static_cast<uint32_t>(raw.size()));
        if (inserted.second) raw.emplace_back();
        raw[inserted.first->second].push_back(node);
    }
    return normalizeComponents(std::move(raw), node_count);
}

SinkProjection computeUniqueSinkProjection(
    const SccPartition& scc,
    const std::vector<CondensationEdge>& dag,
    size_t node_count) {
    const size_t component_count = scc.components.size();
    std::vector<std::vector<uint32_t>> outgoing(component_count);
    std::vector<uint32_t> indegree(component_count, 0);
    for (const CondensationEdge& edge : dag) {
        outgoing[edge.source_component].push_back(edge.target_component);
        ++indegree[edge.target_component];
    }
    std::vector<uint32_t> queue;
    queue.reserve(component_count);
    for (uint32_t component = 0; component < component_count; ++component)
        if (indegree[component] == 0) queue.push_back(component);
    std::vector<uint32_t> topological;
    topological.reserve(component_count);
    for (size_t head = 0; head < queue.size(); ++head) {
        const uint32_t source = queue[head];
        topological.push_back(source);
        for (uint32_t target : outgoing[source])
            if (--indegree[target] == 0) queue.push_back(target);
    }
    if (topological.size() != component_count)
        throw std::logic_error("condensation graph unexpectedly contains a cycle");

    SinkProjection projection;
    projection.sink_of_component.assign(component_count, -1);
    for (auto iterator = topological.rbegin();
         iterator != topological.rend(); ++iterator) {
        const uint32_t component = *iterator;
        if (outgoing[component].empty()) {
            projection.sink_of_component[component] = component;
            continue;
        }
        int64_t unique_sink = -2;
        bool ambiguous = false;
        for (uint32_t target : outgoing[component]) {
            const int64_t target_sink = projection.sink_of_component[target];
            if (target_sink < 0) {
                ambiguous = true;
                break;
            }
            if (unique_sink == -2) unique_sink = target_sink;
            else if (unique_sink != target_sink) {
                ambiguous = true;
                break;
            }
        }
        projection.sink_of_component[component] =
            ambiguous ? -1 : unique_sink;
    }

    phmap::flat_hash_map<uint64_t, uint32_t> label_to_component;
    label_to_component.reserve(component_count);
    std::vector<std::vector<uint32_t>> raw;
    for (uint32_t component = 0; component < component_count; ++component) {
        const int64_t sink = projection.sink_of_component[component];
        const uint64_t label = sink >= 0
            ? static_cast<uint64_t>(sink)
            : static_cast<uint64_t>(component_count) + component;
        auto inserted = label_to_component.try_emplace(
            label, static_cast<uint32_t>(raw.size()));
        if (inserted.second) raw.emplace_back();
        std::vector<uint32_t>& members = raw[inserted.first->second];
        members.insert(members.end(), scc.components[component].begin(),
                       scc.components[component].end());
    }
    projection.partition = normalizeComponents(std::move(raw), node_count);
    return projection;
}

std::vector<uint32_t> representativeNodes(
    const SccPartition& partition,
    const std::vector<FracNode>& nodes) {
    std::vector<uint32_t> representatives(partition.components.size(), 0);
    for (size_t component = 0; component < partition.components.size(); ++component) {
        uint32_t best = partition.components[component].front();
        for (uint32_t node : partition.components[component]) {
            if (nodes[node].length > nodes[best].length ||
                (nodes[node].length == nodes[best].length && node < best))
                best = node;
        }
        representatives[component] = best;
    }
    return representatives;
}

std::string componentRole(size_t component_size,
                          uint32_t in_degree,
                          uint32_t out_degree) {
    if (component_size > 1) {
        if (in_degree > 0 && out_degree > 0) return "scc_intermediate";
        if (out_degree > 0) return "scc_child";
        if (in_degree > 0) return "scc_container";
        return "scc_core";
    }
    if (in_degree == 0 && out_degree == 0) return "isolated";
    if (in_degree == 0 && out_degree > 0) return "contained_leaf";
    if (in_degree > 0 && out_degree == 0) return "container_sink";
    return "containment_intermediate";
}

bool writeClusterFile(const std::string& path,
                      const SccPartition& partition,
                      const std::vector<FracNode>& nodes,
                      bool sketch_by_file,
                      double containment_threshold,
                      std::string& error) {
    std::ofstream output(path);
    if (!output) {
        error = "cannot create cluster output: " + path;
        return false;
    }
    output << "# Directed containment threshold: " << std::setprecision(10)
           << containment_threshold << '\n';
    output << "# Total clusters: " << partition.components.size() << "\n#\n";
    for (size_t component = 0; component < partition.components.size(); ++component) {
        output << "the cluster " << component << " is: \n";
        const std::vector<uint32_t>& members = partition.components[component];
        for (size_t position = 0; position < members.size(); ++position) {
            const FracNode& node = nodes[members[position]];
            if (sketch_by_file) {
                output << '\t' << position << '\t' << node.id << '\t'
                       << node.length << "nt\t" << cleanField(node.file_name)
                       << '\t' << cleanField(node.sequence_name) << '\t'
                       << cleanField(node.comment) << '\n';
            } else {
                output << '\t' << position << '\t' << node.id << '\t'
                       << node.length << "nt\t"
                       << cleanField(node.sequence_name) << '\t'
                       << cleanField(node.comment) << '\n';
            }
        }
        output << '\n';
    }
    return static_cast<bool>(output);
}

bool writeAuxiliaryOutputs(
    const ContainmentSccOptions& options,
    const std::vector<FracNode>& nodes,
    const std::vector<DirectedContainmentEdge>& edges,
    const SccPartition& partition,
    const std::vector<CondensationEdge>& dag,
    const SccPartition& reciprocal_partition,
    const SccPartition& mutual_partition,
    const SccPartition& weak_partition,
    const SinkProjection& sink_projection,
    int kmer_size,
    uint32_t scaled,
    uint64_t seed,
    bool sketch_by_file,
    double containment_threshold,
    const Clock::time_point& output_begin,
    const Clock::time_point& total_begin,
    OutputStats& stats,
    std::string& error) {
    const std::vector<uint32_t> representatives =
        representativeNodes(partition, nodes);
    std::vector<uint32_t> in_degree(partition.components.size(), 0);
    std::vector<uint32_t> out_degree(partition.components.size(), 0);
    std::vector<uint8_t> node_has_edge(nodes.size(), 0);
    for (const CondensationEdge& edge : dag) {
        ++out_degree[edge.source_component];
        ++in_degree[edge.target_component];
    }
    for (const DirectedContainmentEdge& edge : edges) {
        node_has_edge[edge.source] = 1;
        node_has_edge[edge.target] = 1;
    }

    {
        const std::string path = options.output_file + ".scc.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "scc_id\tsize\trepresentative_node\trepresentative\tmin_length"
                  "\tmax_length\tin_degree\tout_degree\trole\n";
        for (size_t component = 0; component < partition.components.size(); ++component) {
            uint64_t min_length = std::numeric_limits<uint64_t>::max();
            uint64_t max_length = 0;
            for (uint32_t node : partition.components[component]) {
                min_length = std::min(min_length, nodes[node].length);
                max_length = std::max(max_length, nodes[node].length);
            }
            const uint32_t representative = representatives[component];
            const std::string& label = sketch_by_file
                ? nodes[representative].file_name
                : nodes[representative].sequence_name;
            output << component << '\t' << partition.components[component].size()
                   << '\t' << representative << '\t' << cleanField(label)
                   << '\t' << min_length << '\t' << max_length << '\t'
                   << in_degree[component] << '\t' << out_degree[component]
                   << '\t' << componentRole(partition.components[component].size(),
                                              in_degree[component],
                                              out_degree[component]) << '\n';
        }
    }
    {
        const std::string path = options.output_file + ".members.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "scc_id\tscc_size\tnode_id\tpath\tsequence_name\tlength"
                  "\tsketch_hashes\tcomponent_in_degree\tcomponent_out_degree\trole\n";
        for (size_t component = 0; component < partition.components.size(); ++component) {
            const std::string role = componentRole(
                partition.components[component].size(), in_degree[component],
                out_degree[component]);
            for (uint32_t node_id : partition.components[component]) {
                const FracNode& node = nodes[node_id];
                output << component << '\t' << partition.components[component].size()
                       << '\t' << node.id << '\t' << cleanField(node.file_name)
                       << '\t' << cleanField(node.sequence_name) << '\t'
                       << node.length << '\t' << node.hashes.size() << '\t'
                       << in_degree[component] << '\t' << out_degree[component]
                       << '\t' << role << '\n';
            }
        }
    }
    {
        const std::string path = options.output_file + ".dag.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "source_scc\ttarget_scc\tsource_size\ttarget_size\tsource_node"
                  "\ttarget_node\tcontainment\treciprocal_containment"
                  "\tcontainment_ani\tcommon_hashes\tsource_length\ttarget_length"
                  "\tlength_ratio\tsource_path\ttarget_path\n";
        output << std::setprecision(8);
        for (const CondensationEdge& edge : dag) {
            const FracNode& source = nodes[edge.source_node];
            const FracNode& target = nodes[edge.target_node];
            const double ani = edge.containment > 0.0f
                ? std::pow(static_cast<double>(edge.containment),
                           1.0 / static_cast<double>(kmer_size))
                : 0.0;
            const double length_ratio = target.length == 0
                ? 0.0 : static_cast<double>(source.length) / target.length;
            output << edge.source_component << '\t' << edge.target_component
                   << '\t' << partition.components[edge.source_component].size()
                   << '\t' << partition.components[edge.target_component].size()
                   << '\t' << edge.source_node << '\t' << edge.target_node
                   << '\t' << edge.containment << '\t'
                   << edge.reciprocal_containment << '\t' << ani << '\t'
                   << edge.common_hashes << '\t' << source.length << '\t'
                   << target.length << '\t' << length_ratio << '\t'
                   << cleanField(source.file_name) << '\t'
                   << cleanField(target.file_name) << '\n';
        }
    }
    if (options.write_edges) {
        const std::string path = options.output_file + ".edges.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "source_node\ttarget_node\tsource_scc\ttarget_scc\tcontainment"
                  "\treciprocal_containment\tcontainment_ani\tcommon_hashes"
                  "\tsource_path\ttarget_path\n";
        output << std::setprecision(8);
        for (const DirectedContainmentEdge& edge : edges) {
            const double ani = edge.containment > 0.0f
                ? std::pow(static_cast<double>(edge.containment),
                           1.0 / static_cast<double>(kmer_size))
                : 0.0;
            output << edge.source << '\t' << edge.target << '\t'
                   << partition.component_of[edge.source] << '\t'
                   << partition.component_of[edge.target] << '\t'
                   << edge.containment << '\t' << edge.reciprocal_containment
                   << '\t' << ani << '\t' << edge.common_hashes << '\t'
                   << cleanField(nodes[edge.source].file_name) << '\t'
                   << cleanField(nodes[edge.target].file_name) << '\n';
        }
    }
    {
        const std::string path = options.output_file + ".assignment.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "scc_id\tstatus\tunique_sink_scc\tflat_cluster\n";
        for (size_t component = 0; component < partition.components.size(); ++component) {
            const int64_t sink = sink_projection.sink_of_component[component];
            const uint32_t node = partition.components[component].front();
            output << component << '\t'
                   << (sink >= 0 ? "unique" : "ambiguous") << '\t';
            if (sink >= 0) output << sink;
            else output << "NA";
            output << '\t'
                   << sink_projection.partition.component_of[node] << '\n';
        }
    }
    {
        stats.output_seconds = elapsedSeconds(output_begin);
        stats.total_seconds = elapsedSeconds(total_begin);
        uint64_t nontrivial_sccs = 0;
        uint64_t nodes_in_nontrivial = 0;
        const uint64_t isolated_nodes = static_cast<uint64_t>(std::count(
            node_has_edge.begin(), node_has_edge.end(), uint8_t{0}));
        uint64_t contained_components = 0;
        for (size_t component = 0; component < partition.components.size(); ++component) {
            if (partition.components[component].size() > 1) {
                ++nontrivial_sccs;
                nodes_in_nontrivial += partition.components[component].size();
            }
            if (out_degree[component] > 0) ++contained_components;
        }
        const std::string path = options.output_file + ".summary.tsv";
        std::ofstream output(path);
        if (!output) { error = "cannot create: " + path; return false; }
        output << "metric\tvalue\n";
        output << "nodes\t" << nodes.size() << '\n';
        output << "kmer_size\t" << kmer_size << '\n';
        output << "scaled\t" << scaled << '\n';
        output << "seed\t" << seed << '\n';
        output << "containment_threshold\t" << std::setprecision(12)
               << containment_threshold << '\n';
        output << "retained_hashes\t" << stats.graph.retained_hashes << '\n';
        output << "unique_hashes\t" << stats.graph.unique_hashes << '\n';
        output << "posting_visits\t" << stats.graph.posting_visits << '\n';
        output << "exact_comparisons\t" << stats.graph.exact_comparisons << '\n';
        output << "passing_candidate_pairs\t" << stats.graph.passing_pairs << '\n';
        output << "directed_edges\t" << edges.size() << '\n';
        output << "sccs\t" << partition.components.size() << '\n';
        output << "reciprocal_components\t"
               << reciprocal_partition.components.size() << '\n';
        output << "mutual_knn_components\t"
               << mutual_partition.components.size() << '\n';
        output << "unique_sink_flat_clusters\t"
               << sink_projection.partition.components.size() << '\n';
        output << "weak_components\t" << weak_partition.components.size() << '\n';
        output << "nontrivial_sccs\t" << nontrivial_sccs << '\n';
        output << "nodes_in_nontrivial_sccs\t" << nodes_in_nontrivial << '\n';
        output << "condensation_dag_edges\t" << dag.size() << '\n';
        output << "components_with_outgoing_containment\t"
               << contained_components << '\n';
        output << "isolated_nodes\t" << isolated_nodes << '\n';
        output << "sketch_seconds\t" << stats.sketch_seconds << '\n';
        output << "graph_seconds\t" << stats.graph_seconds << '\n';
        output << "scc_seconds\t" << stats.scc_seconds << '\n';
        output << "output_seconds\t" << stats.output_seconds << '\n';
        output << "total_seconds\t" << stats.total_seconds << '\n';
    }
    return true;
}

} // namespace

SccPartition computeStronglyConnectedComponents(
    size_t node_count,
    const std::vector<DirectedContainmentEdge>& edges) {
    SccPartition result;
    result.component_of.assign(node_count,
                               std::numeric_limits<uint32_t>::max());
    if (node_count == 0) return result;

    std::vector<uint64_t> offsets;
    std::vector<uint32_t> adjacency;
    buildCsr(node_count, edges, false, offsets, adjacency);
    std::vector<uint8_t> visited(node_count, 0);
    std::vector<uint32_t> finish_order;
    finish_order.reserve(node_count);
    std::vector<std::pair<uint32_t, uint64_t>> dfs;

    for (uint32_t start = 0; start < node_count; ++start) {
        if (visited[start]) continue;
        visited[start] = 1;
        dfs.push_back(std::make_pair(start, offsets[start]));
        while (!dfs.empty()) {
            std::pair<uint32_t, uint64_t>& frame = dfs.back();
            const uint32_t node = frame.first;
            if (frame.second < offsets[node + 1]) {
                const uint32_t next = adjacency[static_cast<size_t>(frame.second++)];
                if (!visited[next]) {
                    visited[next] = 1;
                    dfs.push_back(std::make_pair(next, offsets[next]));
                }
            } else {
                finish_order.push_back(node);
                dfs.pop_back();
            }
        }
    }

    buildCsr(node_count, edges, true, offsets, adjacency);
    std::vector<std::vector<uint32_t>> raw_components;
    std::vector<uint32_t> stack;
    for (auto iterator = finish_order.rbegin();
         iterator != finish_order.rend(); ++iterator) {
        const uint32_t start = *iterator;
        if (result.component_of[start] !=
            std::numeric_limits<uint32_t>::max()) continue;
        const uint32_t raw_id = static_cast<uint32_t>(raw_components.size());
        raw_components.emplace_back();
        result.component_of[start] = raw_id;
        stack.push_back(start);
        while (!stack.empty()) {
            const uint32_t node = stack.back();
            stack.pop_back();
            raw_components.back().push_back(node);
            for (uint64_t position = offsets[node];
                 position < offsets[node + 1]; ++position) {
                const uint32_t next = adjacency[static_cast<size_t>(position)];
                if (result.component_of[next] ==
                    std::numeric_limits<uint32_t>::max()) {
                    result.component_of[next] = raw_id;
                    stack.push_back(next);
                }
            }
        }
    }

    std::vector<uint32_t> component_order(raw_components.size());
    std::iota(component_order.begin(), component_order.end(), 0);
    std::sort(component_order.begin(), component_order.end(),
              [&](uint32_t left, uint32_t right) {
        if (raw_components[left].size() != raw_components[right].size())
            return raw_components[left].size() > raw_components[right].size();
        const uint32_t left_min = *std::min_element(
            raw_components[left].begin(), raw_components[left].end());
        const uint32_t right_min = *std::min_element(
            raw_components[right].begin(), raw_components[right].end());
        return left_min < right_min;
    });
    std::vector<uint32_t> remap(raw_components.size());
    result.components.resize(raw_components.size());
    for (size_t new_id = 0; new_id < component_order.size(); ++new_id) {
        const uint32_t raw_id = component_order[new_id];
        remap[raw_id] = static_cast<uint32_t>(new_id);
        result.components[new_id] = std::move(raw_components[raw_id]);
        std::sort(result.components[new_id].begin(),
                  result.components[new_id].end());
    }
    for (uint32_t& component : result.component_of) component = remap[component];
    return result;
}

std::vector<CondensationEdge> buildCondensationDag(
    const SccPartition& partition,
    const std::vector<DirectedContainmentEdge>& edges) {
    std::vector<CondensationEdge> dag;
    dag.reserve(edges.size());
    for (const DirectedContainmentEdge& edge : edges) {
        const uint32_t source_component = partition.component_of[edge.source];
        const uint32_t target_component = partition.component_of[edge.target];
        if (source_component == target_component) continue;
        dag.push_back(CondensationEdge{
            source_component, target_component, edge.source, edge.target,
            edge.common_hashes, edge.containment,
            edge.reciprocal_containment});
    }
    const auto order = [](const CondensationEdge& left,
                          const CondensationEdge& right) {
        if (left.source_component != right.source_component)
            return left.source_component < right.source_component;
        if (left.target_component != right.target_component)
            return left.target_component < right.target_component;
        if (left.containment != right.containment)
            return left.containment > right.containment;
        if (left.common_hashes != right.common_hashes)
            return left.common_hashes > right.common_hashes;
        if (left.source_node != right.source_node)
            return left.source_node < right.source_node;
        return left.target_node < right.target_node;
    };
    std::sort(dag.begin(), dag.end(), order);
    size_t write = 0;
    for (size_t read = 0; read < dag.size(); ++read) {
        if (write > 0 &&
            dag[read].source_component == dag[write - 1].source_component &&
            dag[read].target_component == dag[write - 1].target_component) {
            continue;
        }
        if (write != read) dag[write] = dag[read];
        ++write;
    }
    dag.resize(write);
    dag.shrink_to_fit();
    return dag;
}

int runContainmentScc(const ContainmentSccOptions& options) {
    const Clock::time_point total_begin = Clock::now();
    ContainmentSccOptions effective = options;
    std::vector<FracNode> nodes;
    bool sketch_by_file = options.sketch_by_file;
    std::string error;
    OutputStats stats;

    const Clock::time_point sketch_begin = Clock::now();
    if (!options.load_sketch_file.empty()) {
        if (!loadSketchCache(options.load_sketch_file, nodes,
                             effective.kmer_size, effective.scaled,
                             effective.seed, sketch_by_file, error)) {
            std::cerr << "ERROR: " << error << '\n';
            return 1;
        }
        std::cerr << "-----loaded " << nodes.size()
                  << " FracMinHash sketches from "
                  << options.load_sketch_file << '\n';
    } else {
        const bool ok = sketch_by_file
            ? sketchFileList(effective, nodes, error)
            : sketchSequences(effective, nodes, error);
        if (!ok) {
            std::cerr << "ERROR: " << error << '\n';
            return 1;
        }
        if (nodes.empty()) {
            std::cerr << "ERROR: no input passed --min-length\n";
            return 1;
        }
        std::string cache_path = options.save_sketch_file;
        if (cache_path.empty() && !options.no_save)
            cache_path = options.output_file + ".fmh";
        if (!cache_path.empty()) {
            if (!saveSketchCache(cache_path, nodes, effective.kmer_size,
                                 effective.scaled, effective.seed,
                                 sketch_by_file, error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            std::cerr << "-----saved FracMinHash sketch cache to: "
                      << cache_path << '\n';
        }
    }
    stats.sketch_seconds = elapsedSeconds(sketch_begin);
    std::cerr << "========time of computing/loading FracMinHash sketches is: "
              << stats.sketch_seconds << "========\n";

    const double containment_threshold =
        effective.containment_threshold >= 0.0
        ? effective.containment_threshold
        : std::exp(-static_cast<double>(effective.kmer_size) *
                   effective.distance_threshold);
    if (containment_threshold < 0.0 || containment_threshold > 1.0) {
        std::cerr << "ERROR: effective containment threshold is outside [0,1]\n";
        return 1;
    }
    std::cerr << "-----FracMinHash parameters: k=" << effective.kmer_size
              << ", scaled=" << effective.scaled
              << ", seed=" << effective.seed << '\n';
    std::cerr << "-----directed containment threshold: "
              << std::setprecision(10) << containment_threshold
              << " (containment ANI "
              << std::pow(containment_threshold,
                          1.0 / effective.kmer_size) << ")\n";

    const Clock::time_point graph_begin = Clock::now();
    std::vector<DirectedContainmentEdge> edges;
    try {
        edges = buildContainmentGraph(nodes, effective,
                                      containment_threshold, stats.graph);
    } catch (const std::exception& exception) {
        std::cerr << "ERROR: containment graph construction failed: "
                  << exception.what() << '\n';
        return 1;
    }
    stats.graph_seconds = elapsedSeconds(graph_begin);
    std::cerr << "-----retained directed containment edges: "
              << edges.size() << '\n';
    std::cerr << "========time of containment graph is: "
              << stats.graph_seconds << "========\n";

    const Clock::time_point scc_begin = Clock::now();
    SccPartition partition = computeStronglyConnectedComponents(
        nodes.size(), edges);
    std::vector<CondensationEdge> dag = buildCondensationDag(partition, edges);
    SccPartition reciprocal_partition = computeReciprocalComponents(
        nodes.size(), edges, containment_threshold);
    SccPartition mutual_partition = computeMutualKnnComponents(
        nodes.size(), edges);
    SccPartition weak_partition = computeWeakComponents(nodes.size(), edges);
    SinkProjection sink_projection = computeUniqueSinkProjection(
        partition, dag, nodes.size());
    stats.scc_seconds = elapsedSeconds(scc_begin);
    std::cerr << "-----strongly connected components: "
              << partition.components.size() << '\n';
    std::cerr << "-----condensation DAG edges: " << dag.size() << '\n';
    std::cerr << "-----reciprocal-core clusters: "
              << reciprocal_partition.components.size() << '\n';
    std::cerr << "-----mutual-kNN clusters: "
              << mutual_partition.components.size() << '\n';
    std::cerr << "-----unique-sink flat clusters: "
              << sink_projection.partition.components.size() << '\n';
    std::cerr << "-----weak containment components: "
              << weak_partition.components.size() << '\n';
    std::cerr << "========time of SCC and condensation DAG is: "
              << stats.scc_seconds << "========\n";

    const Clock::time_point output_begin = Clock::now();
    if (!writeClusterFile(options.output_file, partition, nodes,
                          sketch_by_file, containment_threshold, error)) {
        std::cerr << "ERROR: " << error << '\n';
        return 1;
    }
    if (!writeClusterFile(options.output_file + ".flat.cluster",
                          sink_projection.partition, nodes, sketch_by_file,
                          containment_threshold, error) ||
        !writeClusterFile(options.output_file + ".reciprocal.cluster",
                          reciprocal_partition, nodes, sketch_by_file,
                          containment_threshold, error) ||
        !writeClusterFile(options.output_file + ".mutual.cluster",
                          mutual_partition, nodes, sketch_by_file,
                          containment_threshold, error) ||
        !writeClusterFile(options.output_file + ".weak.cluster",
                          weak_partition, nodes, sketch_by_file,
                          containment_threshold, error)) {
        std::cerr << "ERROR: " << error << '\n';
        return 1;
    }
    if (!writeAuxiliaryOutputs(effective, nodes, edges, partition, dag,
                               reciprocal_partition, mutual_partition,
                               weak_partition,
                               sink_projection,
                               effective.kmer_size, effective.scaled,
                               effective.seed, sketch_by_file,
                               containment_threshold, output_begin,
                               total_begin, stats, error)) {
        std::cerr << "ERROR: " << error << '\n';
        return 1;
    }
    stats.output_seconds = elapsedSeconds(output_begin);
    stats.total_seconds = elapsedSeconds(total_begin);
    std::cerr << "-----wrote SCC clusters to: " << options.output_file << '\n';
    std::cerr << "-----wrote component table to: "
              << options.output_file << ".scc.tsv\n";
    std::cerr << "-----wrote containment DAG to: "
              << options.output_file << ".dag.tsv\n";
    std::cerr << "========total time of clust-containment is: "
              << stats.total_seconds << "========\n";
    return 0;
}

} // namespace rabbittclust
