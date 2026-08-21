#include "TreeFormat.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace tree_format {
namespace {

constexpr std::size_t kPhylipMaxLabelLength = 10;

bool is_phylip_safe_label(const std::string& label) {
  if (label.empty() || label.size() > kPhylipMaxLabelLength) return false;

  for (unsigned char c : label) {
    // A deliberately conservative subset avoids Newick delimiters and the
    // PHYLIP convention that converts underscores into spaces.
    if (!(std::isalnum(c) || c == '.' || c == '-')) return false;
  }
  return true;
}

std::string make_phylip_alias(std::size_t value) {
  static constexpr char kBase36Digits[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  std::string alias(kPhylipMaxLabelLength, '0');
  alias[0] = 'T';
  for (std::size_t position = kPhylipMaxLabelLength;
       position > 1;
       --position) {
    alias[position - 1] = kBase36Digits[value % 36];
    value /= 36;
  }
  if (value != 0) {
    throw std::invalid_argument(
        "too many taxa to generate unique 10-character PHYLIP labels");
  }
  return alias;
}

struct DisjointSet {
  explicit DisjointSet(std::size_t size) : parent(size), rank(size, 0) {
    std::iota(parent.begin(), parent.end(), 0);
  }

  std::size_t find(std::size_t value) {
    if (parent[value] != value) parent[value] = find(parent[value]);
    return parent[value];
  }

  bool unite(std::size_t lhs, std::size_t rhs) {
    lhs = find(lhs);
    rhs = find(rhs);
    if (lhs == rhs) return false;
    if (rank[lhs] < rank[rhs]) std::swap(lhs, rhs);
    parent[rhs] = lhs;
    if (rank[lhs] == rank[rhs]) ++rank[lhs];
    return true;
  }

  std::vector<std::size_t> parent;
  std::vector<unsigned char> rank;
};

using Adjacency = std::vector<std::vector<std::pair<int, double>>>;

Adjacency validate_and_build_adjacency(
    std::size_t taxon_count,
    const std::vector<EdgeInfo>& mst) {
  if (taxon_count == 0) {
    if (!mst.empty()) {
      throw std::invalid_argument("an empty taxon set cannot contain MST edges");
    }
    return {};
  }

  const std::size_t expected_edges = taxon_count - 1;
  if (mst.size() != expected_edges) {
    std::ostringstream message;
    message << "the Newick/PHYLIP tree requires one connected MST with "
            << expected_edges << " edges for " << taxon_count
            << " taxa, but received " << mst.size()
            << "; refusing to omit disconnected components";
    throw std::invalid_argument(message.str());
  }

  Adjacency adjacency(taxon_count);
  DisjointSet components(taxon_count);
  for (std::size_t edge_index = 0; edge_index < mst.size(); ++edge_index) {
    const EdgeInfo& edge = mst[edge_index];
    if (edge.preNode < 0 || edge.sufNode < 0 ||
        static_cast<std::size_t>(edge.preNode) >= taxon_count ||
        static_cast<std::size_t>(edge.sufNode) >= taxon_count) {
      std::ostringstream message;
      message << "MST edge " << edge_index << " has an out-of-range endpoint";
      throw std::invalid_argument(message.str());
    }
    if (edge.preNode == edge.sufNode) {
      std::ostringstream message;
      message << "MST edge " << edge_index << " is a self-loop";
      throw std::invalid_argument(message.str());
    }
    if (!std::isfinite(edge.dist) || edge.dist < 0.0) {
      std::ostringstream message;
      message << "MST edge " << edge_index
              << " has a non-finite or negative branch length";
      throw std::invalid_argument(message.str());
    }
    if (!components.unite(
            static_cast<std::size_t>(edge.preNode),
            static_cast<std::size_t>(edge.sufNode))) {
      std::ostringstream message;
      message << "MST edge " << edge_index << " introduces a cycle";
      throw std::invalid_argument(message.str());
    }

    const double distance = edge.dist == 0.0 ? 0.0 : edge.dist;
    adjacency[edge.preNode].push_back({edge.sufNode, distance});
    adjacency[edge.sufNode].push_back({edge.preNode, distance});
  }

  const std::size_t root = components.find(0);
  for (std::size_t node = 1; node < taxon_count; ++node) {
    if (components.find(node) != root) {
      throw std::invalid_argument(
          "MST edges are disconnected; refusing to write a partial tree");
    }
  }
  return adjacency;
}

void append_branch_length(std::string& output, double distance) {
  char value[64];
  const int length = std::snprintf(value, sizeof(value), "%.6f", distance);
  if (length < 0 || static_cast<std::size_t>(length) >= sizeof(value)) {
    throw std::invalid_argument("cannot format an MST branch length");
  }
  output.append(value, static_cast<std::size_t>(length));
}

void append_mst_tree(
    const Adjacency& adjacency,
    const std::vector<std::string>& labels,
    std::string& output) {
  struct Frame {
    int node;
    int parent;
    std::size_t next_neighbor;
    bool opened;
    bool has_incoming_length;
    double incoming_length;
  };

  std::vector<Frame> frames;
  frames.reserve(labels.size());
  frames.push_back({0, -1, 0, false, false, 0.0});

  while (!frames.empty()) {
    Frame& frame = frames.back();
    if (!frame.opened) {
      const std::size_t child_count = adjacency[frame.node].size() -
          (frame.parent >= 0 ? 1U : 0U);
      if (child_count == 0) {
        output += labels[frame.node];
        const bool has_length = frame.has_incoming_length;
        const double length = frame.incoming_length;
        frames.pop_back();
        if (has_length) {
          output += ':';
          append_branch_length(output, length);
        }
        continue;
      }

      output += '(';
      // PHYLIP and conventional phylogenetic readers treat only leaf labels
      // as taxa.  This zero-length leaf places the sampled genome at its
      // original internal MST vertex without changing any genome-to-genome
      // path length.
      output += labels[frame.node];
      output += ":0.000000";
      frame.opened = true;
    }

    while (frame.next_neighbor < adjacency[frame.node].size() &&
           adjacency[frame.node][frame.next_neighbor].first == frame.parent) {
      ++frame.next_neighbor;
    }

    if (frame.next_neighbor == adjacency[frame.node].size()) {
      output += ')';
      const bool has_length = frame.has_incoming_length;
      const double length = frame.incoming_length;
      frames.pop_back();
      if (has_length) {
        output += ':';
        append_branch_length(output, length);
      }
      continue;
    }

    const auto child = adjacency[frame.node][frame.next_neighbor++];
    output += ',';
    frames.push_back({
        child.first, frame.node, 0, false, true, child.second});
  }
}

}  // namespace

std::string quote_newick_label(const std::string& label) {
  if (label.empty()) {
    throw std::invalid_argument("Newick taxon labels cannot be empty");
  }

  bool needs_quotes = false;
  for (unsigned char c : label) {
    if (c == '\0' || c == '\n' || c == '\r') {
      throw std::invalid_argument(
          "Newick taxon labels cannot contain NUL or line breaks");
    }
    if (std::isspace(c) || c == '(' || c == ')' || c == '[' || c == ']' ||
        c == ',' || c == ':' || c == ';' || c == '\'') {
      needs_quotes = true;
    }
  }
  if (!needs_quotes) return label;

  std::string quoted;
  quoted.reserve(label.size() + 2);
  quoted += '\'';
  for (char c : label) {
    if (c == '\'') quoted += '\'';
    quoted += c;
  }
  quoted += '\'';
  return quoted;
}

PhylipLabelSet make_phylip_labels(
    const std::vector<std::string>& original_labels) {
  PhylipLabelSet result;
  result.encoded.resize(original_labels.size());

  std::unordered_map<std::string, std::size_t> frequencies;
  frequencies.reserve(original_labels.size());
  for (const std::string& label : original_labels) ++frequencies[label];

  std::unordered_set<std::string> used;
  used.reserve(original_labels.size() * 2 + 1);
  for (const std::string& label : original_labels) {
    if (frequencies[label] == 1 && is_phylip_safe_label(label)) {
      used.insert(label);
    }
  }

  std::size_t next_alias = 1;
  for (std::size_t index = 0; index < original_labels.size(); ++index) {
    const std::string& original = original_labels[index];
    if (frequencies[original] == 1 && is_phylip_safe_label(original)) {
      result.encoded[index] = original;
      continue;
    }

    std::string alias;
    do {
      alias = make_phylip_alias(next_alias++);
    } while (used.find(alias) != used.end());
    used.insert(alias);
    result.encoded[index] = std::move(alias);
    result.changed = true;
  }
  return result;
}

std::string build_mst_newick(
    const std::vector<std::string>& labels,
    const std::vector<EdgeInfo>& mst) {
  const Adjacency adjacency =
      validate_and_build_adjacency(labels.size(), mst);
  if (labels.empty()) return ";";

  std::string output;
  std::size_t estimated_size = 1 + mst.size() * 16;
  for (const std::string& label : labels) estimated_size += label.size() + 12;
  output.reserve(estimated_size);
  append_mst_tree(adjacency, labels, output);
  output += ';';
  return output;
}

}  // namespace tree_format
