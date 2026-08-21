#include "TreeFormat.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

template <typename Callback>
void expect_invalid_argument(Callback callback) {
  bool threw = false;
  try {
    callback();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  require(threw, "expected std::invalid_argument");
}

void test_internal_mst_vertices_become_terminal_taxa() {
  const std::vector<std::string> labels = {"A", "B", "C"};
  const std::vector<EdgeInfo> edges = {
      {0, 1, 0.1},
      {1, 2, 0.2},
  };

  const std::string newick = tree_format::build_mst_newick(labels, edges);
  require(
      newick == "(A:0.000000,(B:0.000000,C:0.200000):0.100000);",
      "internal MST vertex was not emitted as a zero-length leaf taxon");
}

void test_root_leaf_prevents_a_phylip_unifurcation() {
  const std::vector<std::string> labels = {"A", "B"};
  const std::vector<EdgeInfo> edges = {{0, 1, 0.25}};

  const std::string newick = tree_format::build_mst_newick(labels, edges);
  require(
      newick == "(A:0.000000,B:0.250000);",
      "two-taxon MST produced a root unifurcation");
}

void test_newick_label_quoting() {
  require(
      tree_format::quote_newick_label("plain.fa") == "plain.fa",
      "plain Newick label changed");
  require(
      tree_format::quote_newick_label("A B") == "'A B'",
      "Newick label containing a space was not quoted");
  require(
      tree_format::quote_newick_label("O'Brien") == "'O''Brien'",
      "single quote in Newick label was not escaped");
  expect_invalid_argument([] {
    tree_format::quote_newick_label("line\nbreak");
  });
}

void test_phylip_labels_are_unique_and_bounded() {
  const std::vector<std::string> original = {
      "duplicate",
      "short.fa",
      "with_underscore",
      "/a/very/long/common/path/genome.fa",
      "duplicate",
  };
  const tree_format::PhylipLabelSet labels =
      tree_format::make_phylip_labels(original);

  require(labels.changed, "unsafe PHYLIP labels were not aliased");
  require(
      labels.encoded.size() == original.size(),
      "PHYLIP label count changed");
  require(
      labels.encoded[1] == "short.fa",
      "safe PHYLIP label changed");

  std::unordered_set<std::string> unique;
  for (const std::string& label : labels.encoded) {
    require(!label.empty(), "empty PHYLIP label generated");
    require(label.size() <= 10, "PHYLIP label exceeds 10 characters");
    require(
        label.find('_') == std::string::npos,
        "PHYLIP alias contains an underscore");
    unique.insert(label);
  }
  require(
      unique.size() == labels.encoded.size(),
      "PHYLIP labels are not unique");
}

void test_invalid_or_disconnected_graphs_are_rejected() {
  const std::vector<std::string> four_labels = {"A", "B", "C", "D"};
  expect_invalid_argument([&] {
    tree_format::build_mst_newick(
        four_labels, {{0, 1, 0.1}, {1, 2, 0.2}});
  });

  const std::vector<std::string> three_labels = {"A", "B", "C"};
  expect_invalid_argument([&] {
    tree_format::build_mst_newick(
        three_labels, {{0, 1, 0.1}, {1, 0, 0.2}});
  });
  expect_invalid_argument([&] {
    tree_format::build_mst_newick(
        three_labels, {{0, 1, -0.1}, {1, 2, 0.2}});
  });
  expect_invalid_argument([&] {
    tree_format::build_mst_newick(
        three_labels, {{0, 1, 0.1}, {1, 3, 0.2}});
  });
}

}  // namespace

int main() {
  test_internal_mst_vertices_become_terminal_taxa();
  test_root_leaf_prevents_a_phylip_unifurcation();
  test_newick_label_quoting();
  test_phylip_labels_are_unique_and_bounded();
  test_invalid_or_disconnected_graphs_are_rejected();
  std::cout << "MST tree-format tests passed\n";
  return 0;
}
