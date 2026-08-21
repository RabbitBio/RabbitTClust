#ifndef RABBITTCLUST_TREE_FORMAT_H
#define RABBITTCLUST_TREE_FORMAT_H

#include <string>
#include <vector>

#include "MSTTypes.h"

namespace tree_format {

// PHYLIP data programs use a 10-character species-name field.  Labels that
// cannot be represented losslessly are replaced by deterministic aliases;
// callers can persist encoded/original as a sidecar mapping.
struct PhylipLabelSet {
  std::vector<std::string> encoded;
  bool changed = false;
};

// Quote a taxon label for ordinary Newick output while preserving its text.
std::string quote_newick_label(const std::string& label);

// Return unique, ASCII-only PHYLIP labels no longer than 10 characters.
PhylipLabelSet make_phylip_labels(
    const std::vector<std::string>& original_labels);

// Serialize one connected MST as one Newick tree.  Every input label is
// emitted as a terminal taxon.  When an MST vertex is internal, a zero-length
// pendant edge represents that taxon at the original graph vertex, preserving
// all pairwise MST path lengths.
//
// labels must already be valid Newick tokens (quote_newick_label can be used
// for ordinary Newick; make_phylip_labels returns safe unquoted tokens).
// Throws std::invalid_argument when edges do not form a connected tree.
std::string build_mst_newick(
    const std::vector<std::string>& labels,
    const std::vector<EdgeInfo>& mst);

}  // namespace tree_format

#endif
