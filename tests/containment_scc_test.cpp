#include "containment_scc.h"

#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using rabbittclust::CondensationEdge;
using rabbittclust::DirectedContainmentEdge;
using rabbittclust::SccPartition;

static DirectedContainmentEdge edge(uint32_t source, uint32_t target) {
    DirectedContainmentEdge result;
    result.source = source;
    result.target = target;
    result.common_hashes = 10;
    result.containment = 0.95f;
    result.reciprocal_containment = 0.5f;
    return result;
}

int main() {
    // 0 and 1 are reciprocal; 2 is contained in that SCC; 3 is contained in
    // 2; node 4 is isolated. Direction is child/fragment -> container.
    const std::vector<DirectedContainmentEdge> edges = {
        edge(0, 1), edge(1, 0), edge(2, 0), edge(3, 2)
    };
    const SccPartition partition =
        rabbittclust::computeStronglyConnectedComponents(5, edges);

    assert(partition.components.size() == 4);
    assert(partition.components[0].size() == 2);
    assert(partition.component_of[0] == partition.component_of[1]);
    assert(partition.component_of[2] != partition.component_of[0]);
    assert(partition.component_of[3] != partition.component_of[2]);
    assert(partition.component_of[4] != partition.component_of[0]);

    const std::vector<CondensationEdge> dag =
        rabbittclust::buildCondensationDag(partition, edges);
    assert(dag.size() == 2);
    std::set<std::pair<uint32_t, uint32_t>> arcs;
    for (const CondensationEdge& item : dag)
        arcs.emplace(item.source_component, item.target_component);
    assert(arcs.count({partition.component_of[2], partition.component_of[0]}) == 1);
    assert(arcs.count({partition.component_of[3], partition.component_of[2]}) == 1);

    // Duplicate graph edges must collapse to one condensation arc, retaining
    // the best containment witness.
    std::vector<DirectedContainmentEdge> duplicate_edges = edges;
    DirectedContainmentEdge weaker = edge(2, 1);
    weaker.containment = 0.80f;
    duplicate_edges.push_back(weaker);
    const std::vector<CondensationEdge> duplicate_dag =
        rabbittclust::buildCondensationDag(partition, duplicate_edges);
    assert(duplicate_dag.size() == 2);

    // End-to-end fixture: identical full genomes form an SCC, a strict
    // fragment points to that SCC, and an unrelated sequence stays isolated.
    const std::string data_dir = RABBITTCLUST_TEST_DATA_DIR;
    const std::string list_path = "containment-scc-test-input.list";
    const std::string output_path = "containment-scc-test-output.cluster";
    {
        std::ofstream list(list_path);
        assert(list);
        list << data_dir << "/full_a.fa\n"
             << data_dir << "/full_b.fa\n"
             << data_dir << "/fragment.fa\n"
             << data_dir << "/unrelated.fa\n";
    }
    rabbittclust::ContainmentSccOptions options;
    options.input_file = list_path;
    options.output_file = output_path;
    options.sketch_by_file = true;
    options.no_save = true;
    options.threads = 2;
    options.min_length = 50;
    options.kmer_size = 15;
    options.scaled = 1;
    options.containment_threshold = 0.8;
    options.candidate_hashes = 0;
    options.max_candidates = 0;
    options.max_neighbors = 0;
    options.max_posting = 0;
    options.posting_scan_limit = 0;
    options.min_shared_hashes = 3;
    assert(rabbittclust::runContainmentScc(options) == 0);

    std::map<std::string, std::string> summary;
    {
        std::ifstream input(output_path + ".summary.tsv");
        assert(input);
        std::string line;
        std::getline(input, line);
        while (std::getline(input, line)) {
            const size_t tab = line.find('\t');
            assert(tab != std::string::npos);
            summary[line.substr(0, tab)] = line.substr(tab + 1);
        }
    }
    assert(summary["nodes"] == "4");
    assert(summary["sccs"] == "3");
    assert(summary["condensation_dag_edges"] == "1");
    assert(summary["unique_sink_flat_clusters"] == "2");
    assert(summary["weak_components"] == "2");

    {
        std::ifstream input(output_path + ".dag.tsv");
        assert(input);
        std::string header;
        std::string relation;
        assert(static_cast<bool>(std::getline(input, header)));
        assert(static_cast<bool>(std::getline(input, relation)));
        assert(relation.find("\t1\t0.373") != std::string::npos);
        assert(!std::getline(input, relation));
    }

    std::cout << "containment_scc_test: OK\n";
    return 0;
}
