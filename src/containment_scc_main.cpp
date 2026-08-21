#include "CLI11.hpp"
#include "containment_scc.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <sys/sysinfo.h>

int main(int argc, char** argv) {
    using rabbittclust::ContainmentSccOptions;

    CLI::App app{
        "clust-containment v2.5.0, directed FracMinHash containment graph, "
        "strongly connected components, and condensation DAG"};

    ContainmentSccOptions options;
    options.threads = std::max(1, get_nprocs_conf());

    app.add_option("-i,--input", options.input_file,
                   "FASTA/FASTA.GZ input, or a genome-file list with -l");
    app.add_flag("-l,--list", options.sketch_by_file,
                 "input is a genome list, one FASTA/FASTA.GZ path per line");
    app.add_option("-o,--output", options.output_file,
                   "RabbitTClust-compatible SCC cluster output prefix/file")
        ->required();
    app.add_option("-t,--threads", options.threads,
                   "worker threads (default: all configured CPUs)");
    app.add_option("-m,--min-length", options.min_length,
                   "ignore genomes/sequences shorter than this many nt");
    app.add_option("-k,--kmer-size", options.kmer_size,
                   "FracMinHash canonical k-mer size (1..32)");
    app.add_option("--scaled", options.scaled,
                   "FracMinHash scaled factor; retain about 1/scaled k-mers");
    app.add_option("--seed", options.seed,
                   "coordinated FracMinHash seed");
    app.add_option("-d,--distance", options.distance_threshold,
                   "Mash-distance-equivalent cutoff; containment cutoff is exp(-k*d) (RefSeq assembly default: 0.02)");
    app.add_option("--containment-threshold", options.containment_threshold,
                   "direct containment cutoff in [0,1], overriding --distance");

    app.add_option("--candidate-hashes", options.candidate_hashes,
                   "frequency-stratified hashes used for candidate recall (0: all)");
    app.add_option("--max-candidates", options.max_candidates,
                   "exact candidates evaluated per node (0: all recalled)");
    app.add_option("--max-neighbors", options.max_neighbors,
                   "passing candidate pairs retained per query node (0: all)");
    app.add_option("--max-posting", options.max_posting,
                   "ignore hashes occurring in more nodes than this (0: unlimited)");
    app.add_option("--posting-scan-limit", options.posting_scan_limit,
                   "deterministically sample at most this many nodes per posting (0: all)");
    app.add_option("--min-candidate-hits", options.min_candidate_hits,
                   "minimum selected-hash votes before exact comparison");
    app.add_option("--min-shared-hashes", options.min_shared_hashes,
                   "minimum full-sketch intersection for a directed edge");

    app.add_option("--load-sketch,--presketched", options.load_sketch_file,
                   "load a native .fmh sketch cache instead of reading FASTA");
    app.add_option("--save-sketch", options.save_sketch_file,
                   "write a reusable native .fmh sketch cache");
    app.add_flag("-e,--no-save", options.no_save,
                 "do not write an implicit <output>.fmh sketch cache");
    app.add_flag("--write-edges", options.write_edges,
                 "also write every retained directed edge to <output>.edges.tsv");

    CLI11_PARSE(app, argc, argv);

    if (options.input_file.empty() && options.load_sketch_file.empty()) {
        std::cerr << "ERROR: provide -i/--input or --load-sketch\n";
        return 1;
    }
    if (!options.input_file.empty() && !options.load_sketch_file.empty()) {
        std::cerr << "ERROR: -i/--input and --load-sketch are mutually exclusive\n";
        return 1;
    }
    if (options.threads < 1) {
        std::cerr << "ERROR: --threads must be positive\n";
        return 1;
    }
    if (options.kmer_size < 1 || options.kmer_size > 32) {
        std::cerr << "ERROR: --kmer-size must be in [1,32]\n";
        return 1;
    }
    if (options.scaled == 0) {
        std::cerr << "ERROR: --scaled must be positive\n";
        return 1;
    }
    if (options.distance_threshold < 0.0 || options.distance_threshold > 1.0) {
        std::cerr << "ERROR: --distance must be in [0,1]\n";
        return 1;
    }
    if (options.containment_threshold > 1.0) {
        std::cerr << "ERROR: --containment-threshold must be in [0,1]\n";
        return 1;
    }
    if (options.min_candidate_hits == 0 || options.min_shared_hashes == 0) {
        std::cerr << "ERROR: minimum hit counts must be positive\n";
        return 1;
    }

    return rabbittclust::runContainmentScc(options);
}
