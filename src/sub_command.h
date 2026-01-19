#include <iostream>
#include "SketchInfo.h"
#include "MST.h"
#include "Sketch_IO.h"
#include "common.hpp"
#include "MST_IO.h"
#include "greedy.h"
#include "leiden.h"
#include <vector>
#include <string>
using namespace std;

void compute_sketches(vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, bool isSave, int threads);

void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads, bool use_inverted_index);

void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool sketchByFile, bool no_dense, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads, bool use_inverted_index);

bool tune_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, bool& isJaccard, int& kmerSize, double& threshold, int& containCompress, int& sketchSize);
bool tune_kssd_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, int& kmerSize, double& threshold, int &drlevel);

void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads, bool use_inverted_index);

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, double threshold, int threads);
void clust_from_mst_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, double threshold, int threads);

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);
void append_clust_mst_fast(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, bool isContainment, int min_len, bool no_save, double threshold, int threads);

void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);
void append_clust_greedy_fast(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);

void compute_kssd_sketches(vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads);
void compute_kssd_clusters(vector<KssdSketchInfo>& sketches, const KssdParameters info, bool sketchByFile, bool no_dense, bool isContainment, const string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, double threshold, bool isSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index);

void clust_from_genome_fast(const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index);
void clust_from_sketch_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool isContainment, double threshold, int threads, double dedup_dist, int reps_per_cluster, bool use_inverted_index, bool save_rep_index);

// Build a reusable KSSD sketch + inverted index "database" into db_folder and exit.
// input_file can be either:
// - a plain genome file list (-l), one path per line
// - a RabbitTClust cluster output (e.g. *.cluster / *.cluster.dedup), from which genome paths will be extracted.
void build_kssd_db_fast(const string input_file, const string db_folder, bool isSetKmer, bool& isContainment, int minLen, int& kmerSize, int& drlevel, int threads);

// Graph-based clustering functions (Leiden/Louvain)
void clust_from_genome_leiden(const string inputFile, string outputFile, string folder_path, bool sketchByFile, const int kmerSize, const int drlevel, const int minLen, bool noSave, double threshold, double resolution, bool use_leiden, bool use_parallel_louvain, bool use_edge_parallel, bool use_warm_start, int knn_k, int threads);
void clust_from_sketch_leiden(string folder_path, string outputFile, double threshold, double resolution, bool use_leiden, bool use_parallel_louvain, bool use_edge_parallel, bool use_warm_start, int knn_k, int threads);
void clust_from_pregraph_leiden(string folder_path, string outputFile, double resolution, bool use_leiden, int threads);

