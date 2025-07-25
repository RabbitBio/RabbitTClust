#include <iostream>
#include "SketchInfo.h"
#include "MST.h"
#include "Sketch_IO.h"
#include "common.hpp"
#include "MST_IO.h"
#include "greedy.h"
#include <vector>
#include <string>
using namespace std;

void compute_sketches(vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, bool isSave, int threads);

void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, bool no_dense, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads);

void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool sketchByFile, bool no_dense, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads);

bool tune_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, bool& isJaccard, int& kmerSize, double& threshold, int& containCompress, int& sketchSize);
bool tune_kssd_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, int& kmerSize, double& threshold, int &drlevel);

void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads);

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads);
void clust_from_mst_fast(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads);

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool no_dense, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);
void append_clust_mst_fast(string folder_path, string input_file, string output_file, bool is_newick_tree, bool no_dense, bool sketch_by_file, bool isContainment, int min_len, bool no_save, double threshold, int threads);

void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);

void compute_kssd_sketches(vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads);
void compute_kssd_clusters(vector<KssdSketchInfo>& sketches, const KssdParameters info, bool sketchByFile, bool no_dense, bool isContainment, const string folder_path, string outputFile, bool is_newick_tree, double threshold, bool isSave, int threads);

void clust_from_genome_fast(const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads);
void clust_from_sketch_fast(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, bool isContainment, double threshold, int threads);
void clust_from_genomes_fast_MPI(int my_rank, int comm_sz, const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads);
void compute_kssd_sketches_mpi(int my_rank, vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, bool isSave, int threads); 
void distribute_compute_clusters(int my_rank, int comm_sz, vector<KssdSketchInfo>& sketches, const KssdParameters info,int start_index, int end_index, bool sketchByFile, string output_file, bool is_newick_tree, string folder_path, double threshold, bool isSave, int threads, bool no_dense, bool isContainment);
