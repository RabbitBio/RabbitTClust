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

void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads);

void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool sketchByFile, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool isSave, int threads);

bool tune_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, bool& isJaccard, int& kmerSize, double& threshold, int& containCompress, int& sketchSize);

void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, double threshold, int threads);

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, double threshold, int threads);

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);

void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads);

