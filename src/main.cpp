/*
 * This version is for result checking using single thread.
 * The input parameter and the final output is the same with the final version.
 *
 * The input can be a single file with numbers of sequences to create sketches by sequences.
 * And it can alse be a single file with list of genome files to create sketches by files(genomes).
 * This is strategy of clustering based distance computing by sequences or by genomes.
 *
 * The program includes several sections:
 * section 1: read the input arguments and init the parameters.
 * section 2: read genome files and create the sketches.
 * section 3: compute the distance matrix and generate the Minimum Spanning Tree(MST) or greedy incremental clustering.
 * section 4: generate the clusters with the MST using different distance threshold.
 *
 * Author: Xiaoming Xu
 * Mar 5, 2021
 *
 */
#include <iostream>
#include "SketchInfo.h"
#include "Sketch.h"// need to add the include path in Makefile.
#include <zlib.h>
#include "MST.h"
#include <omp.h>
#include "UnionFind.h"
#include <algorithm>
#include "common.hpp"
#include "MST_IO.h"
#include <math.h>
#include "Sketch_IO.h"

#ifdef GREEDY_CLUST
#include "greedy.h"
#endif

#include <fstream>
#include <sstream>
#include <sys/sysinfo.h>
#include <sys/stat.h>

#include "CLI11.hpp"
#include "sub_command.h"

#ifdef GREEDY_CLUST
#else
#endif


using namespace std;

int main(int argc, char * argv[]){
	#ifdef GREEDY_CLUST
		CLI::App app{"clust-greedy v.2.2.1, greedy incremental clustering module for RabbitTClust"};
	#else
		CLI::App app{"clust-mst v.2.2.1, minimum-spanning-tree-based module for RabbitTClust"};
	#endif
	//section 1: init parameters
	int argIndex = 1;
	string inputFile = "genome.fna";
	string inputFile1 = "genome.info";
	string sketchFunc = "MinHash";
	string outputFile = "result.out";
	int threads = 1;
	threads = get_nprocs_conf();
	bool sketchByFile = false;
	bool isContainment = false;
	bool isJaccard = false;
	bool useFile = false;
	double threshold = 0.05;
	int kmerSize = 21;
	int sketchSize = 1000;
	int containCompress = 1000;
	bool mstLoadSketch = false;
	string mstSketchFile = "sketch.info";
	bool isSetKmer = false;
	uint64_t minLen = 10000;
	string folder_path;
	bool is_newick_tree = false;

	bool noSave = false;

	auto option_threads = app.add_option("-t, --threads", threads,  "set the thread number, default all CPUs of the platform");
	auto option_min_len = app.add_option("-m, --min-length", minLen, "set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000");

	auto option_containment = app.add_option("-c, --containment", containCompress, "use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress");
	auto option_kmer_size = app.add_option("-k, --kmer-size", kmerSize, "set the kmer size");
	auto option_sketch_size = app.add_option("-s, --sketch-size", sketchSize, "set the sketch size for Jaccard Index and Mash distance, default 1000");

	auto flag_input_list = app.add_flag("-l, --list", sketchByFile, "input is genome list, one genome per line");
	auto flag_no_save = app.add_flag("-e, --no-save", noSave, "not save the intermediate files, such as sketches or MST");
	auto option_threshold = app.add_option("-d, --threshold", threshold, "set the distance threshold for clustering");
	auto option_function = app.add_option("-F, --function", sketchFunc, "set the sketch function, such as MinHash, KSSD, default MinHash");
	auto option_output = app.add_option("-o, --output", outputFile, "set the output name of cluster result");
	auto option_input = app.add_option("-i, --input", inputFile, "set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)");
	auto option_presketched = app.add_option("--presketched", folder_path, "clustering by the pre-generated sketch files rather than genomes");
#ifndef GREEDY_CLUST
	auto option_premsted = app.add_option("--premsted", folder_path, "clustering by the pre-generated mst files rather than genomes for clust-mst");
	auto flag_newick_tree = app.add_flag("--newick-tree", is_newick_tree, "output the newick tree format file for clust-mst");
#endif
	auto option_append = app.add_option("--append", inputFile, "append genome file or file list with the pre-generated sketch or MST files");

	option_output->required();
	option_append->excludes(option_input);

	CLI11_PARSE(app, argc, argv);

	if(*option_function){
		fprintf(stderr, "-----set the sketch function: %s\n", sketchFunc.c_str());
	}
	if(threads < 1){
		fprintf(stderr, "-----Invalid thread number %d\n", threads);
		return 1;
	}
	if(option_threads){
		fprintf(stderr, "-----set the thread number %d\n", threads);
	}
	if(*option_min_len){
		fprintf(stderr, "-----set the filter minimum length: %ld\n", minLen);
	}
	if(*option_containment){
		isContainment = true;
		fprintf(stderr, "-----use AAF distance with containment coefficient, the sketch size is in porportion with 1/%d\n", containCompress);
	}
	if(*option_kmer_size){
		isSetKmer = true;
		fprintf(stderr, "-----set kmerSize: %d\n", kmerSize);
	}
	if(*option_sketch_size){
		isJaccard = true;
		fprintf(stderr, "-----set sketchSize:  %d\n", sketchSize);
	}
	if(*option_threshold){
		fprintf(stderr, "-----set threshold:  %d\n", threshold);
	}


#ifndef GREEDY_CLUST
//======clust-mst=========================================================================
	if(*option_premsted && !*option_append){
		clust_from_mst(folder_path, outputFile, is_newick_tree, threshold, threads);
		return 0;
	}
	if(*option_append && !*option_presketched && !*option_premsted){
		cerr << "ERROR option --append, option --presketched or --premsted needed" << endl;
		return 1;
	}
	if(*option_append && (*option_premsted || *option_presketched)){
		append_clust_mst(folder_path, inputFile, outputFile, is_newick_tree, sketchByFile, minLen, noSave, threshold, threads);
		return 0;
	}
//======clust-mst=========================================================================
#else
//======clust-greedy======================================================================
	if(*option_append && !*option_presketched){
		cerr << "ERROR option --append, option --presketched needed" << endl;
		return 1;
	}
	if(*option_append && *option_presketched){
		append_clust_greedy(folder_path, inputFile, outputFile, sketchByFile, minLen, noSave, threshold, threads);
		return 0;
	}
//======clust-greedy======================================================================
#endif
	
	if(*option_presketched && !*option_append){
		clust_from_sketches(folder_path, outputFile, is_newick_tree, threshold, threads);
		return 0;
	}

	if(!tune_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, isJaccard, kmerSize, threshold, containCompress, sketchSize)){
		return 1;
	}
	
	clust_from_genomes(inputFile, outputFile, is_newick_tree, sketchByFile, kmerSize, sketchSize, threshold,sketchFunc, isContainment, containCompress, minLen, folder_path, noSave, threads);

	return 0;
}//end main













