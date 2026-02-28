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

#ifdef LEIDEN_CLUST
#include "leiden.h"
#endif

#ifdef DBSCAN_CLUST
#include "dbscan.h"
#endif

#include <fstream>
#include <sstream>
#include <sys/sysinfo.h>
#include <sys/stat.h>

#include "CLI11.hpp"
#include "sub_command.h"
#ifdef USE_MPI
#include <mpi.h>
#include <sched.h>
#endif

#ifdef GREEDY_CLUST
#else
#endif


using namespace std;

int main(int argc, char * argv[]){
	#ifdef GREEDY_CLUST
		CLI::App app{"clust-greedy v.2.2.1, greedy incremental clustering module for RabbitTClust"};
	#elif defined(LEIDEN_CLUST)
		CLI::App app{"clust-leiden v.2.2.1, Graph-based community detection (Louvain) clustering module for RabbitTClust"};
	#elif defined(DBSCAN_CLUST)
		CLI::App app{"clust-dbscan v.2.2.1, DBSCAN density-based clustering module for RabbitTClust"};
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
	int kmerSize = 19;
	int sketchSize = 1000;
	int containCompress = 1000;
	int drlevel = 3;
	bool mstLoadSketch = false;
	string mstSketchFile = "sketch.info";
	bool isSetKmer = false;
	uint64_t minLen = 10000;
	string folder_path;
	bool is_newick_tree = false;
	bool is_linkage_matrix = false;
	bool is_auto_threshold = false;
	bool is_stability = false;
	bool is_fast = false;
	bool no_dense = false;
	double dedup_dist = -1.0;
	int reps_per_cluster = 0;
	string buildDB_folder;

	bool noSave = false;
	bool use_inverted_index = true;  // Default: use inverted index optimization

	auto option_threads = app.add_option("-t, --threads", threads,  "set the thread number, default all CPUs of the platform");
	auto option_min_len = app.add_option("-m, --min-length", minLen, "set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000");

	auto option_containment = app.add_option("-c, --containment", containCompress, "use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress");
	auto option_kmer_size = app.add_option("-k, --kmer-size", kmerSize, "set the kmer size");
	auto option_sketch_size = app.add_option("-s, --sketch-size", sketchSize, "set the sketch size for Jaccard Index and Mash distance, default 1000");

	auto flag_input_list = app.add_flag("-l, --list", sketchByFile, "input is genome list, one genome per line");
	auto flag_no_save = app.add_flag("-e, --no-save", noSave, "not save the intermediate files, such as sketches or MST");
	bool save_rep_index = false;
	auto flag_save_rep = app.add_flag("--save-rep", save_rep_index, "save representative inverted index for incremental clustering (greedy only)");
	auto option_threshold = app.add_option("-d, --threshold", threshold, "set the distance threshold for clustering");
	auto option_output = app.add_option("-o, --output", outputFile, "set the output name of cluster result");
	auto option_input = app.add_option("-i, --input", inputFile, "set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)");
	auto option_presketched = app.add_option("--presketched", folder_path, "clustering by the pre-generated sketch files rather than genomes");
	auto flag_is_fast = app.add_flag("--fast", is_fast, "use the kssd algorithm for sketching and distance computing");
#ifdef USE_MPI
	bool is_mpi = false;
	auto flag_mpi = app.add_flag("--mpi", is_mpi, "use MPI for distributed sketching and MST (clust-mst, works with both MinHash and KSSD/--fast)");
#endif
	auto flag_inverted_index = app.add_flag("--inverted-index", use_inverted_index, "use inverted index optimization for greedy clustering (MinHash only)");
	
#ifdef DBSCAN_CLUST
	// DBSCAN parameters (only for DBSCAN mode)
	double dbscan_eps = 0.05;
	int dbscan_minpts = 5;
	int dbscan_knn = 0;  // k-NN parameter: 0 = disabled, >0 = keep k nearest neighbors per point
	int dbscan_max_posting = 0;  // max posting list size to keep (0 = disabled)
	auto option_dbscan_eps = app.add_option("--eps", dbscan_eps, "DBSCAN epsilon parameter (distance threshold, default 0.05)");
	auto option_dbscan_minpts = app.add_option("--minpts", dbscan_minpts, "DBSCAN minPts parameter (minimum points to form cluster, default 5)");
	auto option_dbscan_knn = app.add_option("--knn", dbscan_knn, "k-NN pre-filtering: keep only k nearest neighbors per point (0=disabled, recommended: 500-1000, default 0)");
	auto option_dbscan_max_posting = app.add_option("--max-posting", dbscan_max_posting, "DBSCAN inverted-index pruning: drop hash keys with posting size > max-posting (0=disabled)");
	auto option_drlevel = app.add_option("--drlevel", drlevel, "set the dimention reduction level for Kssd sketches, default 3 with a dimention reduction of 1/4096");
#elif defined(LEIDEN_CLUST)
	double leiden_resolution = 1.0;
	bool use_louvain = false;  // Default: use Leiden (this is clust-leiden!)
	int knn_k = 0;  // k-NN parameter: if > 0, keep only k nearest neighbors per node; if 0, use smart defaults
	string pregraph_path;
	
	// Internal flags (auto-controlled, not exposed as CLI options)
	bool use_edge_parallel = false;
	bool use_warm_start = false;
	bool use_parallel_louvain = false;  // DEPRECATED: old graph-partitioning approach, no longer used
	auto option_leiden_resolution = app.add_option("--resolution", leiden_resolution, "Resolution parameter (higher = more clusters, default 1.0)");
	auto flag_use_louvain = app.add_flag("--louvain", use_louvain, "Use Louvain algorithm (RECOMMENDED: fast & optimized, auto-enables all optimizations: edge-parallel + warm-start + knn=1000)");
	auto option_knn = app.add_option("--knn", knn_k, "k-NN filtering: keep only k nearest neighbors per node (default: 1000 for --louvain, 500 for --leiden; use 0 to disable)");
	auto option_pregraph = app.add_option("--pregraph", pregraph_path, "Cluster from pre-built graph (for fast resolution adjustment)");
	auto option_drlevel = app.add_option("--drlevel", drlevel, "set the dimention reduction level for Kssd sketches, default 3 with a dimention reduction of 1/4096");
#elif !defined(GREEDY_CLUST)
	auto option_premsted = app.add_option("--premsted", folder_path, "clustering by the pre-generated mst files rather than genomes for clust-mst");
	auto flag_newick_tree = app.add_flag("--newick-tree", is_newick_tree, "output the newick tree format file for clust-mst");
	auto flag_linkage_matrix = app.add_flag("--linkage-matrix", is_linkage_matrix, "output the single-linkage linkage matrix for clust-mst");
	auto flag_auto_threshold = app.add_flag("--auto-threshold", is_auto_threshold, "automatically select optimal threshold based on MST edge length distribution");
	auto flag_stability = app.add_flag("--stability", is_stability, "evaluate threshold stability by measuring clustering consistency under small perturbations (works with --auto-threshold or user-specified threshold)");
	auto option_drlevel = app.add_option("--drlevel", drlevel, "set the dimention reduction level for Kssd sketches, default 3 with a dimention reduction of 1/4096");
	auto flag_no_dense = app.add_flag("--no-dense", no_dense, "not calculate the density and ANI datas");
	auto option_dedup_dist = app.add_option("--dedup-dist", dedup_dist, "within each cluster, collapse near-duplicate nodes connected by forest edges with dist <= dedup-dist; output to <output>.dedup");
	auto option_reps_per_cluster = app.add_option("--reps-per-cluster", reps_per_cluster, "select up to k representatives per cluster (after optional dedup); output to <output>.reps");
	auto option_buildDB = app.add_option("--buildDB", buildDB_folder, "build a reusable KSSD sketch+index database into the given folder and exit (clust-mst --fast only)");
#endif
	auto option_append = app.add_option("--append", inputFile, "append genome file or file list with the pre-generated sketch or MST files");

	option_append->excludes(option_input);

	CLI11_PARSE(app, argc, argv);

#ifdef USE_MPI
	int my_rank = 0, comm_sz = 1;
	if (is_mpi) {
		int provided;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
		MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		// mpirun often binds each rank to 1 core. When threads > 1, detect
		// this and expand the affinity so OpenMP threads can use multiple cores.
		// On clusters, prefer: mpirun --bind-to none, or
		//   srun --cpu-bind=none, or --map-by slot:PE=<threads> --bind-to core
		{
			cpu_set_t current;
			CPU_ZERO(&current);
			sched_getaffinity(0, sizeof(current), &current);
			int allowed = CPU_COUNT(&current);
			if (allowed < threads) {
				cpu_set_t all;
				CPU_ZERO(&all);
				long ncpus = sysconf(_SC_NPROCESSORS_ONLN);
				for (long i = 0; i < ncpus; i++)
					CPU_SET(i, &all);
				sched_setaffinity(0, sizeof(all), &all);
				#pragma omp parallel
				{
					sched_setaffinity(0, sizeof(all), &all);
				}
				if (my_rank == 0)
					fprintf(stderr, "-----MPI: rank was bound to %d core(s) but needs %d threads, "
							"reset affinity to all %ld cores\n", allowed, threads, ncpus);
			}
		}
	}
#endif

	// Manual validation: output is required unless we're in --buildDB mode.
	if(buildDB_folder.empty() && option_output->count() == 0){
		cerr << "ERROR: option -o/--output is required (unless --buildDB is used)" << endl;
		return 1;
	}

	if(threads < 1){
		fprintf(stderr, "-----Invalid thread number %d\n", threads);
		return 1;
	}
	if(option_threads){
		fprintf(stderr, "-----set the thread number %d\n", threads);
	}
#ifdef _OPENMP
	omp_set_num_threads(threads);
#endif
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
		fprintf(stderr, "-----set threshold:  %g\n", threshold);
	}


#ifdef GREEDY_CLUST
//======clust-greedy======================================================================
	// Set default threshold for greedy clustering if not provided
	if (!*option_threshold) {
		threshold = 0.05;  // Default threshold for greedy clustering
		cerr << "-----use default threshold: " << threshold << endl;
	}
	
 if (is_fast && *option_presketched && !*option_append) {
    clust_from_sketch_fast(folder_path, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, no_dense, isContainment, threshold, threads, dedup_dist, reps_per_cluster, use_inverted_index, save_rep_index);
    return 0;
} 

  
  if(*option_append && !*option_presketched){
		cerr << "ERROR option --append, option --presketched needed" << endl;
		return 1;
	}
	if(*option_append && *option_presketched){
		if(is_fast){
			append_clust_greedy_fast(folder_path, inputFile, outputFile, sketchByFile, minLen, noSave, threshold, threads, save_rep_index);
		} else {
			append_clust_greedy(folder_path, inputFile, outputFile, sketchByFile, minLen, noSave, threshold, threads, save_rep_index);
		}
		return 0;
	}
//======clust-greedy======================================================================
#elif defined(LEIDEN_CLUST)
//======clust-leiden======================================================================
	// Handle algorithm selection
	if (!*option_threshold) {
		threshold = 0.05;  // Default threshold for graph construction
		cerr << "-----use default threshold: " << threshold << endl;
	}
	
	if (*option_leiden_resolution) {
		cerr << "-----Resolution parameter: " << leiden_resolution << endl;
	}
	
	// Simplified Louvain: --louvain auto-enables all optimizations
	if (use_louvain && !use_edge_parallel) {
		// Auto-enable optimizations for --louvain
		use_edge_parallel = true;
		use_warm_start = true;
		
		// Auto-select knn if not specified
		if (knn_k == 0) {
			knn_k = 1000;  // Default k-NN for large datasets
			cerr << "-----Auto-enabled: edge-parallel + warm-start + knn=" << knn_k << endl;
		}
	}
	
	// Smart knn defaults based on expected dataset size
	if (knn_k == 0) {
		// If knn not specified and not using louvain, use conservative defaults
		knn_k = 500;  // Moderate default for Leiden
		cerr << "-----Auto-selecting k-NN: k=" << knn_k << " (use --knn 0 to disable)" << endl;
	}
	
	if (knn_k > 0 && knn_k < 10) {
		cerr << "WARNING: --knn value too small (" << knn_k << "), recommend at least 50. Using 50." << endl;
		knn_k = 50;
	}
	
	// Report algorithm selection
	if (use_louvain) {
		cerr << "-----Algorithm: Louvain (Optimized)" << endl;
		cerr << "  - Edge-parallel construction" << endl;
		if (use_warm_start) {
			cerr << "  - Warm-start initialization" << endl;
		}
		if (knn_k > 0) {
			cerr << "  - k-NN filtering (k=" << knn_k << ")" << endl;
		}
	} else {
		cerr << "-----Algorithm: Leiden" << endl;
		if (knn_k > 0) {
			cerr << "  - k-NN filtering (k=" << knn_k << ")" << endl;
		}
	}
	
	if (*option_pregraph) {
		cerr << "-----Clustering from pre-built graph (fast resolution adjustment)" << endl;
		clust_from_pregraph_leiden(pregraph_path, outputFile, leiden_resolution, !use_louvain, threads);
		return 0;
	}
	
	if (is_fast && *option_presketched && !*option_append) {
		clust_from_sketch_leiden(folder_path, outputFile, threshold, leiden_resolution, !use_louvain, use_parallel_louvain, use_edge_parallel, use_warm_start, knn_k, threads);
		return 0;
	}
	
	if(*option_presketched && !*option_append){
		cerr << "ERROR: clust-leiden requires --fast option" << endl;
		return 1;
	}
	
	if(!isSetKmer){
		kmerSize = 19;
		cerr << "-----use default kmerSize: " << kmerSize << endl;
	}
	
	if(is_fast){
		if(drlevel < 0 || drlevel > 8){
			cerr << "ERROR: invalid drlevel " << drlevel << ", should be in [0, 8]" << endl;
			return 1;
		}
		clust_from_genome_leiden(inputFile, outputFile, folder_path, sketchByFile, kmerSize, drlevel, minLen, noSave, threshold, leiden_resolution, !use_louvain, use_parallel_louvain, use_edge_parallel, use_warm_start, knn_k, threads);
		return 0;
	} else {
		cerr << "ERROR: clust-leiden requires --fast option" << endl;
		return 1;
	}
//======clust-leiden======================================================================
#elif defined(DBSCAN_CLUST)
//======clust-dbscan======================================================================
	// DBSCAN requires --fast (KSSD)
	if(!is_fast){
		cerr << "ERROR: clust-dbscan requires --fast option" << endl;
		return 1;
	}
	
	cerr << "-----Using DBSCAN clustering" << endl;
	cerr << "-----DBSCAN parameters: eps=" << dbscan_eps << ", minPts=" << dbscan_minpts;
	if(dbscan_knn > 0) {
		cerr << ", knn=" << dbscan_knn;
	}
	if(dbscan_max_posting > 0) {
		cerr << ", max-posting=" << dbscan_max_posting;
	}
	cerr << endl;
	
	if(!isSetKmer){
		kmerSize = 19;
		cerr << "-----use default kmerSize: " << kmerSize << endl;
	}
	
	if(drlevel < 0 || drlevel > 8){
		cerr << "ERROR: invalid drlevel " << drlevel << ", should be in [0, 8]" << endl;
		return 1;
	}
	
	if(*option_presketched && !*option_append){
		clust_from_sketch_dbscan(folder_path, outputFile, sketchByFile, dbscan_eps, dbscan_minpts, threads, dbscan_knn, dbscan_max_posting);
		return 0;
	}
	
	if(*option_append){
		cerr << "ERROR: --append not supported for DBSCAN clustering" << endl;
		return 1;
	}
	
	if(!tune_kssd_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, kmerSize, threshold, drlevel)){
		return 1;
	}
	
	clust_from_genome_dbscan(inputFile, outputFile, folder_path, sketchByFile, kmerSize, drlevel, minLen, noSave, dbscan_eps, dbscan_minpts, threads, dbscan_knn, dbscan_max_posting);
	return 0;
//======clust-dbscan======================================================================
#else
//======clust-mst=========================================================================
	// Set default threshold for MST clustering if not provided
	if (!*option_threshold) {
		threshold = 0.05;  // Default threshold for MST clustering
		cerr << "-----use default threshold: " << threshold << endl;
	}

#ifdef USE_MPI
	if (is_mpi) {
		if (is_fast) {
			if (*option_presketched && !*option_append) {
				clust_from_sketches_fast_MPI(my_rank, comm_sz, 10, drlevel, outputFile, folder_path, is_newick_tree, no_dense, true, isContainment, threshold, noSave, threads);
			} else if (!*option_append) {
				if (!tune_kssd_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, kmerSize, threshold, drlevel)) { MPI_Finalize(); return 1; }
				clust_from_genomes_fast_MPI(my_rank, comm_sz, inputFile, outputFile, folder_path, is_newick_tree, no_dense, sketchByFile, isContainment, kmerSize, threshold, drlevel, minLen, noSave, threads);
			} else {
				cerr << "ERROR: --mpi does not support --append" << endl;
				MPI_Finalize();
				return 1;
			}
		} else {
			if (*option_presketched && !*option_append) {
				clust_from_sketches_MPI(my_rank, comm_sz, outputFile, folder_path, is_newick_tree, no_dense, isContainment, threshold, noSave, threads);
			} else if (!*option_append) {
				if (!tune_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, isJaccard, kmerSize, threshold, containCompress, sketchSize)) { MPI_Finalize(); return 1; }
				clust_from_genomes_MPI(my_rank, comm_sz, inputFile, outputFile, folder_path, is_newick_tree, no_dense, sketchByFile, isContainment, kmerSize, sketchSize, threshold, sketchFunc, containCompress, minLen, noSave, threads);
			} else {
				cerr << "ERROR: --mpi does not support --append" << endl;
				MPI_Finalize();
				return 1;
			}
		}
		MPI_Finalize();
		return 0;
	}
#endif
	
	if(is_fast){
		if(!buildDB_folder.empty()){
			if(!sketchByFile){
				cerr << "ERROR: --buildDB currently requires -l/--list (genome file list or cluster output with file paths)" << endl;
				return 1;
			}
			if(option_input->count() == 0){
				cerr << "ERROR: --buildDB requires -i/--input (a genome list or a *.cluster/*.cluster.dedup file)" << endl;
				return 1;
			}
			build_kssd_db_fast(inputFile, buildDB_folder, isSetKmer, isContainment, (int)minLen, kmerSize, drlevel, threads);
			return 0;
		}
		if(*option_premsted && !*option_append){
			clust_from_mst_fast(folder_path, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, threshold, threads);
			return 0;
		}
		if(*option_presketched && !*option_append){
			clust_from_sketch_fast(folder_path, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, no_dense, isContainment, threshold, threads, dedup_dist, reps_per_cluster, use_inverted_index, save_rep_index);
			return 0;
		}
		if(*option_append && !*option_premsted && !*option_presketched){
			cerr << "ERROR: option --append, option --presketched or --premsted needed" << endl;
			return 1;
		}
		if(*option_append && (*option_presketched || *option_premsted)){
			append_clust_mst_fast(folder_path, inputFile, outputFile, is_newick_tree, is_linkage_matrix, no_dense, sketchByFile, isContainment, minLen, noSave, threshold, threads);
			return 0;
		}
		if(!tune_kssd_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, kmerSize, threshold, drlevel)){
			return 1;
		}
		clust_from_genome_fast(inputFile, outputFile, folder_path, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, sketchByFile, isContainment, kmerSize, threshold, drlevel, minLen, noSave, threads, dedup_dist, reps_per_cluster, save_rep_index);
		return 0;
	}

	if(*option_premsted && !*option_append){
		clust_from_mst(folder_path, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, threshold, threads);
		return 0;
	}
	if(*option_append && !*option_presketched && !*option_premsted){
		cerr << "ERROR option --append, option --presketched or --premsted needed" << endl;
		return 1;
	}
	if(*option_append && (*option_premsted || *option_presketched)){
		append_clust_mst(folder_path, inputFile, outputFile, is_newick_tree, is_linkage_matrix, no_dense, sketchByFile, minLen, noSave, threshold, threads);
		return 0;
	}
//======clust-mst=========================================================================
#endif
	
	if(*option_presketched && !*option_append){
		clust_from_sketches(folder_path, outputFile, is_newick_tree, is_auto_threshold, is_stability, no_dense, threshold, threads, use_inverted_index, save_rep_index);
		return 0;
	}

	if(!tune_parameters(sketchByFile, isSetKmer, inputFile, threads, minLen, isContainment, isJaccard, kmerSize, threshold, containCompress, sketchSize)){
		return 1;
	}
  if(is_fast){
  
    clust_from_genome_fast(inputFile, outputFile, folder_path, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, sketchByFile, isContainment, kmerSize, threshold, drlevel, minLen, noSave, threads, dedup_dist, reps_per_cluster, save_rep_index);
    return 0;

  }
	
	clust_from_genomes(inputFile, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, sketchByFile, no_dense, kmerSize, sketchSize, threshold,sketchFunc, isContainment, containCompress, minLen, folder_path, noSave, threads, use_inverted_index, save_rep_index);

	return 0;
}//end main



