#include "sub_command.h"
#include <assert.h>
#include "cluster_postprocess.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <sys/stat.h>  // For stat()
#include <omp.h>  // For OpenMP functions

#ifdef DBSCAN_CLUST
#include "dbscan.h"
#endif

using namespace std;

#ifdef GREEDY_CLUST
void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads, bool save_rep_index){
	bool isSave = !no_save;
	
	// Check for existing cluster state (MinHash incremental update)
	string state_file = folder_path + "/cluster_state.bin";
	MinHashClusterState minhash_state;
	struct stat buffer;
	bool has_state = (stat(state_file.c_str(), &buffer) == 0);
	
	if (has_state) {
		cerr << "===== Incremental Update Mode (MinHash) =====" << endl;
		cerr << "Found existing cluster state, loading..." << endl;
		has_state = minhash_state.load(state_file);
	}
	
	if (!has_state) {
		// No state file, perform initial clustering (old behavior for backward compatibility)
		int sketch_func_id_0; 
		vector<SketchInfo> pre_sketches; 
		bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0); 
		if(pre_sketch_by_file != sketch_by_file){
			cerr << "ERROR: append_clust_greedy(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
			exit(1);
		}
		int sketch_func_id_1, kmer_size, contain_compress, sketch_size, half_k, half_subk, drlevel;
		bool is_containment;
		read_sketch_parameters(folder_path, sketch_func_id_1, kmer_size, is_containment, contain_compress, sketch_size, half_k, half_subk, drlevel);
		assert(sketch_func_id_0 == sketch_func_id_1);
		cerr << "-----use the same sketch parameters with pre-generated sketches" << endl;
		if(sketch_func_id_0 == 0){
			cerr << "---the kmer size is: " << kmer_size << endl;
			if(is_containment)
				cerr << "---use the AAF distance (variable-sketch-size), the sketch size is in proportion with 1/" << contain_compress << endl;
			else 
				cerr << "---use the Mash distance (fixed-sketch-size), the sketch size is: " << sketch_size << endl;
		}
		else if(sketch_func_id_0 == 1){
			cerr << "---use the KSSD sketches" << endl;
			cerr << "---the half_k is: " << half_k << endl;
			cerr << "---the half_subk is: " << half_subk << endl;
			cerr << "---the drlevel is: " << drlevel << endl;
		}
		cerr << "---the thread number is: " << threads << endl;
		cerr << "---the threshold is: " << threshold << endl;
		string sketch_func;
		if(sketch_func_id_0 == 0) 
			sketch_func = "MinHash";
		else if(sketch_func_id_0 == 1)
			sketch_func = "KSSD";
		vector<SketchInfo> append_sketches;
		string append_folder_path;
		compute_sketches(append_sketches, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, false, threads, nullptr);
		vector<SketchInfo> final_sketches;
		final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
		final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
		if(sketch_by_file)
			sort(final_sketches.begin(), final_sketches.end(), cmpGenomeSize);
		else
			sort(final_sketches.begin(), final_sketches.end(), cmpSeqSize);
		vector<SketchInfo>().swap(pre_sketches);
		vector<SketchInfo>().swap(append_sketches);
		string new_folder_path = currentDataTime();
		if(!no_save){
			string command = "mkdir -p " + new_folder_path;
			system(command.c_str());
			saveSketches(final_sketches, new_folder_path, sketch_by_file, sketch_func, is_containment, contain_compress, sketch_size, kmer_size);
		}
		vector<vector<int>> cluster;
		cluster = greedyCluster(final_sketches, sketch_func_id_0, threshold, threads);
		printResult(cluster, final_sketches, sketch_by_file, output_file);
		cerr << "-----write the cluster result into: " << output_file << endl;
		cerr << "-----the cluster number of " << output_file << " is: " << cluster.size() << endl;
	} else {
		// Incremental update mode (MinHash)
		cerr << "---the threshold is: " << threshold << endl;
		cerr << "---the thread number is: " << threads << endl;
		
		// Reload all existing sketches from files (MinHash objects cannot be easily reconstructed)
		int sketch_func_id_0;
		vector<SketchInfo> pre_sketches;
		bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0);
		if(pre_sketch_by_file != sketch_by_file){
			cerr << "Warning: the input format of append genomes and pre-sketched genome is not same" << endl;
		}
		
		// Rebuild representatives from loaded sketches
		minhash_state.all_sketches = pre_sketches;
		minhash_state.representatives.clear();
		minhash_state.representatives.reserve(minhash_state.representative_ids.size());
		
		// Build mapping from old rep_id to new rep_idx
		vector<int> valid_rep_ids;
		valid_rep_ids.reserve(minhash_state.representative_ids.size());
		
		for (int rep_id : minhash_state.representative_ids) {
			if (rep_id >= 0 && rep_id < pre_sketches.size()) {
				if (pre_sketches[rep_id].minHash != nullptr) {
					minhash_state.representatives.push_back(pre_sketches[rep_id]);
					valid_rep_ids.push_back(rep_id);
				} else {
					cerr << "ERROR: Representative " << rep_id << " has null minHash, cannot continue" << endl;
					exit(1);
				}
			} else {
				cerr << "ERROR: Representative ID " << rep_id << " out of range [0, " << pre_sketches.size() << "), cannot continue" << endl;
				exit(1);
			}
		}
		
		// Update representative_ids to match the new representatives
		minhash_state.representative_ids = valid_rep_ids;
		
		cerr << "Successfully loaded " << minhash_state.representatives.size() << " representatives" << endl;
		
		// Rebuild inverted index from representatives
		cerr << "Rebuilding inverted index..." << endl;
		minhash_state.build_inverted_index();
		cerr << "Inverted index rebuilt: " << minhash_state.inverted_index.size() << " unique hashes" << endl;
		
		// Load sketch parameters to get contain_compress
		int sketch_func_id_1, kmer_size, contain_compress, sketch_size, half_k, half_subk, drlevel;
		bool is_containment;
		read_sketch_parameters(folder_path, sketch_func_id_1, kmer_size, is_containment, contain_compress, sketch_size, half_k, half_subk, drlevel);
		
		// Validate parameters
		if (minhash_state.kmer_size <= 0 || minhash_state.sketch_size <= 0) {
			cerr << "ERROR: Invalid kmer_size or sketch_size: kmer=" << minhash_state.kmer_size 
			     << ", sketch=" << minhash_state.sketch_size << endl;
			exit(1);
		}
		
		// Use contain_compress from loaded parameters (not 0)
		if (minhash_state.is_containment && contain_compress <= 0) {
			cerr << "ERROR: is_containment is true but contain_compress is " << contain_compress << endl;
			exit(1);
		}
		
		string sketch_func = "MinHash";
		vector<SketchInfo> new_sketches;
		string new_folder_path;
		cerr << "Computing sketches for new genomes..." << endl;
		cerr << "  Parameters: kmer_size=" << minhash_state.kmer_size 
		     << ", sketch_size=" << minhash_state.sketch_size 
		     << ", is_containment=" << minhash_state.is_containment 
		     << ", contain_compress=" << contain_compress << endl;
		compute_sketches(new_sketches, input_file, new_folder_path, sketch_by_file, min_len, minhash_state.kmer_size, minhash_state.sketch_size, sketch_func, minhash_state.is_containment, contain_compress, false, threads, nullptr);
		
		cerr << "New genomes sketched: " << new_sketches.size() << endl;
		
		if (new_sketches.empty()) {
			cerr << "ERROR: No new sketches generated" << endl;
			return;
		}
		
		cerr << "Starting incremental clustering..." << endl;
		vector<vector<int>> cluster = MinHashIncrementalCluster(minhash_state, new_sketches, threads);
		cerr << "Incremental clustering completed" << endl;
		
		if (!no_save && save_rep_index) {
			minhash_state.save(folder_path + "/cluster_state.bin");
		}
		
		printResult(cluster, minhash_state.all_sketches, sketch_by_file, output_file);
		cerr << "-----write the cluster result into: " << output_file << endl;
		cerr << "-----the cluster number of " << output_file << " is: " << cluster.size() << endl;
		cerr << "-----updated cluster state saved" << endl;
	}
}

void append_clust_greedy_fast(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads, bool save_rep_index){
	bool isSave = !no_save;
	
	// --fast means KSSD, so we only handle KSSD here
	string state_file = folder_path + "/cluster_state.bin";
	KssdClusterState kssd_state;
	struct stat buffer;
	bool has_state = (stat(state_file.c_str(), &buffer) == 0);
	
	if (has_state) {
		cerr << "===== Incremental Update Mode (KSSD) =====" << endl;
		cerr << "Found existing cluster state, loading..." << endl;
		has_state = kssd_state.load(state_file);
	}
	
	if (!has_state) {
			cerr << "===== Initial State Building Mode (KSSD) =====" << endl;
			cerr << "No existing state found, building state from pre-sketched genomes..." << endl;
			
			vector<KssdSketchInfo> pre_sketches; 
			KssdParameters pre_info;
			bool pre_sketch_by_file = loadKssdSketches(folder_path, threads, pre_sketches, pre_info); 
			if(pre_sketch_by_file != sketch_by_file){
				cerr << "Warning: the input format of append genomes and pre-sketched genome is not same" << endl;
			}
			
			int kmer_size = pre_info.half_k * 2;
			int drlevel = pre_info.drlevel;
			
			cerr << "-----use the same sketch parameters with pre-generated sketches" << endl;
			cerr << "---use the KSSD sketches" << endl;
			cerr << "---the half_k is: " << pre_info.half_k << endl;
			cerr << "---the half_subk is: " << pre_info.half_subk << endl;
			cerr << "---the drlevel is: " << drlevel << endl;
			cerr << "---the threshold is: " << threshold << endl;
			
			vector<KssdSketchInfo> append_sketches;
			KssdParameters append_info;
			string append_folder_path;
			compute_kssd_sketches(append_sketches, append_info, isSave, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, drlevel, threads);
			
			kssd_state = KssdInitialClusterWithState(pre_sketches, pre_info, threshold, threads, kmer_size);
			
			if (!no_save && save_rep_index) {
				kssd_state.save(folder_path + "/cluster_state.bin");
			}
			
			vector<vector<int>> cluster = KssdIncrementalCluster(kssd_state, append_sketches, threads);
			
			if (!no_save && save_rep_index) {
				kssd_state.save(folder_path + "/cluster_state.bin");
			}
			
			printKssdResult(cluster, kssd_state.all_sketches, sketch_by_file, output_file);
			cerr << "-----write the cluster result into: " << output_file << endl;
			cerr << "-----the cluster number of " << output_file << " is: " << cluster.size() << endl;
			cerr << "-----saved cluster state (with inverted index) for future incremental updates" << endl;
			
		} else {
			cerr << "---the threshold is: " << threshold << endl;
			cerr << "---the thread number is: " << threads << endl;
			
			vector<KssdSketchInfo> new_sketches;
			KssdParameters new_info;
			string new_folder_path;
			compute_kssd_sketches(new_sketches, new_info, isSave, input_file, new_folder_path, sketch_by_file, min_len, kssd_state.kmer_size, kssd_state.params.drlevel, threads);
			
			cerr << "New genomes sketched: " << new_sketches.size() << endl;
			
			vector<vector<int>> cluster = KssdIncrementalCluster(kssd_state, new_sketches, threads);
			
			if (!no_save && save_rep_index) {
				kssd_state.save(folder_path + "/cluster_state.bin");
			}
			
			printKssdResult(cluster, kssd_state.all_sketches, sketch_by_file, output_file);
			cerr << "-----write the cluster result into: " << output_file << endl;
			cerr << "-----the cluster number of " << output_file << " is: " << cluster.size() << endl;
			cerr << "-----updated cluster state saved" << endl;
	}
}

// =====================================================================
// RepDB workflow functions
// =====================================================================

void repdb_build_from_sketch(string folder_path, string db_path, string output_file, double threshold, int threads) {
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	bool sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);

	int kmer_size = info.half_k * 2;
	cerr << "===== RepDB Build (from pre-sketched) =====" << endl;
	cerr << "  Genomes:    " << sketches.size() << endl;
	cerr << "  Threshold:  " << threshold << endl;
	cerr << "  Kmer size:  " << kmer_size << endl;

	KssdClusterState state = KssdInitialClusterWithState(sketches, info, threshold, threads, kmer_size);

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printKssdResult(state.clusters, state.all_sketches, sketchByFile, output_file, threshold);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	cerr << "\n===== RepDB Build Summary =====" << endl;
	cerr << "  Total genomes:    " << sketches.size() << endl;
	cerr << "  Representatives:  " << state.representatives.size() << endl;
	cerr << "  Compression:      " << std::fixed << std::setprecision(2)
		 << (1.0 - (double)state.representatives.size() / sketches.size()) * 100.0 << "%" << endl;
	cerr << "  RepDB saved to:   " << db_path << endl;
	cerr << "===============================" << endl;
}

void repdb_build_from_genome(string input_file, string db_path, string output_file, bool sketch_by_file, int min_len, int kmer_size, int drlevel, double threshold, int threads) {
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	string folder_path;
	compute_kssd_sketches(sketches, info, true, input_file, folder_path, sketch_by_file, min_len, kmer_size, drlevel, threads);

	int actual_kmer_size = info.half_k * 2;
	cerr << "===== RepDB Build (from genomes) =====" << endl;
	cerr << "  Genomes:    " << sketches.size() << endl;
	cerr << "  Threshold:  " << threshold << endl;
	cerr << "  Kmer size:  " << actual_kmer_size << endl;

	KssdClusterState state = KssdInitialClusterWithState(sketches, info, threshold, threads, actual_kmer_size);

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printKssdResult(state.clusters, state.all_sketches, sketch_by_file, output_file, threshold);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	cerr << "\n===== RepDB Build Summary =====" << endl;
	cerr << "  Total genomes:    " << sketches.size() << endl;
	cerr << "  Representatives:  " << state.representatives.size() << endl;
	cerr << "  Compression:      " << std::fixed << std::setprecision(2)
		 << (1.0 - (double)state.representatives.size() / sketches.size()) * 100.0 << "%" << endl;
	cerr << "  RepDB saved to:   " << db_path << endl;
	cerr << "===============================" << endl;
}

void repdb_query(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int topk, int threads) {
	KssdClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<KssdSketchInfo> queries;
	KssdParameters qinfo;
	string qfolder;
	compute_kssd_sketches(queries, qinfo, false, input_file, qfolder, sketch_by_file, min_len, state.kmer_size, state.params.drlevel, threads);

	cerr << "===== RepDB Query =====" << endl;
	cerr << "  Query genomes:  " << queries.size() << endl;
	cerr << "  Top-k:          " << topk << endl;
	cerr << "  DB reps:        " << state.representatives.size() << endl;

	FILE* fp = fopen(output_file.c_str(), "w");
	if (!fp) {
		cerr << "ERROR: Cannot open output file: " << output_file << endl;
		exit(1);
	}

	fprintf(fp, "#query\trank\trep_name\tdistance\tcluster_id\tcluster_size\n");

	for (size_t i = 0; i < queries.size(); i++) {
		auto results = state.query_topk(queries[i], topk, threads);
		string qname = queries[i].fileName;
		if (qname.empty()) qname = "query_" + std::to_string(i);

		if (results.empty()) {
			fprintf(fp, "%s\t0\tno_match\t-1\t-1\t0\n", qname.c_str());
		} else {
			for (int r = 0; r < (int)results.size(); r++) {
				fprintf(fp, "%s\t%d\t%s\t%.6f\t%d\t%d\n",
					qname.c_str(),
					r + 1,
					results[r].genome_name.c_str(),
					results[r].distance,
					results[r].cluster_id,
					results[r].cluster_size);
			}
		}

		if ((i + 1) % 10000 == 0) {
			cerr << "---queried: " << (i + 1) << " / " << queries.size() << endl;
		}
	}
	fclose(fp);

	cerr << "===== Query Results =====" << endl;
	cerr << "  Output: " << output_file << endl;
	cerr << "=========================" << endl;
}

void repdb_assign(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int threads) {
	KssdClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<KssdSketchInfo> queries;
	KssdParameters qinfo;
	string qfolder;
	compute_kssd_sketches(queries, qinfo, false, input_file, qfolder, sketch_by_file, min_len, state.kmer_size, state.params.drlevel, threads);

	cerr << "===== RepDB Assignment =====" << endl;
	cerr << "  Query genomes:  " << queries.size() << endl;
	cerr << "  DB reps:        " << state.representatives.size() << endl;
	cerr << "  Threshold:      " << state.threshold << endl;

	FILE* fp = fopen(output_file.c_str(), "w");
	if (!fp) {
		cerr << "ERROR: Cannot open output file: " << output_file << endl;
		exit(1);
	}

	fprintf(fp, "#query\tassigned_cluster\trep_name\tdistance\tcluster_size\tstatus\n");

	int assigned = 0, unassigned = 0;

	for (size_t i = 0; i < queries.size(); i++) {
		auto result = state.assign(queries[i], threads);
		string qname = queries[i].fileName;
		if (qname.empty()) qname = "query_" + std::to_string(i);

		if (result.rep_idx >= 0) {
			fprintf(fp, "%s\t%d\t%s\t%.6f\t%d\tassigned\n",
				qname.c_str(),
				result.cluster_id,
				result.genome_name.c_str(),
				result.distance,
				result.cluster_size);
			assigned++;
		} else {
			fprintf(fp, "%s\t-1\tunassigned\t-1\t0\tnovel\n", qname.c_str());
			unassigned++;
		}

		if ((i + 1) % 10000 == 0) {
			cerr << "---assigned: " << (i + 1) << " / " << queries.size() << endl;
		}
	}
	fclose(fp);

	cerr << "===== Assignment Results =====" << endl;
	cerr << "  Assigned:    " << assigned << " (" << std::fixed << std::setprecision(1)
		 << (100.0 * assigned / queries.size()) << "%)" << endl;
	cerr << "  Novel:       " << unassigned << " (" << std::fixed << std::setprecision(1)
		 << (100.0 * unassigned / queries.size()) << "%)" << endl;
	cerr << "  Output:      " << output_file << endl;
	cerr << "==============================" << endl;
}

void repdb_append(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int threads) {
	KssdClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<KssdSketchInfo> new_sketches;
	KssdParameters new_info;
	string new_folder;
	compute_kssd_sketches(new_sketches, new_info, false, input_file, new_folder, sketch_by_file, min_len, state.kmer_size, state.params.drlevel, threads);

	int old_rep_count = state.representatives.size();
	int old_total = state.all_sketches.size();

	cerr << "===== RepDB Append =====" << endl;
	cerr << "  Existing reps:    " << old_rep_count << endl;
	cerr << "  Existing genomes: " << old_total << endl;
	cerr << "  New genomes:      " << new_sketches.size() << endl;

	vector<vector<int>> clusters = KssdIncrementalCluster(state, new_sketches, threads);

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printKssdResult(clusters, state.all_sketches, sketch_by_file, output_file, state.threshold);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	int new_rep_count = state.representatives.size();
	cerr << "\n===== Append Summary =====" << endl;
	cerr << "  New reps added:   " << (new_rep_count - old_rep_count) << endl;
	cerr << "  Total reps now:   " << new_rep_count << endl;
	cerr << "  Total genomes:    " << state.all_sketches.size() << endl;
	cerr << "  RepDB updated:    " << db_path << endl;
	cerr << "==========================" << endl;
}

void repdb_stats(string db_path) {
	KssdClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load RepDB from: " << db_path << endl;
		exit(1);
	}
	state.print_stats(std::cout);
}

// =====================================================================
// MinHash RepDB workflow functions
// =====================================================================

void mh_repdb_build_from_sketch(string folder_path, string db_path, string output_file, double threshold, int threads, bool use_inverted_index) {
	vector<SketchInfo> sketches;
	int sketch_func_id;
	bool sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);

	int kmer_size = sketches[0].minHash->getKmerSize();
	int sketch_size = sketches[0].minHash->getSketchSize();
	bool is_containment = sketches[0].isContainment;

	cerr << "===== MinHash RepDB Build (from pre-sketched) =====" << endl;
	cerr << "  Genomes:       " << sketches.size() << endl;
	cerr << "  Threshold:     " << threshold << endl;
	cerr << "  Kmer size:     " << kmer_size << endl;
	cerr << "  Sketch size:   " << sketch_size << endl;
	cerr << "  Containment:   " << (is_containment ? "yes" : "no") << endl;

	MinHashClusterState state = MinHashInitialClusterWithState(sketches, threshold, threads, kmer_size, sketch_size, is_containment);

	state.rep_hash_arrays.resize(state.representatives.size());
	for (size_t i = 0; i < state.representatives.size(); i++) {
		if (state.representatives[i].minHash)
			state.rep_hash_arrays[i] = state.representatives[i].minHash->storeMinHashes();
	}

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printResult(state.clusters, state.all_sketches, sketchByFile, output_file);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	cerr << "\n===== MinHash RepDB Build Summary =====" << endl;
	cerr << "  Total genomes:    " << sketches.size() << endl;
	cerr << "  Representatives:  " << state.representatives.size() << endl;
	cerr << "  Compression:      " << std::fixed << std::setprecision(2)
		 << (1.0 - (double)state.representatives.size() / sketches.size()) * 100.0 << "%" << endl;
	cerr << "  RepDB saved to:   " << db_path << endl;
	cerr << "===============================" << endl;
}

void mh_repdb_build_from_genome(string input_file, string db_path, string output_file, bool sketch_by_file, int min_len, int kmer_size, int sketch_size, string sketch_func, bool is_containment, int contain_compress, double threshold, int threads) {
	vector<SketchInfo> sketches;
	string folder_path;
	compute_sketches(sketches, input_file, folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, true, threads, nullptr);

	int actual_kmer_size = sketches[0].minHash->getKmerSize();
	int actual_sketch_size = sketches[0].minHash->getSketchSize();
	bool actual_containment = sketches[0].isContainment;

	cerr << "===== MinHash RepDB Build (from genomes) =====" << endl;
	cerr << "  Genomes:       " << sketches.size() << endl;
	cerr << "  Threshold:     " << threshold << endl;
	cerr << "  Kmer size:     " << actual_kmer_size << endl;
	cerr << "  Sketch size:   " << actual_sketch_size << endl;

	MinHashClusterState state = MinHashInitialClusterWithState(sketches, threshold, threads, actual_kmer_size, actual_sketch_size, actual_containment);

	state.rep_hash_arrays.resize(state.representatives.size());
	for (size_t i = 0; i < state.representatives.size(); i++) {
		if (state.representatives[i].minHash)
			state.rep_hash_arrays[i] = state.representatives[i].minHash->storeMinHashes();
	}

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printResult(state.clusters, state.all_sketches, sketch_by_file, output_file);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	cerr << "\n===== MinHash RepDB Build Summary =====" << endl;
	cerr << "  Total genomes:    " << sketches.size() << endl;
	cerr << "  Representatives:  " << state.representatives.size() << endl;
	cerr << "  Compression:      " << std::fixed << std::setprecision(2)
		 << (1.0 - (double)state.representatives.size() / sketches.size()) * 100.0 << "%" << endl;
	cerr << "  RepDB saved to:   " << db_path << endl;
	cerr << "===============================" << endl;
}

void mh_repdb_query(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int topk, int threads) {
	MinHashClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load MinHash RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<SketchInfo> queries;
	string qfolder;
	compute_sketches(queries, input_file, qfolder, sketch_by_file, min_len, state.kmer_size, state.sketch_size, "MinHash", state.is_containment, 1000, false, threads, nullptr);

	cerr << "===== MinHash RepDB Query =====" << endl;
	cerr << "  Query genomes:  " << queries.size() << endl;
	cerr << "  Top-k:          " << topk << endl;
	cerr << "  DB reps:        " << state.representatives.size() << endl;

	FILE* fp = fopen(output_file.c_str(), "w");
	if (!fp) { cerr << "ERROR: Cannot open output file: " << output_file << endl; exit(1); }

	fprintf(fp, "#query\trank\trep_name\tdistance\tcluster_id\tcluster_size\n");

	for (size_t i = 0; i < queries.size(); i++) {
		auto results = state.query_topk(queries[i], topk, threads);
		string qname = queries[i].fileName;
		if (qname.empty()) qname = "query_" + std::to_string(i);

		if (results.empty()) {
			fprintf(fp, "%s\t0\tno_match\t-1\t-1\t0\n", qname.c_str());
		} else {
			for (int r = 0; r < (int)results.size(); r++) {
				fprintf(fp, "%s\t%d\t%s\t%.6f\t%d\t%d\n",
					qname.c_str(), r + 1,
					results[r].genome_name.c_str(),
					results[r].distance,
					results[r].cluster_id,
					results[r].cluster_size);
			}
		}
		if ((i + 1) % 10000 == 0)
			cerr << "---queried: " << (i + 1) << " / " << queries.size() << endl;
	}
	fclose(fp);
	cerr << "===== Query Results =====" << endl;
	cerr << "  Output: " << output_file << endl;
	cerr << "=========================" << endl;
}

void mh_repdb_assign(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int threads) {
	MinHashClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load MinHash RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<SketchInfo> queries;
	string qfolder;
	compute_sketches(queries, input_file, qfolder, sketch_by_file, min_len, state.kmer_size, state.sketch_size, "MinHash", state.is_containment, 1000, false, threads, nullptr);

	cerr << "===== MinHash RepDB Assignment =====" << endl;
	cerr << "  Query genomes:  " << queries.size() << endl;
	cerr << "  DB reps:        " << state.representatives.size() << endl;
	cerr << "  Threshold:      " << state.threshold << endl;

	FILE* fp = fopen(output_file.c_str(), "w");
	if (!fp) { cerr << "ERROR: Cannot open output file: " << output_file << endl; exit(1); }

	fprintf(fp, "#query\tassigned_cluster\trep_name\tdistance\tcluster_size\tstatus\n");

	int assigned = 0, unassigned = 0;
	for (size_t i = 0; i < queries.size(); i++) {
		auto result = state.assign(queries[i], threads);
		string qname = queries[i].fileName;
		if (qname.empty()) qname = "query_" + std::to_string(i);

		if (result.rep_idx >= 0) {
			fprintf(fp, "%s\t%d\t%s\t%.6f\t%d\tassigned\n",
				qname.c_str(), result.cluster_id,
				result.genome_name.c_str(), result.distance, result.cluster_size);
			assigned++;
		} else {
			fprintf(fp, "%s\t-1\tunassigned\t-1\t0\tnovel\n", qname.c_str());
			unassigned++;
		}
		if ((i + 1) % 10000 == 0)
			cerr << "---assigned: " << (i + 1) << " / " << queries.size() << endl;
	}
	fclose(fp);

	cerr << "===== Assignment Results =====" << endl;
	cerr << "  Assigned:    " << assigned << " (" << std::fixed << std::setprecision(1)
		 << (100.0 * assigned / queries.size()) << "%)" << endl;
	cerr << "  Novel:       " << unassigned << " (" << std::fixed << std::setprecision(1)
		 << (100.0 * unassigned / queries.size()) << "%)" << endl;
	cerr << "  Output:      " << output_file << endl;
	cerr << "==============================" << endl;
}

static void printRepDBClusterResult(vector<vector<int>>& clusters, vector<SketchInfo>& sketches, string outputFile, double threshold) {
	FILE* fp = fopen(outputFile.c_str(), "w");
	if (!fp) { cerr << "ERROR: cannot open file: " << outputFile << endl; exit(1); }

	if (threshold >= 0.0) {
		fprintf(fp, "# Clustering threshold: %.6f\n", threshold);
		fprintf(fp, "# Total clusters: %zu\n", clusters.size());
		fprintf(fp, "#\n");
	}

	for (size_t i = 0; i < clusters.size(); i++) {
		fprintf(fp, "the cluster %zu is: \n", i);
		for (size_t j = 0; j < clusters[i].size(); j++) {
			int curId = clusters[i][j];
			if (curId < 0 || curId >= (int)sketches.size()) continue;
			const char* fname = sketches[curId].fileName.empty() ? "N/A" : sketches[curId].fileName.c_str();
			const char* seqName = "N/A";
			const char* seqComment = "";
			if (!sketches[curId].fileSeqs.empty()) {
				seqName = sketches[curId].fileSeqs[0].name.c_str();
				seqComment = sketches[curId].fileSeqs[0].comment.c_str();
			}
			fprintf(fp, "\t%5zu\t%6d\t%12ldnt\t%20s\t%20s\t%s\n",
				j, curId, (long)sketches[curId].totalSeqLength, fname, seqName, seqComment);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void mh_repdb_append(string db_path, string input_file, string output_file, bool sketch_by_file, int min_len, int threads) {
	MinHashClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load MinHash RepDB from: " << db_path << endl;
		exit(1);
	}

	vector<SketchInfo> new_sketches;
	string nfolder;
	compute_sketches(new_sketches, input_file, nfolder, sketch_by_file, min_len, state.kmer_size, state.sketch_size, "MinHash", state.is_containment, 1000, false, threads, nullptr);

	int old_rep_count = state.representatives.size();
	int old_total = state.all_sketches.size();

	cerr << "===== MinHash RepDB Append =====" << endl;
	cerr << "  Existing reps:    " << old_rep_count << endl;
	cerr << "  Existing genomes: " << old_total << endl;
	cerr << "  New genomes:      " << new_sketches.size() << endl;

	MinHashIncrementalCluster(state, new_sketches, threads);

	state.save_repdb(db_path);

	if (!output_file.empty()) {
		printRepDBClusterResult(state.clusters, state.all_sketches, output_file, state.threshold);
		cerr << "-----write the cluster result into: " << output_file << endl;
	}

	int new_rep_count = state.representatives.size();
	cerr << "\n===== Append Summary =====" << endl;
	cerr << "  New reps added:   " << (new_rep_count - old_rep_count) << endl;
	cerr << "  Total reps now:   " << new_rep_count << endl;
	cerr << "  Total genomes:    " << state.all_sketches.size() << endl;
	cerr << "  RepDB updated:    " << db_path << endl;
	cerr << "==========================" << endl;
}

void mh_repdb_stats(string db_path) {
	MinHashClusterState state;
	if (!state.load_repdb(db_path)) {
		cerr << "ERROR: Failed to load MinHash RepDB from: " << db_path << endl;
		exit(1);
	}
	state.print_stats(std::cout);
}

#endif

#ifndef GREEDY_CLUST
void append_clust_mst_fast(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, bool isContainment, int min_len, bool no_save, double threshold, int threads){
	bool isSave = !no_save;
	int sketch_func_id_0; 
	vector<KssdSketchInfo> pre_sketches; 
	KssdParameters pre_info;
	bool pre_sketch_by_file = loadKssdSketches(folder_path, threads, pre_sketches, pre_info); 
	if(pre_sketch_by_file != sketch_by_file){
		cerr << "Warning: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
		cerr << "the output cluster file may not have the genome file name" << endl;
	}
	vector<EdgeInfo> pre_mst;
	loadMST(folder_path, pre_mst);
	int kmer_size = pre_info.half_k * 2;
	int drlevel = pre_info.drlevel;

	cerr << "---the thread number is: " << threads << endl;
	cerr << "---the threshold is: " << threshold << endl;

	vector<KssdSketchInfo> append_sketches;
	KssdParameters append_info;
	string append_folder_path;
	compute_kssd_sketches(append_sketches, append_info, isSave, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, drlevel, threads);

	vector<KssdSketchInfo> final_sketches;
	int pre_sketch_size = pre_sketches.size();
	final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
	final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
	vector<KssdSketchInfo>().swap(pre_sketches);
	vector<KssdSketchInfo>().swap(append_sketches);
	string new_folder_path = append_folder_path;
	if(!no_save){
		string command = "mkdir -p " + new_folder_path;
		system(command.c_str());
		saveKssdSketches(final_sketches, pre_info, new_folder_path, sketch_by_file);
	}

	bool use64 = append_info.half_k - append_info.drlevel > 8;
	KssdInvertedIndex inverted_index;
	inverted_index.init(use64);

	cerr << "-----loading pre-existing KSSD inverted index from: " << folder_path << endl;
	double t_idx0 = get_sec();
	string pre_index_file = folder_path + '/' + "kssd.sketch.index";
	string pre_dict_file = folder_path + '/' + "kssd.sketch.dict";
	if(use64){
		FILE* fp_index = fopen(pre_index_file.c_str(), "rb");
		if(!fp_index){ cerr << "ERROR: append_clust_mst_fast(), cannot open index file: " << pre_index_file << endl; exit(1); }
		size_t hash_number;
		fread(&hash_number, sizeof(size_t), 1, fp_index);
		uint64_t* hash_arr = new uint64_t[hash_number];
		uint32_t* hash_size_arr = new uint32_t[hash_number];
		fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);

		FILE* fp_dict = fopen(pre_dict_file.c_str(), "rb");
		if(!fp_dict){ cerr << "ERROR: append_clust_mst_fast(), cannot open dict file: " << pre_dict_file << endl; exit(1); }
		uint64_t total_size = 0;
		for(size_t i = 0; i < hash_number; i++) total_size += hash_size_arr[i];
		uint32_t* dict_buf = new uint32_t[total_size];
		fread(dict_buf, sizeof(uint32_t), total_size, fp_dict);
		fclose(fp_dict);

		inverted_index.hash_map_64.reserve(hash_number);
		size_t off = 0;
		for(size_t i = 0; i < hash_number; i++){
			size_t len = hash_size_arr[i];
			inverted_index.hash_map_64[hash_arr[i]] = vector<uint32_t>(dict_buf + off, dict_buf + off + len);
			off += len;
		}
		delete[] dict_buf;
		delete[] hash_arr;
		delete[] hash_size_arr;
	} else {
		FILE* fp_index = fopen(pre_index_file.c_str(), "rb");
		if(!fp_index){ cerr << "ERROR: append_clust_mst_fast(), cannot open index file: " << pre_index_file << endl; exit(1); }
		size_t hash_number;
		fread(&hash_number, sizeof(size_t), 1, fp_index);
		uint32_t* hash_arr = new uint32_t[hash_number];
		uint32_t* hash_size_arr = new uint32_t[hash_number];
		fread(hash_arr, sizeof(uint32_t), hash_number, fp_index);
		fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);

		FILE* fp_dict = fopen(pre_dict_file.c_str(), "rb");
		if(!fp_dict){ cerr << "ERROR: append_clust_mst_fast(), cannot open dict file: " << pre_dict_file << endl; exit(1); }
		uint64_t total_size = 0;
		for(size_t i = 0; i < hash_number; i++) total_size += hash_size_arr[i];
		uint32_t* dict_buf = new uint32_t[total_size];
		fread(dict_buf, sizeof(uint32_t), total_size, fp_dict);
		fclose(fp_dict);

		inverted_index.hash_map_32.reserve(hash_number);
		size_t off = 0;
		for(size_t i = 0; i < hash_number; i++){
			size_t len = hash_size_arr[i];
			inverted_index.hash_map_32[hash_arr[i]] = vector<uint32_t>(dict_buf + off, dict_buf + off + len);
			off += len;
		}
		delete[] dict_buf;
		delete[] hash_arr;
		delete[] hash_size_arr;
	}
	double t_idx1 = get_sec();
	cerr << "-----loaded pre-existing inverted index in " << t_idx1 - t_idx0 << " seconds" << endl;

	cerr << "-----appending " << final_sketches.size() - pre_sketch_size << " new genomes to inverted index..." << endl;
	if(use64){
		for(size_t i = pre_sketch_size; i < final_sketches.size(); i++){
			for(size_t j = 0; j < final_sketches[i].hash64_arr.size(); j++){
				inverted_index.hash_map_64[final_sketches[i].hash64_arr[j]].push_back((uint32_t)i);
			}
		}
	} else {
		for(size_t i = pre_sketch_size; i < final_sketches.size(); i++){
			for(size_t j = 0; j < final_sketches[i].hash32_arr.size(); j++){
				inverted_index.hash_map_32[final_sketches[i].hash32_arr[j]].push_back((uint32_t)i);
			}
		}
	}
	double t_idx2 = get_sec();
	cerr << "-----appended new genomes to inverted index in " << t_idx2 - t_idx1 << " seconds" << endl;

	int ** pre_dense_arr;
	uint64_t* pre_ani_arr;
	int pre_dense_span;
	int pre_genome_number;
	if(!no_dense){
		loadDense(pre_dense_arr, folder_path, pre_dense_span, pre_genome_number);
	}

	int ** dense_arr;
	int dense_span = DENSE_SPAN;
	uint64_t* ani_arr;
	vector<EdgeInfo> append_mst = compute_kssd_mst(final_sketches, append_info, new_folder_path, pre_sketch_size, no_dense, isContainment, threads, dense_arr, dense_span, ani_arr, threshold, &inverted_index);
	vector<EdgeInfo> final_graph;
	final_graph.insert(final_graph.end(), pre_mst.begin(), pre_mst.end());
	final_graph.insert(final_graph.end(), append_mst.begin(), append_mst.end());
	vector<EdgeInfo>().swap(pre_mst);
	vector<EdgeInfo>().swap(append_mst);
	sort(final_graph.begin(), final_graph.end(), cmpEdge);
	vector<EdgeInfo> final_mst = kruskalAlgorithm(final_graph, final_sketches.size());
	vector<EdgeInfo>().swap(final_graph);
	if(is_newick_tree){
		string output_newick_file = output_file + ".newick.tree";
		print_kssd_newick_tree(final_sketches, final_mst, pre_sketch_by_file, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	if(is_phylip_tree){
		string output_phylip_file = output_file + ".phylip.tree";
		print_kssd_phylip_tree(final_sketches, final_mst, pre_sketch_by_file, output_phylip_file);
		cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
	}
	if(is_nexus_tree){
		string output_nexus_file = output_file + ".nexus.tree";
		print_kssd_nexus_tree(final_sketches, final_mst, pre_sketch_by_file, output_nexus_file);
		cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
	}
	if(is_linkage_matrix){
		string output_linkage_file = output_file + ".linkage.txt";
		print_kssd_linkage_matrix(final_sketches, final_mst, output_linkage_file);
		cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
	}

	vector<EdgeInfo> forest = generateForest(final_mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, final_sketches.size());
	printKssdResult(tmpClust, final_sketches, pre_sketch_by_file, output_file, threshold);
	cerr << "-----write the cluster result into: " << output_file << endl;
	cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
	if(!no_save){
		saveKssdMST(final_sketches, final_mst, new_folder_path, sketch_by_file);
		transSketchesFromIndex(inverted_index, append_info, new_folder_path);
	}

	if(!no_dense){
		loadANI(folder_path, pre_ani_arr, sketch_func_id_0);
		for(int i = 0; i < 101; i++)
			ani_arr[i] += pre_ani_arr[i];
		for(int i = 0; i < pre_dense_span; i++){
			for(int j = 0; j < pre_genome_number; j++){
				dense_arr[i][j] += pre_dense_arr[i][j];
			}
		}

		if(!no_save){
			saveANI(new_folder_path, ani_arr, sketch_func_id_0);
			saveDense(new_folder_path, dense_arr, dense_span, final_sketches.size());
		}

		int alpha = 2;
		int denseIndex = threshold / 0.01;
		vector<int> totalNoiseArr;
		for(int i = 0; i < tmpClust.size(); i++){
			if(tmpClust[i].size() == 1) continue;
			vector<PairInt> curDenseArr;
			set<int> denseSet;
			for(int j = 0; j < tmpClust[i].size(); j++){
				int element = tmpClust[i][j];
				PairInt p(element, dense_arr[denseIndex][element]);
				denseSet.insert(dense_arr[denseIndex][element]);
				curDenseArr.push_back(p);
			}
			vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
			totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
		}
		cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
		forest = modifyForest(forest, totalNoiseArr, threads);
		vector<vector<int>> cluster = generateClusterWithBfs(forest, final_sketches.size());
		string outputFileNew = output_file + ".removeNoise";
		printKssdResult(cluster, final_sketches, pre_sketch_by_file, outputFileNew);
		cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
		cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
	}

}

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
	int sketch_func_id_0; 
	vector<SketchInfo> pre_sketches; 
	bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0); 
	if(pre_sketch_by_file != sketch_by_file){
		cerr << "Warning: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
		cerr << "the output cluster file may not have the genome file name" << endl;
	}
	vector<EdgeInfo> pre_mst;
	loadMST(folder_path, pre_mst);
	int sketch_func_id_1, kmer_size, contain_compress, sketch_size, half_k, half_subk, drlevel;
	bool is_containment;
	read_sketch_parameters(folder_path, sketch_func_id_1, kmer_size, is_containment, contain_compress, sketch_size, half_k, half_subk, drlevel);
	assert(sketch_func_id_0 == sketch_func_id_1);
	cerr << "-----use the same sketch parameters with pre-generated sketches" << endl;
	if(sketch_func_id_0 == 0){
		cerr << "---the kmer size is: " << kmer_size << endl;
		if(is_containment)
			cerr << "---use the AAF distance (variable-sketch-size), the sketch size is in proportion with 1/" << contain_compress << endl;
		else 
			cerr << "---use the Mash distance (fixed-sketch-size), the sketch size is: " << sketch_size << endl;
	}
	else if(sketch_func_id_0 == 1){
		cerr << "---use the KSSD sketches" << endl;
		cerr << "---the half_k is: " << half_k << endl;
		cerr << "---the half_subk is: " << half_subk << endl;
		cerr << "---the drlevel is: " << drlevel << endl;
	}
	cerr << "---the thread number is: " << threads << endl;
	cerr << "---the threshold is: " << threshold << endl;
	string sketch_func;
	if(sketch_func_id_0 == 0) 
		sketch_func = "MinHash";
	else if(sketch_func_id_0 == 1)
		sketch_func = "KSSD";
	vector<SketchInfo> append_sketches;
	string append_folder_path;
	compute_sketches(append_sketches, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, false, threads);

	vector<SketchInfo> final_sketches;
	int pre_sketch_size = pre_sketches.size();
	final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
	final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
	vector<SketchInfo>().swap(pre_sketches);
	vector<SketchInfo>().swap(append_sketches);
	string new_folder_path = currentDataTime();
	if(!no_save){
		string command = "mkdir -p " + new_folder_path;
		system(command.c_str());
		saveSketches(final_sketches, new_folder_path, sketch_by_file, sketch_func, is_containment, contain_compress, sketch_size, kmer_size);
	}

	int ** pre_dense_arr;
	uint64_t* pre_ani_arr;
	int pre_dense_span;
	int pre_genome_number;
	if(!no_dense){
		loadDense(pre_dense_arr, folder_path, pre_dense_span, pre_genome_number);
	}

	int ** dense_arr;
	int dense_span = DENSE_SPAN;
	uint64_t* ani_arr;
	
	vector<EdgeInfo> append_mst;
	if(sketch_func_id_0 == 0 && final_sketches.size() > 0 && final_sketches[0].minHash != nullptr) {
		MinHashInvertedIndex inverted_index;
		double t_idx0 = get_sec();

		bool loaded = loadMinHashIndex(folder_path, inverted_index);
		if(loaded){
			double t_idx1 = get_sec();
			cerr << "-----loaded pre-existing MinHash inverted index in " << t_idx1 - t_idx0 << " seconds" << endl;
			cerr << "-----appending " << final_sketches.size() - pre_sketch_size << " new genomes to MinHash inverted index..." << endl;
			for(size_t i = pre_sketch_size; i < final_sketches.size(); i++){
				if(final_sketches[i].minHash != nullptr) {
					vector<uint64_t> hashes = final_sketches[i].minHash->storeMinHashes();
					for(uint64_t hash : hashes) {
						inverted_index.hash_map[hash].push_back((uint32_t)i);
					}
				}
			}
			double t_idx2 = get_sec();
			cerr << "-----appended new genomes to MinHash inverted index in " << t_idx2 - t_idx1 << " seconds" << endl;
		} else {
			cerr << "-----no pre-existing MinHash inverted index found, building from scratch..." << endl;
			vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);
			#pragma omp parallel num_threads(threads)
			{
				int tid = omp_get_thread_num();
				auto& local_index = thread_local_indices[tid];
				#pragma omp for schedule(dynamic, 100)
				for(size_t i = 0; i < final_sketches.size(); i++){
					if(final_sketches[i].minHash != nullptr) {
						final_sketches[i].id = i;
						vector<uint64_t> hashes = final_sketches[i].minHash->storeMinHashes();
						for(uint64_t hash : hashes) {
							local_index[hash].push_back(i);
						}
					}
				}
			}
			for(int tid = 0; tid < threads; tid++){
				for(auto& entry : thread_local_indices[tid]){
					inverted_index.hash_map[entry.first].insert(
						inverted_index.hash_map[entry.first].end(),
						entry.second.begin(), entry.second.end()
					);
				}
			}
			cerr << "-----MinHash inverted index built: " << inverted_index.hash_map.size() << " unique hashes" << endl;
		}

		append_mst = compute_minhash_mst(final_sketches, pre_sketch_size, no_dense, is_containment, threads, dense_arr, dense_span, ani_arr, threshold, kmer_size, &inverted_index);

		if(!no_save){
			inverted_index.save_to_file(new_folder_path);
		}
	} else {
		append_mst = modifyMST(final_sketches, pre_sketch_size, sketch_func_id_0, threads, no_dense, dense_arr, dense_span, ani_arr);
	}
	vector<EdgeInfo> final_graph;
	final_graph.insert(final_graph.end(), pre_mst.begin(), pre_mst.end());
	final_graph.insert(final_graph.end(), append_mst.begin(), append_mst.end());
	vector<EdgeInfo>().swap(pre_mst);
	vector<EdgeInfo>().swap(append_mst);
	sort(final_graph.begin(), final_graph.end(), cmpEdge);
	vector<EdgeInfo> final_mst = kruskalAlgorithm(final_graph, final_sketches.size());
	vector<EdgeInfo>().swap(final_graph);
	if(is_newick_tree){
		string output_newick_file = output_file + ".newick.tree";
		print_newick_tree(final_sketches, final_mst, pre_sketch_by_file, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	if(is_phylip_tree){
		string output_phylip_file = output_file + ".phylip.tree";
		print_phylip_tree(final_sketches, final_mst, pre_sketch_by_file, output_phylip_file);
		cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
	}
	if(is_nexus_tree){
		string output_nexus_file = output_file + ".nexus.tree";
		print_nexus_tree(final_sketches, final_mst, pre_sketch_by_file, output_nexus_file);
		cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
	}
	if(is_linkage_matrix){
		string output_linkage_file = output_file + ".linkage.txt";
		print_linkage_matrix(final_sketches, final_mst, output_linkage_file);
		cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
	}

	vector<EdgeInfo> forest = generateForest(final_mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, final_sketches.size());
	printResult(tmpClust, final_sketches, pre_sketch_by_file, output_file, threshold);
	cerr << "-----write the cluster result into: " << output_file << endl;
	cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
	if(!no_save){
		saveMST(final_sketches, final_mst, new_folder_path, sketch_by_file);
	}

	if(!no_dense){
		loadANI(folder_path, pre_ani_arr, sketch_func_id_0);
		for(int i = 0; i < 101; i++)
			ani_arr[i] += pre_ani_arr[i];
		for(int i = 0; i < pre_dense_span; i++){
			for(int j = 0; j < pre_genome_number; j++){
				dense_arr[i][j] += pre_dense_arr[i][j];
			}
		}
		if(!no_save){
			saveANI(new_folder_path, ani_arr, sketch_func_id_0);
			saveDense(new_folder_path, dense_arr, dense_span, final_sketches.size());
		}

		int alpha = 2;
		int denseIndex = threshold / 0.01;
		vector<int> totalNoiseArr;
		for(int i = 0; i < tmpClust.size(); i++){
			if(tmpClust[i].size() == 1) continue;
			vector<PairInt> curDenseArr;
			set<int> denseSet;
			for(int j = 0; j < tmpClust[i].size(); j++){
				int element = tmpClust[i][j];
				PairInt p(element, dense_arr[denseIndex][element]);
				denseSet.insert(dense_arr[denseIndex][element]);
				curDenseArr.push_back(p);
			}
			vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
			totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
		}
		cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
		forest = modifyForest(forest, totalNoiseArr, threads);
		vector<vector<int>> cluster = generateClusterWithBfs(forest, final_sketches.size());
		string outputFileNew = output_file + ".removeNoise";
		printResult(cluster, final_sketches, pre_sketch_by_file, outputFileNew);
		cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
		cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
	}

}
void clust_from_mst_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads){
	vector<KssdSketchInfo> sketches;
	vector<EdgeInfo> mst;
	vector<vector<int>> cluster;
	bool sketchByFile = load_kssd_genome_info(folder_path, "mst", sketches);
	loadMST(folder_path, mst);

	if(is_newick_tree){
		string output_newick_file = outputFile + ".newick.tree";
		print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	if(is_phylip_tree){
		string output_phylip_file = outputFile + ".phylip.tree";
		print_kssd_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
		cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
	}
	if(is_nexus_tree){
		string output_nexus_file = outputFile + ".nexus.tree";
		print_kssd_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
		cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
	}
	if(is_linkage_matrix){
		string output_linkage_file = outputFile + ".linkage.txt";
		print_kssd_linkage_matrix(sketches, mst, output_linkage_file);
		cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
	}

	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printKssdResult(tmpClust, sketches, sketchByFile, outputFile, threshold);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;
	if(!no_dense){
		int **denseArr;
		int genome_number = sketches.size();
		int denseSpan = DENSE_SPAN;
		loadDense(denseArr, folder_path, denseSpan, genome_number);
		int alpha = 2;
		int denseIndex = threshold / 0.01;
		vector<int> totalNoiseArr;
		for(int i = 0; i < tmpClust.size(); i++){
			if(tmpClust[i].size() == 1) continue;
			vector<PairInt> curDenseArr;
			set<int> denseSet;
			for(int j = 0; j < tmpClust[i].size(); j++){
				int element = tmpClust[i][j];
				PairInt p(element, denseArr[denseIndex][element]);
				denseSet.insert(denseArr[denseIndex][element]);
				curDenseArr.push_back(p);
			}
			vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
			totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
		}
		cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
		forest = modifyForest(forest, totalNoiseArr, threads);
		cluster = generateClusterWithBfs(forest, sketches.size());
		string outputFileNew = outputFile + ".removeNoise";
		printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
		cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
		cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
	}
}

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads){
	vector<SketchInfo> sketches;
	vector<EdgeInfo> mst;
	vector<vector<int>> cluster;
	bool sketchByFile = load_genome_info(folder_path, "mst", sketches);
	loadMST(folder_path, mst);

	if(is_newick_tree){
		string output_newick_file = outputFile + ".newick.tree";
		print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	if(is_phylip_tree){
		string output_phylip_file = outputFile + ".phylip.tree";
		print_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
		cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
	}
	if(is_nexus_tree){
		string output_nexus_file = outputFile + ".nexus.tree";
		print_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
		cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
	}
	if(is_linkage_matrix){
		string output_linkage_file = outputFile + ".linkage.txt";
		print_linkage_matrix(sketches, mst, output_linkage_file);
		cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
	}
	
	// Automatic threshold selection based on MST edge length distribution
	if(is_auto_threshold){
		if(mst.empty()){
			cerr << "-----WARNING: MST is empty, cannot perform automatic threshold selection" << endl;
		} else if(mst.size() < 2){
			cerr << "-----WARNING: MST has only " << mst.size() << " edge(s), cannot perform automatic threshold selection" << endl;
		} else {
			cerr << "-----analyzing MST edge length distribution for automatic threshold selection..." << endl;
			if(is_stability){
				cerr << "-----stability evaluation enabled (this may take longer)..." << endl;
			}
			EdgeLengthStats stats = analyzeEdgeLengthDistribution(mst);
			// Use smaller min_gap_ratio (0.05 instead of 0.1) to detect more gaps
			// This helps find natural breakpoints even when distribution is dense
			vector<ThresholdCandidate> candidates = findThresholdCandidates(mst, 5, 0.05, is_stability, sketches.size());
			ThresholdCandidate optimal = selectOptimalThreshold(candidates, mst);
			
			string threshold_analysis_file = outputFile + ".threshold_analysis.txt";
			printThresholdAnalysis(mst, stats, candidates, optimal, threshold_analysis_file);
			
			cerr << "-----optimal threshold: " << optimal.threshold << " (confidence: " << optimal.confidence 
			     << ", suggested level: " << optimal.level << ")" << endl;
			if(is_stability){
				cerr << "-----stability: " << optimal.stability_score 
				     << " (split: " << optimal.stability_split << ", merge: " << optimal.stability_merge << ")" << endl;
				cerr << "-----near edges: " << optimal.near_edge_count << ", clusters: " << optimal.cluster_count << endl;
			}
			cerr << "-----threshold analysis written to: " << threshold_analysis_file << endl;
		}
	} else if(is_stability){
		// Evaluate stability for user-specified threshold (without auto-threshold)
		if(mst.empty()){
			cerr << "-----WARNING: MST is empty, cannot evaluate threshold stability" << endl;
		} else {
			cerr << "-----evaluating stability for threshold: " << threshold << "..." << endl;
			StabilityResult stability = computeThresholdStability(mst, threshold, sketches.size(), 0.01, 5, 100);
			vector<EdgeInfo> forest = generateForest(mst, threshold);
			vector<vector<int>> clusters = generateClusterWithBfs(forest, sketches.size());
			cerr << "-----threshold stability: " << stability.overall 
			     << " (split: " << stability.split << ", merge: " << stability.merge << ")" << endl;
			cerr << "-----near edges evaluated: " << stability.near_edge_count << ", clusters: " << clusters.size() << endl;
		}
	}

	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printResult(tmpClust, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

	if(!no_dense){
		int **denseArr;
		int genome_number = sketches.size();
		int denseSpan = DENSE_SPAN;
		loadDense(denseArr, folder_path, denseSpan, genome_number);
		int alpha = 2;
		int denseIndex = threshold / 0.01;
		vector<int> totalNoiseArr;
		for(int i = 0; i < tmpClust.size(); i++){
			if(tmpClust[i].size() == 1) continue;
			vector<PairInt> curDenseArr;
			set<int> denseSet;
			for(int j = 0; j < tmpClust[i].size(); j++){
				int element = tmpClust[i][j];
				PairInt p(element, denseArr[denseIndex][element]);
				denseSet.insert(denseArr[denseIndex][element]);
				curDenseArr.push_back(p);
			}
			vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
			totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
		}
		cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
		forest = modifyForest(forest, totalNoiseArr, threads);
		cluster = generateClusterWithBfs(forest, sketches.size());
		string outputFileNew = outputFile + ".removeNoise";
		printResult(cluster, sketches, sketchByFile, outputFileNew);
		cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
		cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
	}
}
#endif

void clust_from_genome_fast(const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index){
	bool isSave = !noSave;
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	KssdInvertedIndex inverted_index;  // Build index while generating sketches
	
	// Generate sketches and build inverted index simultaneously (pipeline optimization)
	compute_kssd_sketches_with_index(sketches, info, inverted_index, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
	
#ifndef GREEDY_CLUST	
	// Write index to files (index already built during sketch generation)
	transSketchesFromIndex(inverted_index, info, folder_path);
#else
#endif
	// Pass memory index to compute_kssd_clusters, which will pass it to compute_kssd_mst
	compute_kssd_clusters(sketches, info, sketchByFile, no_dense, isContainment, folder_path, outputFile, is_newick_tree, is_phylip_tree, is_nexus_tree, is_linkage_matrix, is_auto_threshold, is_stability, threshold, isSave, threads, dedup_dist, reps_per_cluster, save_rep_index, &inverted_index);

}

void compute_kssd_clusters(vector<KssdSketchInfo>& sketches, const KssdParameters info, bool sketchByFile, bool no_dense, bool isContainment, const string folder_path, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, double threshold, bool isSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index, KssdInvertedIndex* inverted_index){
	vector<vector<int>> cluster;
	double t2 = get_sec();

#ifdef GREEDY_CLUST
	//======clust-greedy====================================================================
	int sketch_func_id = 0;
	int kmer_size = info.half_k * 2; 
	
	if (save_rep_index && isSave) {
		KssdClusterState state = KssdInitialClusterWithState(sketches, info, threshold, threads, kmer_size);
		
		string state_file = folder_path + "/cluster_state.bin";
		state.save(state_file);
		cerr << "-----saved cluster state (with inverted index) to: " << state_file << endl;
		
		printKssdResult(state.clusters, sketches, sketchByFile, outputFile);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of " << outputFile << " is: " << state.clusters.size() << endl;
	} else {
		cluster = KssdGreedyClusterWithInvertedIndex(sketches, sketch_func_id, threshold, threads, kmer_size);
		printKssdResult(cluster, sketches, sketchByFile, outputFile);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
	}
	double t3 = get_sec();
#ifdef Timer
	cerr << "========time of greedyCluster is: " << t3 - t2 << "========" << endl;
#endif
	//======clust-greedy====================================================================
#else



	//======clust-mst=======================================================================
	int **denseArr;
	uint64_t* aniArr; //= new uint64_t[101];
	int denseSpan = DENSE_SPAN;
	// Use memory inverted index if available, otherwise load from file
	vector<EdgeInfo> mst = compute_kssd_mst(sketches, info, folder_path, 0, no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold, inverted_index);
	double t3 = get_sec();
#ifdef Timer
	cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
#endif
	//string sketch_func_id = "fast_kssd_sketch";
	if(isSave){
		if(!no_dense){
			saveANI(folder_path, aniArr, 1);
			saveDense(folder_path, denseArr, denseSpan, sketches.size());
		}
		saveKssdMST(sketches, mst, folder_path, sketchByFile);
	}
	double t4 = get_sec();
#ifdef Timer
	cerr << "========time of saveMST is: " << t4 - t3 << "========" << endl;
#endif

	//generate the Newick tree format
	if(is_newick_tree){
		string output_newick_file = outputFile + ".newick.tree";
		print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	if(is_phylip_tree){
		string output_phylip_file = outputFile + ".phylip.tree";
		print_kssd_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
		cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
	}
	if(is_nexus_tree){
		string output_nexus_file = outputFile + ".nexus.tree";
		print_kssd_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
		cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
	}
	if(is_linkage_matrix){
		string output_linkage_file = outputFile + ".linkage.txt";
		print_kssd_linkage_matrix(sketches, mst, output_linkage_file);
		cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
	}
	
	// Automatic threshold selection based on MST edge length distribution
	if(is_auto_threshold){
		if(mst.empty()){
			cerr << "-----WARNING: MST is empty, cannot perform automatic threshold selection" << endl;
		} else if(mst.size() < 2){
			cerr << "-----WARNING: MST has only " << mst.size() << " edge(s), cannot perform automatic threshold selection" << endl;
		} else {
			cerr << "-----analyzing MST edge length distribution for automatic threshold selection..." << endl;
			if(is_stability){
				cerr << "-----stability evaluation enabled (this may take longer)..." << endl;
			}
			EdgeLengthStats stats = analyzeEdgeLengthDistribution(mst);
			// Use smaller min_gap_ratio (0.05 instead of 0.1) to detect more gaps
			// This helps find natural breakpoints even when distribution is dense
			vector<ThresholdCandidate> candidates = findThresholdCandidates(mst, 5, 0.05, is_stability, sketches.size());
			ThresholdCandidate optimal = selectOptimalThreshold(candidates, mst);
			
			string threshold_analysis_file = outputFile + ".threshold_analysis.txt";
			printThresholdAnalysis(mst, stats, candidates, optimal, threshold_analysis_file);
			
			cerr << "-----optimal threshold: " << optimal.threshold << " (confidence: " << optimal.confidence 
			     << ", suggested level: " << optimal.level << ")" << endl;
			if(is_stability){
				cerr << "-----stability: " << optimal.stability_score 
				     << " (split: " << optimal.stability_split << ", merge: " << optimal.stability_merge << ")" << endl;
				cerr << "-----near edges: " << optimal.near_edge_count << ", clusters: " << optimal.cluster_count << endl;
			}
			cerr << "-----threshold analysis written to: " << threshold_analysis_file << endl;
		}
	} else if(is_stability){
		// Evaluate stability for user-specified threshold (without auto-threshold)
		if(mst.empty()){
			cerr << "-----WARNING: MST is empty, cannot evaluate threshold stability" << endl;
		} else {
			cerr << "-----evaluating stability for threshold: " << threshold << "..." << endl;
			StabilityResult stability = computeThresholdStability(mst, threshold, sketches.size(), 0.01, 5, 100);
			vector<EdgeInfo> forest = generateForest(mst, threshold);
			vector<vector<int>> clusters = generateClusterWithBfs(forest, sketches.size());
			cerr << "-----threshold stability: " << stability.overall 
			     << " (split: " << stability.split << ", merge: " << stability.merge << ")" << endl;
			cerr << "-----near edges evaluated: " << stability.near_edge_count << ", clusters: " << clusters.size() << endl;
		}
	}

	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printKssdResult(tmpClust, sketches, sketchByFile, outputFile, threshold);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

	if(dedup_dist > 0 || reps_per_cluster > 0){
		vector<int> node_to_rep;
		vector<vector<int>> candidates = build_dedup_candidates_per_cluster(tmpClust, forest, sketches, sketchByFile, dedup_dist, node_to_rep);
		if(dedup_dist > 0){
			string outputDedup = outputFile + ".dedup";
			printKssdResult(candidates, sketches, sketchByFile, outputDedup);
			cerr << "-----write the deduped cluster result into: " << outputDedup << endl;
		}
		if(reps_per_cluster > 0){
			vector<vector<int>> reps = select_k_reps_per_cluster_tree(tmpClust, candidates, forest, sketches.size(), node_to_rep, reps_per_cluster);
			string outputReps = outputFile + ".reps";
			printKssdResult(reps, sketches, sketchByFile, outputReps);
			cerr << "-----write the reps-per-cluster result into: " << outputReps << endl;
		}
	}
	//tune cluster by noise cluster
	if(!no_dense){
		int alpha = 2;
		int denseIndex = threshold / 0.01;
		vector<int> totalNoiseArr;
		for(int i = 0; i < tmpClust.size(); i++){
			if(tmpClust[i].size() == 1) continue;
			vector<PairInt> curDenseArr;
			set<int> denseSet;
			for(int j = 0; j < tmpClust[i].size(); j++){
				int element = tmpClust[i][j];
				PairInt p(element, denseArr[denseIndex][element]);
				denseSet.insert(denseArr[denseIndex][element]);
				curDenseArr.push_back(p);
			}
			vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
			totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
		}
		cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
		forest = modifyForest(forest, totalNoiseArr, threads);
		cluster = generateClusterWithBfs(forest, sketches.size());
		string outputFileNew = outputFile + ".removeNoise";
		printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
		cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
		cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;

		if(dedup_dist > 0 || reps_per_cluster > 0){
			vector<int> node_to_rep;
			vector<vector<int>> candidates = build_dedup_candidates_per_cluster(cluster, forest, sketches, sketchByFile, dedup_dist, node_to_rep);
			if(dedup_dist > 0){
				string outputDedup = outputFileNew + ".dedup";
				printKssdResult(candidates, sketches, sketchByFile, outputDedup);
				cerr << "-----write the deduped cluster result into: " << outputDedup << endl;
			}
			if(reps_per_cluster > 0){
				vector<vector<int>> reps = select_k_reps_per_cluster_tree(cluster, candidates, forest, sketches.size(), node_to_rep, reps_per_cluster);
				string outputReps = outputFileNew + ".reps";
				printKssdResult(reps, sketches, sketchByFile, outputReps);
				cerr << "-----write the reps-per-cluster result into: " << outputReps << endl;
			}
		}
		double t5 = get_sec();
#ifdef Timer
		cerr << "========time of tuning cluster is: " << t5 - t4 << "========" << endl;
#endif
	}
	//======clust-mst=======================================================================
#endif//endif GREEDY_CLUST
}

void compute_kssd_sketches(vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads){
	double t0 = get_sec();
	if(sketchByFile){
		//if(!sketchFiles(inputFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){
		if(!sketchFileWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads)){
			cerr << "ERROR: sketchFileWithKssd(), cannot finish the sketch generation by genome files" << endl;
			exit(1);
		}
	}//end sketch by sequence
		else{
			cerr << "use the sketch seuqnce with kssd " << endl;
			if(!sketchSequencesWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads)){
				cerr << "ERROR: sketchSequencesWithKssd (), cannot finish the sketch generation by genome sequences" << endl;
				exit(1);
			}
		}//end sketch by file
		cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;
		double t1 = get_sec();
#ifdef Timer
		cerr << "========time of computing sketch is: " << t1 - t0 << "========" << endl;
#endif
		if(folder_path.empty()){
			folder_path = currentDataTime();
		}
		if(isSave){
			string command = "mkdir -p " + folder_path;
			system(command.c_str());
			saveKssdSketches(sketches, info, folder_path, sketchByFile);
			double t2 = get_sec();
#ifdef Timer
			cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
		}

	}

// Pipeline version: generate sketches and build inverted index simultaneously
void compute_kssd_sketches_with_index(vector<KssdSketchInfo>& sketches, KssdParameters& info, KssdInvertedIndex& inverted_index, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads){
	double t0 = get_sec();
	if(sketchByFile){
		if(!sketchFileWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads, &inverted_index)){
			cerr << "ERROR: sketchFileWithKssd(), cannot finish the sketch generation by genome files" << endl;
			exit(1);
		}
	} else {
		cerr << "use the sketch sequence with kssd " << endl;
		if(!sketchSequencesWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads, &inverted_index)){
			cerr << "ERROR: sketchSequencesWithKssd (), cannot finish the sketch generation by genome sequences" << endl;
			exit(1);
		}
	}
	cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;
	double t1 = get_sec();
#ifdef Timer
	cerr << "========time of computing sketch (with index) is: " << t1 - t0 << "========" << endl;
#endif
	if(folder_path.empty()){
		folder_path = currentDataTime();
	}
	if(isSave){
		string command = "mkdir -p " + folder_path;
		system(command.c_str());
		saveKssdSketches(sketches, info, folder_path, sketchByFile);
		double t2 = get_sec();
#ifdef Timer
		cerr << "========time of saveKssdSketches is: " << t2 - t1 << "========" << endl;
#endif
	}
}

	static bool looks_like_cluster_result_file(const string& filePath){
		std::ifstream in(filePath);
		if(!in.good()) return false;
		string line;
		while(std::getline(in, line)){
			bool allSpace = true;
			for(char c : line){
				if(!isspace((unsigned char)c)){ allSpace = false; break; }
			}
			if(allSpace) continue;
			if(line.rfind("the cluster", 0) == 0) return true;
			return false;
		}
		return false;
	}

	static string materialize_file_list_from_cluster(const string& clusterFile, const string& outListFile){
		std::ifstream in(clusterFile);
		if(!in.good()){
			cerr << "ERROR: build_kssd_db_fast(), cannot open cluster file: " << clusterFile << endl;
			exit(1);
		}
		std::ofstream out(outListFile);
		if(!out.good()){
			cerr << "ERROR: build_kssd_db_fast(), cannot write list file: " << outListFile << endl;
			exit(1);
		}

		std::unordered_set<string> seen;
		string line;
		while(std::getline(in, line)){
			if(line.empty()) continue;
			if(line.rfind("the cluster", 0) == 0) continue;
			if(!(line[0] == '\t' || line[0] == ' ')) continue;

			std::istringstream iss(line);
			string localIdx, globalId, lenNt, filePath;
			if(!(iss >> localIdx >> globalId >> lenNt >> filePath)) continue;
			if(seen.insert(filePath).second){
				out << filePath << "\n";
			}
		}
		return outListFile;
	}

	void build_kssd_db_fast(const string input_file, const string db_folder, bool isSetKmer, bool& isContainment, int minLen, int& kmerSize, int& drlevel, int threads){
		string command = "mkdir -p " + db_folder;
		system(command.c_str());

		string listFile = input_file;
		if(looks_like_cluster_result_file(input_file)){
			string outList = db_folder + "/builddb.list";
			listFile = materialize_file_list_from_cluster(input_file, outList);
			cerr << "-----buildDB: extracted genome paths from cluster file into: " << listFile << endl;
		}else{
			cerr << "-----buildDB: using input as genome file list: " << listFile << endl;
		}

		// Tune kmerSize/drlevel using the materialized list (cluster files are not compatible with tuning directly)
		{
			double dummyThreshold = 0.0;
			bool sketchByFile = true;
			if(!tune_kssd_parameters(sketchByFile, isSetKmer, listFile, threads, minLen, isContainment, kmerSize, dummyThreshold, drlevel)){
				cerr << "ERROR: build_kssd_db_fast(), failed to tune KSSD parameters" << endl;
				exit(1);
			}
		}

		vector<KssdSketchInfo> sketches;
		KssdParameters info;
		string folder_path = db_folder;
		bool sketchByFile = true;
		bool isSave = true;
		compute_kssd_sketches(sketches, info, isSave, listFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
		transSketches(sketches, info, folder_path, threads);
		cerr << "-----buildDB: finished building KSSD DB at: " << folder_path << endl;
	}

	void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool sketchByFile, bool no_dense, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads, bool use_inverted_index, bool save_rep_index){
		bool isSave = !noSave;
		vector<SketchInfo> sketches;
		int sketch_func_id;
		if(sketchFunc == "MinHash")	sketch_func_id = 0;
		else if(sketchFunc == "KSSD") sketch_func_id = 1;

		MinHashInvertedIndex minhash_index;
		MinHashInvertedIndex* idx_ptr = (sketch_func_id == 0) ? &minhash_index : nullptr;

		compute_sketches(sketches, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, isSave, threads, idx_ptr);

		compute_clusters(sketches, sketchByFile, outputFile, is_newick_tree, is_phylip_tree, is_nexus_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, folder_path, sketch_func_id, threshold, isSave, threads, use_inverted_index, save_rep_index, idx_ptr);
	}

	bool tune_kssd_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, int& kmerSize, double& threshold, int &drlevel){
		uint64_t maxSize, minSize, averageSize;
		calSize(sketchByFile, inputFile, threads, minLen, maxSize, minSize, averageSize);

		//======tune the sketch_size===============
		int compression = 1 << (4 * drlevel);
		int sketchSize = averageSize / compression;
		//=====tune the kmer_size===============
		double warning_rate = 0.01;
		double recommend_rate = 0.0001;
		int alphabetSize = 4;//for "AGCT"
		int recommendedKmerSize = ceil(log(maxSize * (1 - recommend_rate) / recommend_rate) / log(4));
		int warningKmerSize = ceil(log(maxSize * (1 - warning_rate) / warning_rate) / log(4));
		if(!isSetKmer){
			kmerSize = recommendedKmerSize;
		}
		else{
			if(kmerSize < warningKmerSize){
				cerr << "the kmerSize " << kmerSize << " is too small for the maximum genome size of " << maxSize << endl;
				cerr << "replace the kmerSize to the: " << recommendedKmerSize << " for reducing the random collision of kmers" << endl;
				kmerSize = recommendedKmerSize;
			}
			else if(kmerSize > recommendedKmerSize + 3){
				cerr << "the kmerSize " << kmerSize << " maybe too large for the maximum genome size of " << maxSize << endl;
				cerr << "replace the kmerSize to the " << recommendedKmerSize << " for increasing the sensitivity of genome comparison" << endl;
				kmerSize = recommendedKmerSize;
			}
		}

		//=====tune the distance threshold===============
		double minJaccard = 0.001;
		if(!isContainment){
			minJaccard = 1.0 / sketchSize;
		}
		else{
			//minJaccard = 1.0 / (averageSize / containCompress);
			minJaccard = 1.0 / (minSize / compression);
		}

		double maxDist;
		if(minJaccard >= 1.0)	
			maxDist = 1.0;
		else
			maxDist = -1.0 / kmerSize * log(2*minJaccard / (1.0 + minJaccard));
		cerr << "-----the max recommand distance threshold is: " << maxDist << endl;
		if(threshold > maxDist){
			cerr << "ERROR: tune_parameters(), the threshold: " << threshold << " is out of the valid distance range estimated by Mash distance or AAF distance" << endl;
			cerr << "Please set a distance threshold with -d option" << endl;
			return false;
		}

#ifdef DEBUG
		if(sketchByFile) cerr << "-----sketch by file!" << endl;
		else cerr << "-----sketch by sequence!" << endl;
		cerr << "-----the kmerSize is: " << kmerSize << endl;
		cerr << "-----the thread number is: " << threads << endl;
		cerr << "-----the threshold is: " << threshold << endl;
		if(isContainment)
			cerr << "-----use the AAF distance, the sketchSize is approximately in proportion with 1/" << compression << endl;
		else
			cerr << "-----use the Mash distance, the sketchSize is about: " << sketchSize << endl;
#endif

		return true;
	}

	bool tune_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, bool& isJaccard, int& kmerSize, double& threshold, int& containCompress, int& sketchSize){
		uint64_t maxSize, minSize, averageSize;
		calSize(sketchByFile, inputFile, threads, minLen, maxSize, minSize, averageSize);

		//======tune the sketch_size===============
		if(isContainment && isJaccard){
			cerr << "ERROR: tune_parameters(), conflict distance measurement of Mash distance (fixed-sketch-size) and AAF distance (variable-sketch-size) " << endl;
			return false;
		}
#ifdef GREEDY_CLUST
		//======clust-greedy====================================================================
		if(!isContainment && !isJaccard){
			containCompress = averageSize / 1000;
			isContainment = true;
		}
		else if(!isContainment && isJaccard){
			//do nothing
		}
		else{
			if(averageSize / containCompress < 10){
				cerr << "the containCompress " << containCompress << " is too large and the sketch size is too small" << endl;
				containCompress = averageSize / 1000;
				cerr << "set the containCompress to: " << containCompress << endl;
			}
		}
		//=======clust-greedy===================================================================
#endif
		//=====tune the kmer_size===============
		double warning_rate = 0.01;
		double recommend_rate = 0.0001;
		int alphabetSize = 4;//for "AGCT"
		int recommendedKmerSize = ceil(log(maxSize * (1 - recommend_rate) / recommend_rate) / log(4));
		int warningKmerSize = ceil(log(maxSize * (1 - warning_rate) / warning_rate) / log(4));
		if(!isSetKmer){
			kmerSize = recommendedKmerSize;
		}
		else{
			if(kmerSize < warningKmerSize){
				cerr << "the kmerSize " << kmerSize << " is too small for the maximum genome size of " << maxSize << endl;
				cerr << "replace the kmerSize to the: " << recommendedKmerSize << " for reducing the random collision of kmers" << endl;
				kmerSize = recommendedKmerSize;
			}
			else if(kmerSize > recommendedKmerSize + 3){
				cerr << "the kmerSize " << kmerSize << " maybe too large for the maximum genome size of " << maxSize << endl;
				cerr << "replace the kmerSize to the " << recommendedKmerSize << " for increasing the sensitivity of genome comparison" << endl;
				kmerSize = recommendedKmerSize;
			}
		}

		//=====tune the distance threshold===============
		double minJaccard = 0.001;
		if(!isContainment){
			minJaccard = 1.0 / sketchSize;
		}
		else{
			//minJaccard = 1.0 / (averageSize / containCompress);
			minJaccard = 1.0 / (minSize / containCompress);
		}

		double maxDist;
		if(minJaccard >= 1.0)	
			maxDist = 1.0;
		else
			maxDist = -1.0 / kmerSize * log(2*minJaccard / (1.0 + minJaccard));
		cerr << "-----the max recommand distance threshold is: " << maxDist << endl;
		if(threshold > maxDist){
			cerr << "ERROR: tune_parameters(), the threshold: " << threshold << " is out of the valid distance range estimated by Mash distance or AAF distance" << endl;
			cerr << "Please set a distance threshold with -d option" << endl;
			return false;
		}

#ifdef DEBUG
		if(sketchByFile) cerr << "-----sketch by file!" << endl;
		else cerr << "-----sketch by sequence!" << endl;
		cerr << "-----the kmerSize is: " << kmerSize << endl;
		cerr << "-----the thread number is: " << threads << endl;
		cerr << "-----the threshold is: " << threshold << endl;
		if(isContainment)
			cerr << "-----use the AAF distance (variable-sketch-size), the sketchSize is in proportion with 1/" << containCompress << endl;
		else
			cerr << "-----use the Mash distance (fixed-sketch-size), the sketchSize is: " << sketchSize << endl;
#endif

		return true;
	}

	void clust_from_sketch_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool no_dense, bool isContainment, double threshold, int threads, double dedup_dist, int reps_per_cluster, bool use_inverted_index, bool save_rep_index){
		vector<KssdSketchInfo> sketches;
		vector<vector<int>> cluster;
		bool sketchByFile;
		KssdParameters info;
#ifdef GREEDY_CLUST
		//======clust-greedy====================================================================
		double time0 = get_sec();
		sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
		cerr << "-----the size of sketches is: " << sketches.size() << endl;
		double time1 = get_sec();
#ifdef Timer
		cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
		int kmer_size = info.half_k * 2;  
		
		if (save_rep_index) {
			KssdClusterState state = KssdInitialClusterWithState(sketches, info, threshold, threads, kmer_size);
			
			string state_file = folder_path + "/cluster_state.bin";
			state.save(state_file);
			cerr << "-----saved cluster state (with inverted index) to: " << state_file << endl;
			
			printKssdResult(state.clusters, sketches, sketchByFile, outputFile);
			cerr << "-----write the cluster result into: " << outputFile << endl;
			cerr << "-----the cluster number of " << outputFile << " is: " << state.clusters.size() << endl;
		} else {
			cluster = KssdGreedyClusterWithInvertedIndex(sketches, 0, threshold, threads, kmer_size);
			printKssdResult(cluster, sketches, sketchByFile, outputFile);
			cerr << "-----write the cluster result into: " << outputFile << endl;
			cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
		}
		double time2 = get_sec();
#ifdef Timer
		cerr << "========time of greedy incremental cluster is: " << time2 - time1 << endl;
#endif
		//======clust-greedy====================================================================
#else     


		//======clust-mst=======================================================================
		double time0 = get_sec();
		//sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
		sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);

		cerr << "-----the size of sketches is: " << sketches.size() << endl;
		double time1 = get_sec();
#ifdef Timer
		cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
		int** denseArr;
		uint64_t* aniArr; //= new uint64_t[101];
		int denseSpan = DENSE_SPAN;
		//vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, denseArr, denseSpan, aniArr);
		vector<EdgeInfo> mst = compute_kssd_mst(sketches, info, folder_path, 0, no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold);
		double time2 = get_sec();
#ifdef Timer
		cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
#endif
		if(is_newick_tree){
			string output_newick_file = outputFile + ".newick.tree";
			print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
			cerr << "-----write the newick tree into: " << output_newick_file << endl;
		}
		if(is_phylip_tree){
			string output_phylip_file = outputFile + ".phylip.tree";
			print_kssd_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
			cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
		}
		if(is_nexus_tree){
			string output_nexus_file = outputFile + ".nexus.tree";
			print_kssd_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
			cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
		}
		if(is_linkage_matrix){
			string output_linkage_file = outputFile + ".linkage.txt";
			print_kssd_linkage_matrix(sketches, mst, output_linkage_file);
			cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
		}
		
		// Automatic threshold selection based on MST edge length distribution
		if(is_auto_threshold){
			if(mst.empty()){
				cerr << "-----WARNING: MST is empty, cannot perform automatic threshold selection" << endl;
			} else if(mst.size() < 2){
				cerr << "-----WARNING: MST has only " << mst.size() << " edge(s), cannot perform automatic threshold selection" << endl;
			} else {
				cerr << "-----analyzing MST edge length distribution for automatic threshold selection..." << endl;
				EdgeLengthStats stats = analyzeEdgeLengthDistribution(mst);
				vector<ThresholdCandidate> candidates = findThresholdCandidates(mst, 5, 0.1, false, sketches.size());
				ThresholdCandidate optimal = selectOptimalThreshold(candidates, mst);
				
				string threshold_analysis_file = outputFile + ".threshold_analysis.txt";
				printThresholdAnalysis(mst, stats, candidates, optimal, threshold_analysis_file);
				
				cerr << "-----optimal threshold: " << optimal.threshold << " (confidence: " << optimal.confidence 
				     << ", suggested level: " << optimal.level << ")" << endl;
				cerr << "-----threshold analysis written to: " << threshold_analysis_file << endl;
			}
		}
		
		vector<EdgeInfo> forest = generateForest(mst, threshold);
		vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
		printKssdResult(tmpClust, sketches, sketchByFile, outputFile, threshold);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

		if(dedup_dist > 0 || reps_per_cluster > 0){
			vector<int> node_to_rep;
			vector<vector<int>> candidates = build_dedup_candidates_per_cluster(tmpClust, forest, sketches, sketchByFile, dedup_dist, node_to_rep);
			if(dedup_dist > 0){
				string outputDedup = outputFile + ".dedup";
				printKssdResult(candidates, sketches, sketchByFile, outputDedup);
				cerr << "-----write the deduped cluster result into: " << outputDedup << endl;
			}
			if(reps_per_cluster > 0){
				vector<vector<int>> reps = select_k_reps_per_cluster_tree(tmpClust, candidates, forest, sketches.size(), node_to_rep, reps_per_cluster);
				string outputReps = outputFile + ".reps";
				printKssdResult(reps, sketches, sketchByFile, outputReps);
				cerr << "-----write the reps-per-cluster result into: " << outputReps << endl;
			}
		}

		if(!no_dense){
			int alpha = 2;
			int denseIndex = threshold / 0.01;
			vector<int> totalNoiseArr;
			for(int i = 0; i < tmpClust.size(); i++){
				if(tmpClust[i].size() == 1) continue;
				vector<PairInt> curDenseArr;
				set<int> denseSet;
				for(int j = 0; j < tmpClust[i].size(); j++){
					int element = tmpClust[i][j];
					PairInt p(element, denseArr[denseIndex][element]);
					denseSet.insert(denseArr[denseIndex][element]);
					curDenseArr.push_back(p);
				}
				vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
				totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
			}
			cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
			forest = modifyForest(forest, totalNoiseArr, threads);
			cluster = generateClusterWithBfs(forest, sketches.size());
			string outputFileNew = outputFile + ".removeNoise";
			printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
			cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
			cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;

			if(dedup_dist > 0 || reps_per_cluster > 0){
				vector<int> node_to_rep;
				vector<vector<int>> candidates = build_dedup_candidates_per_cluster(cluster, forest, sketches, sketchByFile, dedup_dist, node_to_rep);
				if(dedup_dist > 0){
					string outputDedup = outputFileNew + ".dedup";
					printKssdResult(candidates, sketches, sketchByFile, outputDedup);
					cerr << "-----write the deduped cluster result into: " << outputDedup << endl;
				}
				if(reps_per_cluster > 0){
					vector<vector<int>> reps = select_k_reps_per_cluster_tree(cluster, candidates, forest, sketches.size(), node_to_rep, reps_per_cluster);
					string outputReps = outputFileNew + ".reps";
					printKssdResult(reps, sketches, sketchByFile, outputReps);
					cerr << "-----write the reps-per-cluster result into: " << outputReps << endl;
				}
			}
		}
		double time3 = get_sec();
#ifdef Timer
		cerr << "========time of generator forest and cluster is: " << time3 - time2 << "========" << endl;
#endif
#endif
}

	void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads, bool use_inverted_index, bool save_rep_index){
		vector<SketchInfo> sketches;
		vector<vector<int>> cluster;
		int sketch_func_id;
		bool sketchByFile;
#ifdef GREEDY_CLUST
		//======clust-greedy====================================================================
		double time0 = get_sec();
		sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
		cerr << "-----the size of sketches is: " << sketches.size() << endl;
		if(sketchByFile)
			sort(sketches.begin(), sketches.end(), cmpGenomeSize);
		else
			sort(sketches.begin(), sketches.end(), cmpSeqSize);
		double time1 = get_sec();
#ifdef Timer
		cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
		// Use inverted index optimization for MinHash if requested
		if (sketch_func_id == 0 && use_inverted_index) {
			// Get actual kmer_size from first sketch
			int kmer_size = sketches[0].minHash->getKmerSize();
			int sketch_size = sketches[0].minHash->getSketchSize();
			bool is_containment = sketches[0].isContainment;
			
			// Save state if requested
			if (save_rep_index) {
				MinHashClusterState state = MinHashInitialClusterWithState(sketches, threshold, threads, kmer_size, sketch_size, is_containment);
				string state_file = folder_path + "/cluster_state.bin";
				state.save(state_file);
				cerr << "-----saved cluster state (with inverted index) to: " << state_file << endl;
				cluster = state.clusters;
			} else {
				cluster = MinHashGreedyClusterWithInvertedIndex(sketches, sketch_func_id, threshold, threads, kmer_size);
			}
		} else {
			cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
		}
		printResult(cluster, sketches, sketchByFile, outputFile);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
		double time2 = get_sec();
#ifdef Timer
		cerr << "========time of greedy incremental cluster is: " << time2 - time1 << endl;
#endif
		//======clust-greedy====================================================================
#else
		//======clust-mst=======================================================================
		double time0 = get_sec();
		sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
		cerr << "-----the size of sketches is: " << sketches.size() << endl;
		double time1 = get_sec();
#ifdef Timer
		cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
		int** denseArr;
		uint64_t* aniArr; //= new uint64_t[101];
		int denseSpan = DENSE_SPAN;
		
		// Use inverted index optimization for MinHash MST
		vector<EdgeInfo> mst;
		if(sketch_func_id == 0 && sketches.size() > 0 && sketches[0].minHash != nullptr) {
			// Read kmerSize from sketch parameters file (more reliable than from sketch object)
			int kmer_size_from_file = 21;  // Default for MinHash
			int sketch_func_id_check, contain_compress, sketch_size, half_k, half_subk, drlevel;
			bool is_containment_from_file;
			read_sketch_parameters(folder_path, sketch_func_id_check, kmer_size_from_file, is_containment_from_file, contain_compress, sketch_size, half_k, half_subk, drlevel);

			MinHashInvertedIndex inverted_index;

			// Try to load pre-built index from sketch folder; rebuild only if not found
			if (!inverted_index.load_from_file(folder_path)) {
				cerr << "-----building MinHash inverted index for MST computation..." << endl;

				vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);

				#pragma omp parallel num_threads(threads)
				{
					int tid = omp_get_thread_num();
					auto& local_index = thread_local_indices[tid];

					#pragma omp for schedule(dynamic, 100)
					for(size_t i = 0; i < sketches.size(); i++){
						if(sketches[i].minHash != nullptr) {
							sketches[i].id = i;
							vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
							for(uint64_t hash : hashes)
								local_index[hash].push_back(i);
						}
					}
				}

				cerr << "-----merging thread-local indices..." << endl;
				for(int tid = 0; tid < threads; tid++){
					for(auto& entry : thread_local_indices[tid]){
						inverted_index.hash_map[entry.first].insert(
							inverted_index.hash_map[entry.first].end(),
							entry.second.begin(), entry.second.end()
						);
					}
				}
				cerr << "-----MinHash inverted index built: " << inverted_index.hash_map.size() << " unique hashes" << endl;
				// Save for future reuse
				inverted_index.save_to_file(folder_path);
			}

			int kmer_size = kmer_size_from_file;
			bool is_containment = is_containment_from_file;
			if(kmer_size <= 0) {
				kmer_size = sketches[0].minHash->getKmerSize();
				is_containment = sketches[0].isContainment;
			}

			// Use compute_minhash_mst with inverted index
			mst = compute_minhash_mst(sketches, 0, no_dense, is_containment, threads, denseArr, denseSpan, aniArr, threshold, kmer_size, &inverted_index);
		} else {
			// Fallback to traditional modifyMST for non-MinHash or when MinHash is not available
			mst = modifyMST(sketches, 0, sketch_func_id, threads, no_dense, denseArr, denseSpan, aniArr);
		}
		double time2 = get_sec();
#ifdef Timer
		cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
#endif
		if(is_newick_tree){
			string output_newick_file = outputFile + ".newick.tree";
			print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
			cerr << "-----write the newick tree into: " << output_newick_file << endl;
		}
		if(is_phylip_tree){
			string output_phylip_file = outputFile + ".phylip.tree";
			print_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
			cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
		}
		if(is_nexus_tree){
			string output_nexus_file = outputFile + ".nexus.tree";
			print_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
			cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
		}
		
		// Automatic threshold selection based on MST edge length distribution
		if(is_auto_threshold){
			if(mst.empty()){
				cerr << "-----WARNING: MST is empty, cannot perform automatic threshold selection" << endl;
			} else if(mst.size() < 2){
				cerr << "-----WARNING: MST has only " << mst.size() << " edge(s), cannot perform automatic threshold selection" << endl;
			} else {
				cerr << "-----analyzing MST edge length distribution for automatic threshold selection..." << endl;
				EdgeLengthStats stats = analyzeEdgeLengthDistribution(mst);
				vector<ThresholdCandidate> candidates = findThresholdCandidates(mst, 5, 0.1, false, sketches.size());
				ThresholdCandidate optimal = selectOptimalThreshold(candidates, mst);
				
				string threshold_analysis_file = outputFile + ".threshold_analysis.txt";
				printThresholdAnalysis(mst, stats, candidates, optimal, threshold_analysis_file);
				
				cerr << "-----optimal threshold: " << optimal.threshold << " (confidence: " << optimal.confidence 
				     << ", suggested level: " << optimal.level << ")" << endl;
				cerr << "-----threshold analysis written to: " << threshold_analysis_file << endl;
			}
		}
		
		vector<EdgeInfo> forest = generateForest(mst, threshold);
		vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
		printResult(tmpClust, sketches, sketchByFile, outputFile, threshold);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

		if(!no_dense){
			int alpha = 2;
			int denseIndex = threshold / 0.01;
			vector<int> totalNoiseArr;
			for(int i = 0; i < tmpClust.size(); i++){
				if(tmpClust[i].size() == 1) continue;
				vector<PairInt> curDenseArr;
				set<int> denseSet;
				for(int j = 0; j < tmpClust[i].size(); j++){
					int element = tmpClust[i][j];
					PairInt p(element, denseArr[denseIndex][element]);
					denseSet.insert(denseArr[denseIndex][element]);
					curDenseArr.push_back(p);
				}
				vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
				totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
			}
			cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
			forest = modifyForest(forest, totalNoiseArr, threads);
			cluster = generateClusterWithBfs(forest, sketches.size());
			string outputFileNew = outputFile + ".removeNoise";
			printResult(cluster, sketches, sketchByFile, outputFileNew);
			cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
			cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
		}
		double time3 = get_sec();
#ifdef Timer
		cerr << "========time of generator forest and cluster is: " << time3 - time2 << "========" << endl;
#endif
		//=======clust-mst======================================================================
#endif
	}

	void compute_sketches(vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress,  bool isSave, int threads, MinHashInvertedIndex* inverted_index){
		double t0 = get_sec();
		if(sketchByFile){
			if(!sketchFiles(inputFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads, inverted_index)){
				cerr << "ERROR: generate_sketches(), cannot finish the sketch generation by genome files" << endl;
				exit(1);
			}
		}//end sketch by sequence
		else{
			if(!sketchSequences(inputFile, kmerSize, sketchSize, minLen, sketchFunc, isContainment, containCompress, sketches, threads, inverted_index)){
				cerr << "ERROR: generate_sketches(), cannot finish the sketch generation by genome sequences" << endl;
				exit(1);
			}
		}//end sketch by file
		cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;
		double t1 = get_sec();
#ifdef Timer
		cerr << "========time of computing sketch is: " << t1 - t0 << "========" << endl;
#endif
		folder_path = currentDataTime();
		if(isSave){
			string command = "mkdir -p " + folder_path;
			system(command.c_str());
			saveSketches(sketches, folder_path, sketchByFile, sketchFunc, isContainment, containCompress, sketchSize, kmerSize);
			if (inverted_index != nullptr && !inverted_index->hash_map.empty())
				inverted_index->save_to_file(folder_path);
			double t2 = get_sec();
#ifdef Timer
			cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
		}
	}

	void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, bool is_phylip_tree, bool is_nexus_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads, bool use_inverted_index, bool save_rep_index, MinHashInvertedIndex* inverted_index){
		vector<vector<int>> cluster;
		double t2 = get_sec();
#ifdef GREEDY_CLUST
		//======clust-greedy====================================================================
		// Use inverted index optimization for MinHash if requested
		if (sketch_func_id == 0 && use_inverted_index) {
			// Get kmer_size from first sketch
			int kmer_size = sketches[0].minHash->getKmerSize();
			int sketch_size = sketches[0].minHash->getSketchSize();
			bool is_containment = sketches[0].isContainment;
			
			// Save state if requested
			if (save_rep_index && isSave) {
				MinHashClusterState state = MinHashInitialClusterWithState(sketches, threshold, threads, kmer_size, sketch_size, is_containment);
				string state_file = folder_path + "/cluster_state.bin";
				state.save(state_file);
				cerr << "-----saved cluster state (with inverted index) to: " << state_file << endl;
				cluster = state.clusters;
			} else {
				cluster = MinHashGreedyClusterWithInvertedIndex(sketches, sketch_func_id, threshold, threads, kmer_size);
			}
		} else {
			cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
		}
		printResult(cluster, sketches, sketchByFile, outputFile);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
		double t3 = get_sec();
#ifdef Timer
		cerr << "========time of greedyCluster is: " << t3 - t2 << "========" << endl;
#endif
		//======clust-greedy====================================================================
#else
		//======clust-mst=======================================================================
		int **denseArr;
		uint64_t* aniArr; //= new uint64_t[101];
		int denseSpan = DENSE_SPAN;
		
		// Use inverted index optimization for MinHash MST
		vector<EdgeInfo> mst;
		if(sketch_func_id == 0 && sketches.size() > 0 && sketches[0].minHash != nullptr) {
			int kmer_size_from_file = 21;
			int sketch_func_id_check, contain_compress, sketch_size, half_k, half_subk, drlevel;
			bool is_containment_from_file;
			read_sketch_parameters(folder_path, sketch_func_id_check, kmer_size_from_file, is_containment_from_file, contain_compress, sketch_size, half_k, half_subk, drlevel);
			
			MinHashInvertedIndex local_index;
			MinHashInvertedIndex* idx_ptr = inverted_index;

			if(idx_ptr == nullptr) {
				// No pre-built index passed: try to load from sketch folder first
				if (!local_index.load_from_file(folder_path)) {
					cerr << "-----building MinHash inverted index for MST computation..." << endl;

					vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);

					#pragma omp parallel num_threads(threads)
					{
						int tid = omp_get_thread_num();
						auto& tl_index = thread_local_indices[tid];

						#pragma omp for schedule(dynamic, 100)
						for(size_t i = 0; i < sketches.size(); i++){
							if(sketches[i].minHash != nullptr) {
								sketches[i].id = i;
								vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
								for(uint64_t hash : hashes)
									tl_index[hash].push_back(i);
							}
						}
					}

					cerr << "-----merging thread-local indices..." << endl;
					for(int tid = 0; tid < threads; tid++){
						for(auto& entry : thread_local_indices[tid]){
							local_index.hash_map[entry.first].insert(
								local_index.hash_map[entry.first].end(),
								entry.second.begin(), entry.second.end()
							);
						}
					}
					// Save for future reuse (only if isSave — but index is useful anyway)
					if(isSave) local_index.save_to_file(folder_path);
				}
				idx_ptr = &local_index;
			}

			cerr << "-----MinHash inverted index: " << idx_ptr->hash_map.size() << " unique hashes" << endl;
			
			int kmer_size = kmer_size_from_file;
			bool is_containment = is_containment_from_file;
			if(kmer_size <= 0) {
				kmer_size = sketches[0].minHash->getKmerSize();
				is_containment = sketches[0].isContainment;
			}
			
			mst = compute_minhash_mst(sketches, 0, no_dense, is_containment, threads, denseArr, denseSpan, aniArr, threshold, kmer_size, idx_ptr);
			
			if(isSave){
				idx_ptr->save_to_file(folder_path);
			}
		} else {
			// Fallback to traditional modifyMST for non-MinHash or when MinHash is not available
			mst = modifyMST(sketches, 0, sketch_func_id, threads, no_dense, denseArr, denseSpan, aniArr);
		}
		double t3 = get_sec();
#ifdef Timer
		cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
#endif
		if(isSave){
			saveMST(sketches, mst, folder_path, sketchByFile);
		}
		double t4 = get_sec();
#ifdef Timer
		cerr << "========time of saveMST is: " << t4 - t3 << "========" << endl;
#endif

		//generate the Newick tree format
		if(is_newick_tree){
			string output_newick_file = outputFile + ".newick.tree";
			print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
			cerr << "-----write the newick tree into: " << output_newick_file << endl;
		}
		if(is_phylip_tree){
			string output_phylip_file = outputFile + ".phylip.tree";
			print_phylip_tree(sketches, mst, sketchByFile, output_phylip_file);
			cerr << "-----write the PHYLIP tree into: " << output_phylip_file << endl;
		}
		if(is_nexus_tree){
			string output_nexus_file = outputFile + ".nexus.tree";
			print_nexus_tree(sketches, mst, sketchByFile, output_nexus_file);
			cerr << "-----write the NEXUS tree into: " << output_nexus_file << endl;
		}
		if(is_linkage_matrix){
			string output_linkage_file = outputFile + ".linkage.txt";
			print_linkage_matrix(sketches, mst, output_linkage_file);
			cerr << "-----write the linkage matrix into: " << output_linkage_file << endl;
		}
		
		// Automatic threshold selection based on MST edge length distribution
		bool use_auto_threshold = false;  // Will be set from parameter
		if(use_auto_threshold){
			cerr << "-----analyzing MST edge length distribution for automatic threshold selection..." << endl;
			EdgeLengthStats stats = analyzeEdgeLengthDistribution(mst);
			// Use smaller min_gap_ratio (0.05 instead of 0.1) to detect more gaps
			// This helps find natural breakpoints even when distribution is dense
			vector<ThresholdCandidate> candidates = findThresholdCandidates(mst, 5, 0.05);
			ThresholdCandidate optimal = selectOptimalThreshold(candidates, mst);
			
			string threshold_analysis_file = outputFile + ".threshold_analysis.txt";
			printThresholdAnalysis(mst, stats, candidates, optimal, threshold_analysis_file);
			
			cerr << "-----optimal threshold: " << optimal.threshold << " (confidence: " << optimal.confidence 
			     << ", suggested level: " << optimal.level << ")" << endl;
			cerr << "-----threshold analysis written to: " << threshold_analysis_file << endl;
		}

		vector<EdgeInfo> forest = generateForest(mst, threshold);
		vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
		printResult(tmpClust, sketches, sketchByFile, outputFile, threshold);
		cerr << "-----write the cluster result into: " << outputFile << endl;
		cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

		//tune cluster by noise cluster
		if(!no_dense){
			if(isSave){
				saveANI(folder_path, aniArr, sketch_func_id);
				saveDense(folder_path, denseArr, denseSpan, sketches.size());
			}
			int alpha = 2;
			int denseIndex = threshold / 0.01;
			vector<int> totalNoiseArr;
			for(int i = 0; i < tmpClust.size(); i++){
				if(tmpClust[i].size() == 1) continue;
				vector<PairInt> curDenseArr;
				set<int> denseSet;
				for(int j = 0; j < tmpClust[i].size(); j++){
					int element = tmpClust[i][j];
					PairInt p(element, denseArr[denseIndex][element]);
					denseSet.insert(denseArr[denseIndex][element]);
					curDenseArr.push_back(p);
				}
				vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
				totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
			}
			cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
			forest = modifyForest(forest, totalNoiseArr, threads);
			cluster = generateClusterWithBfs(forest, sketches.size());
			string outputFileNew = outputFile + ".removeNoise";
			printResult(cluster, sketches, sketchByFile, outputFileNew);
			cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
			cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
			double t5 = get_sec();
#ifdef Timer
			cerr << "========time of tuning cluster is: " << t5 - t4 << "========" << endl;
#endif
		}

		//======clust-mst=======================================================================
#endif//endif GREEDY_CLUST
	}

// ==================== Leiden Clustering Functions ====================
#ifdef LEIDEN_CLUST

void clust_from_genome_leiden(const string inputFile, string outputFile, string folder_path, 
                              bool sketchByFile, const int kmerSize, const int drlevel, 
                              const int minLen, bool noSave, double threshold, 
                              double resolution, bool use_leiden, bool use_parallel_louvain, 
                              bool use_edge_parallel, bool use_warm_start, int knn_k, int threads){
	double time0 = get_sec();
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	bool isSave = !noSave;
	
	compute_kssd_sketches(sketches, info, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
	
	double time1 = get_sec();
#ifdef Timer
	cerr << "========time of generate sketches is: " << time1 - time0 << endl;
#endif
	
	cerr << "-----the size of sketches is: " << sketches.size() << endl;
	
	// Prepare graph save path (empty = don't save, saving 347M edges is very slow!)
	string graph_file = "";  // Changed from folder_path + "/leiden.graph" to disable slow file saving
	
	// Run graph-based clustering with inverted index graph construction
	vector<vector<int>> cluster;
	if (use_edge_parallel) {
		// Use edge-partitioning parallel Louvain (RECOMMENDED - simpler & more efficient)
		cluster = KssdEdgeParallelLouvainCluster(sketches, 0, threshold, threads, kmerSize, resolution, use_warm_start, knn_k, graph_file);
	} else if (use_parallel_louvain) {
		// Use graph-partitioning parallel Louvain (alternative approach)
		cluster = KssdParallelLouvainCluster(sketches, 0, threshold, threads, kmerSize, resolution, graph_file);
	} else {
		// Use standard Leiden or Louvain
		cluster = KssdLeidenCluster(sketches, 0, threshold, threads, kmerSize, resolution, use_leiden, knn_k, graph_file);
	}
	
	printKssdResult(cluster, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
	
	double time2 = get_sec();
#ifdef Timer
	cerr << "========time of Leiden clustering is: " << time2 - time1 << endl;
#endif
}

void clust_from_sketch_leiden(string folder_path, string outputFile, 
                              double threshold, double resolution, 
                              bool use_leiden, bool use_parallel_louvain, 
                              bool use_edge_parallel, bool use_warm_start, int knn_k, int threads){
	vector<KssdSketchInfo> sketches;
	bool sketchByFile;
	KssdParameters info;
	
	double time0 = get_sec();
	sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
	cerr << "-----the size of sketches is: " << sketches.size() << endl;
	double time1 = get_sec();
#ifdef Timer
	cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
	
	int kmer_size = info.half_k * 2;
	
	// Prepare graph save path (empty = don't save, saving millions of edges is very slow!)
	string graph_file = "";  // Changed from folder_path + "/leiden.graph" to disable slow file saving
	
	// Run graph-based clustering with inverted index graph construction
	vector<vector<int>> cluster;
	if (use_edge_parallel) {
		// Use edge-partitioning parallel Louvain (RECOMMENDED - simpler & more efficient)
		cluster = KssdEdgeParallelLouvainCluster(sketches, 0, threshold, threads, kmer_size, resolution, use_warm_start, knn_k, graph_file);
	} else if (use_parallel_louvain) {
		// Use graph-partitioning parallel Louvain (alternative approach)
		cluster = KssdParallelLouvainCluster(sketches, 0, threshold, threads, kmer_size, resolution, graph_file);
	} else {
		// Use standard Leiden or Louvain
		cluster = KssdLeidenCluster(sketches, 0, threshold, threads, kmer_size, resolution, use_leiden, knn_k, graph_file);
	}
	
	printKssdResult(cluster, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
	double time2 = get_sec();
#ifdef Timer
	cerr << "========time of Leiden clustering is: " << time2 - time1 << endl;
#endif
}

void clust_from_pregraph_leiden(string folder_path, string outputFile, 
                                double resolution, bool use_leiden, int threads){
	vector<KssdSketchInfo> sketches;
	bool sketchByFile;
	KssdParameters info;
	
	double time0 = get_sec();
	sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
	cerr << "-----the size of sketches is: " << sketches.size() << endl;
	double time1 = get_sec();
#ifdef Timer
	cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
	
	// Load and cluster from pre-built graph
	string graph_file = folder_path + "/leiden.graph";
	vector<vector<int>> cluster = KssdLeidenClusterFromGraph(graph_file, sketches.size(), resolution, use_leiden);
	
	printKssdResult(cluster, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
	double time2 = get_sec();
#ifdef Timer
	cerr << "========time of Leiden clustering from pre-built graph is: " << time2 - time1 << endl;
#endif
}

#endif // LEIDEN_CLUST

#ifdef DBSCAN_CLUST
// ==================== DBSCAN Clustering Functions ====================

void clust_from_sketch_dbscan(string folder_path, string outputFile, bool sketchByFile, double eps, int minPts, int threads, int knn_k, int max_posting) {
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	
	double time0 = get_sec();
	bool loaded_sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
	if(loaded_sketchByFile != sketchByFile) {
		cerr << "Warning: sketch format mismatch" << endl;
	}
	cerr << "-----the size of sketches is: " << sketches.size() << endl;
	double time1 = get_sec();
#ifdef Timer
	cerr << "========time of load sketches is: " << time1 - time0 << endl;
#endif
	
	int kmer_size = info.half_k * 2;
	
	// Run DBSCAN clustering
	DBSCANResult result = KssdDBSCAN(sketches, eps, minPts, kmer_size, threads, knn_k, max_posting);
	
	// Print results
	printKssdDBSCANResult(result, sketches, sketchByFile, outputFile, eps, minPts);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of " << outputFile << " is: " << result.num_clusters << endl;
	cerr << "-----the noise point number is: " << result.num_noise << endl;
	
	double time2 = get_sec();
#ifdef Timer
	cerr << "========time of DBSCAN clustering is: " << time2 - time1 << endl;
#endif
}

void clust_from_genome_dbscan(const string inputFile, string outputFile, string folder_path, bool sketchByFile, const int kmerSize, const int drlevel, const int minLen, bool noSave, double eps, int minPts, int threads, int knn_k, int max_posting) {
	bool isSave = !noSave;
	
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	string sketch_folder_path;
	
	double time0 = get_sec();
	compute_kssd_sketches(sketches, info, isSave, inputFile, sketch_folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
	cerr << "-----the size of sketches is: " << sketches.size() << endl;
	double time1 = get_sec();
#ifdef Timer
	cerr << "========time of compute sketches is: " << time1 - time0 << endl;
#endif
	
	// Run DBSCAN clustering
	DBSCANResult result = KssdDBSCAN(sketches, eps, minPts, kmerSize, threads, knn_k, max_posting);
	
	// Print results
	printKssdDBSCANResult(result, sketches, sketchByFile, outputFile, eps, minPts);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of " << outputFile << " is: " << result.num_clusters << endl;
	cerr << "-----the noise point number is: " << result.num_noise << endl;
	
	double time2 = get_sec();
#ifdef Timer
	cerr << "========time of DBSCAN clustering is: " << time2 - time1 << endl;
#endif
}

#endif // DBSCAN_CLUST
