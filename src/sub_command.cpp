#include "sub_command.h"
#include <assert.h>
#include "cluster_postprocess.h"
#include <fstream>
#include <sstream>
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
#endif

#ifndef GREEDY_CLUST
void append_clust_mst_fast(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, bool isContainment, int min_len, bool no_save, double threshold, int threads){
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
	transSketches(final_sketches, append_info, new_folder_path, threads);

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
	vector<EdgeInfo> append_mst = compute_kssd_mst(final_sketches, append_info, new_folder_path, pre_sketch_size, no_dense, isContainment, threads, dense_arr, dense_span, ani_arr, threshold);
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

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool is_linkage_matrix, bool no_dense, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
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
	
	// Use inverted index optimization for MinHash MST append
	vector<EdgeInfo> append_mst;
	if(sketch_func_id_0 == 0 && final_sketches.size() > 0 && final_sketches[0].minHash != nullptr) {
		// Build inverted index from all final_sketches (pre + append) using thread-local indices
		MinHashInvertedIndex inverted_index;
		cerr << "-----building MinHash inverted index for append MST computation..." << endl;
		
		// Use thread-local indices to reduce lock contention
		vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);
		
		#pragma omp parallel num_threads(threads)
		{
			int tid = omp_get_thread_num();
			auto& local_index = thread_local_indices[tid];
			
			#pragma omp for schedule(dynamic, 100)
			for(size_t i = 0; i < final_sketches.size(); i++){
				if(final_sketches[i].minHash != nullptr) {
					final_sketches[i].id = i;  // Ensure ID matches position
					vector<uint64_t> hashes = final_sketches[i].minHash->storeMinHashes();
					for(uint64_t hash : hashes) {
						local_index[hash].push_back(i);
					}
				}
			}
		}
		
		// Merge thread-local indices (no lock needed here)
		cerr << "-----merging thread-local indices..." << endl;
		for(int tid = 0; tid < threads; tid++){
			for(auto& entry : thread_local_indices[tid]){
				uint64_t hash = entry.first;
				auto& seq_ids = entry.second;
				inverted_index.hash_map[hash].insert(
					inverted_index.hash_map[hash].end(), 
					seq_ids.begin(), seq_ids.end()
				);
			}
		}
		
		cerr << "-----MinHash inverted index built: " << inverted_index.hash_map.size() << " unique hashes" << endl;
		
		// Use compute_minhash_mst with start_index=pre_sketch_size (only compute edges for new sketches)
		append_mst = compute_minhash_mst(final_sketches, pre_sketch_size, no_dense, is_containment, threads, dense_arr, dense_span, ani_arr, threshold, kmer_size, &inverted_index);
	} else {
		// Fallback to traditional modifyMST for non-MinHash
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
void clust_from_mst_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads){
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

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads){
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

void clust_from_genome_fast(const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index){
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
	compute_kssd_clusters(sketches, info, sketchByFile, no_dense, isContainment, folder_path, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, threshold, isSave, threads, dedup_dist, reps_per_cluster, save_rep_index, &inverted_index);

}

void compute_kssd_clusters(vector<KssdSketchInfo>& sketches, const KssdParameters info, bool sketchByFile, bool no_dense, bool isContainment, const string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, double threshold, bool isSave, int threads, double dedup_dist, int reps_per_cluster, bool save_rep_index, KssdInvertedIndex* inverted_index){
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

	void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool sketchByFile, bool no_dense, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads, bool use_inverted_index, bool save_rep_index){
		bool isSave = !noSave;
		vector<SketchInfo> sketches;
		int sketch_func_id;
		if(sketchFunc == "MinHash")	sketch_func_id = 0;
		else if(sketchFunc == "KSSD") sketch_func_id = 1;

		compute_sketches(sketches, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, isSave, threads, nullptr);

		compute_clusters(sketches, sketchByFile, outputFile, is_newick_tree, is_linkage_matrix, is_auto_threshold, is_stability, no_dense, folder_path, sketch_func_id, threshold, isSave, threads, use_inverted_index, save_rep_index);
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

	void clust_from_sketch_fast(string folder_path, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool no_dense, bool isContainment, double threshold, int threads, double dedup_dist, int reps_per_cluster, bool use_inverted_index, bool save_rep_index){
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

	void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, bool is_auto_threshold, bool is_stability, bool no_dense, double threshold, int threads, bool use_inverted_index, bool save_rep_index){
		vector<SketchInfo> sketches;
		vector<vector<int>> cluster;
		int sketch_func_id;
		bool sketchByFile;
#ifdef GREEDY_CLUST
		//======clust-greedy====================================================================
		double time0 = get_sec();
		sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
		cerr << "-----the size of sketches is: " << sketches.size() << endl;
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
			
			// Build inverted index from loaded sketches using thread-local indices (like KSSD)
			MinHashInvertedIndex inverted_index;
			cerr << "-----building MinHash inverted index for MST computation..." << endl;
			
			// Use thread-local indices to reduce lock contention (similar to KSSD optimization)
			vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);
			
			#pragma omp parallel num_threads(threads)
			{
				int tid = omp_get_thread_num();
				auto& local_index = thread_local_indices[tid];
				
				#pragma omp for schedule(dynamic, 100)
				for(size_t i = 0; i < sketches.size(); i++){
					if(sketches[i].minHash != nullptr) {
						sketches[i].id = i;  // Ensure ID matches position
						vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
						for(uint64_t hash : hashes) {
							local_index[hash].push_back(i);
						}
					}
				}
			}
			
			// Merge thread-local indices (no lock needed here)
			cerr << "-----merging thread-local indices..." << endl;
			for(int tid = 0; tid < threads; tid++){
				for(auto& entry : thread_local_indices[tid]){
					uint64_t hash = entry.first;
					auto& seq_ids = entry.second;
					inverted_index.hash_map[hash].insert(
						inverted_index.hash_map[hash].end(), 
						seq_ids.begin(), seq_ids.end()
					);
				}
			}
			
			cerr << "-----MinHash inverted index built: " << inverted_index.hash_map.size() << " unique hashes" << endl;
			
			// Use kmerSize from file (default 21 for MinHash), fallback to sketch object if file read fails
			int kmer_size = kmer_size_from_file;
			bool is_containment = is_containment_from_file;
			if(kmer_size <= 0) {
				// Fallback: get from sketch object
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
			double t2 = get_sec();
#ifdef Timer
			cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
		}
	}

	void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, bool is_linkage_matrix, bool is_auto_threshold, bool is_stability, bool no_dense, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads, bool use_inverted_index, bool save_rep_index){
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
			// Read kmerSize from sketch parameters file (more reliable than from sketch object)
			int kmer_size_from_file = 21;  // Default for MinHash
			int sketch_func_id_check, contain_compress, sketch_size, half_k, half_subk, drlevel;
			bool is_containment_from_file;
			read_sketch_parameters(folder_path, sketch_func_id_check, kmer_size_from_file, is_containment_from_file, contain_compress, sketch_size, half_k, half_subk, drlevel);
			
			// Build inverted index from sketches using thread-local indices (like KSSD)
			MinHashInvertedIndex inverted_index;
			cerr << "-----building MinHash inverted index for MST computation..." << endl;
			
			// Use thread-local indices to reduce lock contention
			vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);
			
			#pragma omp parallel num_threads(threads)
			{
				int tid = omp_get_thread_num();
				auto& local_index = thread_local_indices[tid];
				
				#pragma omp for schedule(dynamic, 100)
				for(size_t i = 0; i < sketches.size(); i++){
					if(sketches[i].minHash != nullptr) {
						sketches[i].id = i;  // Ensure ID matches position
						vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
						for(uint64_t hash : hashes) {
							local_index[hash].push_back(i);
						}
					}
				}
			}
			
			// Merge thread-local indices (no lock needed here)
			cerr << "-----merging thread-local indices..." << endl;
			for(int tid = 0; tid < threads; tid++){
				for(auto& entry : thread_local_indices[tid]){
					uint64_t hash = entry.first;
					auto& seq_ids = entry.second;
					inverted_index.hash_map[hash].insert(
						inverted_index.hash_map[hash].end(), 
						seq_ids.begin(), seq_ids.end()
					);
				}
			}
			
			cerr << "-----MinHash inverted index built: " << inverted_index.hash_map.size() << " unique hashes" << endl;
			
			// Use kmerSize from file (default 21 for MinHash), fallback to sketch object if file read fails
			int kmer_size = kmer_size_from_file;
			bool is_containment = is_containment_from_file;
			if(kmer_size <= 0) {
				// Fallback: get from sketch object
				kmer_size = sketches[0].minHash->getKmerSize();
				is_containment = sketches[0].isContainment;
			}
			
			// Use compute_minhash_mst with inverted index
			mst = compute_minhash_mst(sketches, 0, no_dense, is_containment, threads, denseArr, denseSpan, aniArr, threshold, kmer_size, &inverted_index);
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

#ifdef USE_MPI
#include <mpi.h>
#include "common.hpp"

static void SafeBcast(char* buf, size_t total_size, int root, MPI_Comm comm) {
	const size_t chunk = 512 * 1024 * 1024;
	for (size_t offset = 0; offset < total_size; offset += chunk) {
		size_t send = (total_size - offset < chunk) ? (total_size - offset) : chunk;
		MPI_Bcast(buf + offset, (int)send, MPI_CHAR, root, comm);
	}
}

static void distribute_compute_clusters(int my_rank, int comm_sz, vector<KssdSketchInfo>& sketches, const KssdParameters& info, bool sketchByFile, string output_file, bool is_newick_tree, string folder_path, double threshold, bool isSave, int threads, bool no_dense, bool isContainment, char* index_buffer, size_t index_size, char* dict_buffer, size_t dict_size) {
	if (my_rank == 0) cerr << "-----MPI: computing local MST (OpenMP within each rank)..." << endl;
	int** denseArr = nullptr;
	uint64_t* aniArr = nullptr;
	int denseSpan = DENSE_SPAN;
	double t_mst_start = get_sec();
	vector<EdgeInfo> my_mst = compute_kssd_mst_mpi(my_rank, comm_sz, sketches, info, folder_path, no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold, index_buffer, index_size, dict_buffer, dict_size);
	double t_mst_end = get_sec();
	if (my_rank == 0) cerr << "========time of MPI local MST computation is: " << t_mst_end - t_mst_start << " s========" << endl;
	if (my_rank == 0) cerr << "-----MPI: gathering edges to rank 0, running Kruskal..." << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank != 0) {
		size_t edge_count = my_mst.size();
		MPI_Send(&edge_count, 1, MPI_UINT64_T, 0, my_rank, MPI_COMM_WORLD);
		if (edge_count > 0)
			MPI_Send(my_mst.data(), (int)(edge_count * sizeof(EdgeInfo)), MPI_CHAR, 0, my_rank + comm_sz, MPI_COMM_WORLD);
		if (denseArr) {
			for (int i = 0; i < denseSpan; i++) delete[] denseArr[i];
			delete[] denseArr;
		}
		if (aniArr) delete[] aniArr;
		return;
	}

	double t_gather_start = get_sec();
	vector<EdgeInfo> sum_mst(my_mst.begin(), my_mst.end());
	for (int id = 1; id < comm_sz; id++) {
		size_t recv_edge_count;
		MPI_Recv(&recv_edge_count, 1, MPI_UINT64_T, id, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (recv_edge_count > 0) {
			vector<EdgeInfo> cur_MST(recv_edge_count);
			MPI_Recv(cur_MST.data(), (int)(recv_edge_count * sizeof(EdgeInfo)), MPI_CHAR, id, id + comm_sz, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum_mst.insert(sum_mst.end(), cur_MST.begin(), cur_MST.end());
		}
	}
	cerr << "========total edges gathered: " << sum_mst.size() << "========" << endl;
	sort(sum_mst.begin(), sum_mst.end(), cmpEdge);
	vector<EdgeInfo> final_mst = kruskalAlgorithm(sum_mst, (int)sketches.size());
	double t_gather_end = get_sec();
	cerr << "========time of gather + Kruskal is: " << t_gather_end - t_gather_start << " s========" << endl;
	if (denseArr) {
		for (int i = 0; i < denseSpan; i++) delete[] denseArr[i];
		delete[] denseArr;
	}
	if (aniArr) delete[] aniArr;

	if (is_newick_tree) {
		string output_newick_file = output_file + ".newick.tree";
		print_kssd_newick_tree(sketches, final_mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	vector<EdgeInfo> forest = generateForest(final_mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, (int)sketches.size());
	printKssdResult(tmpClust, sketches, sketchByFile, output_file);
	cerr << "-----write the cluster result into: " << output_file << endl;
	cerr << "-----the cluster number of " << output_file << " is: " << tmpClust.size() << endl;
}

void compute_kssd_sketches_mpi(int my_rank, int comm_sz, vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads) {
	double t_sketch_start = get_sec();
	if (sketchByFile) {
		vector<string> file_list;
		ifstream ifs(inputFile);
		if (!ifs.is_open()) {
			cerr << "ERROR: cannot open input file list: " << inputFile << endl;
			exit(1);
		}
		string line;
		while (getline(ifs, line)) { if (!line.empty()) file_list.push_back(line); }
		ifs.close();
		if (!sketchFileWithKssd_mpi(file_list, my_rank, comm_sz, minLen, kmerSize, drlevel, sketches, info, threads)) {
			cerr << "ERROR: compute_kssd_sketches_mpi: sketch generation failed" << endl;
			exit(1);
		}
	} else {
		cerr << "ERROR: MPI path currently supports only -l/--list (sketch by file)" << endl;
		exit(1);
	}
	double t_sketch_end = get_sec();
	if (my_rank == 0) cerr << "========time of MPI sketch generation is: " << t_sketch_end - t_sketch_start << " s========" << endl;
	cerr << "-----rank " << my_rank << " sketches size: " << sketches.size() << " (sketch phase used OpenMP)" << endl;

	MPI_Barrier(MPI_COMM_WORLD);

	if (isSave) {
		if (my_rank == 0) {
			folder_path = currentDataTime();
			string command = "mkdir -p " + folder_path;
			if (system(command.c_str()) != 0) {
				cerr << "FATAL: Rank 0 failed to create directory: " << folder_path << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		int path_len;
		if (my_rank == 0) path_len = (int)folder_path.length();
		MPI_Bcast(&path_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

		char folder_path_buffer[path_len + 1];
		if (my_rank == 0) strcpy(folder_path_buffer, folder_path.c_str());
		MPI_Bcast(folder_path_buffer, path_len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
		if (my_rank != 0) folder_path = string(folder_path_buffer);

		if (my_rank == 0) cerr << "-----All sketches will be saved to directory: " << folder_path << endl;
		MPI_Barrier(MPI_COMM_WORLD);

		string info_file = folder_path + "/kssd.info.sketch";

		size_t local_sketch_count = sketches.size();
		size_t total_sketch_count = 0;
		MPI_Reduce(&local_sketch_count, &total_sketch_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_rank == 0) {
			FILE* fp_info = fopen(info_file.c_str(), "wb");
			if (fp_info) {
				fwrite(&sketchByFile, sizeof(bool), 1, fp_info);
				fwrite(&total_sketch_count, sizeof(size_t), 1, fp_info);
				fclose(fp_info);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int rank_to_write = 0; rank_to_write < comm_sz; ++rank_to_write) {
			if (my_rank == rank_to_write) {
				FILE* fp_info = fopen(info_file.c_str(), "ab");
				if (fp_info) {
					append_binary_genome_info(fp_info, sketches, sketchByFile);
					fclose(fp_info);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		string hash_file = folder_path + "/kssd.hash.sketch";

		if (my_rank == 0) {
			FILE* fp_hash = fopen(hash_file.c_str(), "wb");
			if (fp_hash) {
				fwrite(&info, sizeof(KssdParameters), 1, fp_hash);
				fclose(fp_hash);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int rank_to_write = 0; rank_to_write < comm_sz; ++rank_to_write) {
			if (my_rank == rank_to_write) {
				FILE* fp_hash = fopen(hash_file.c_str(), "ab");
				if (fp_hash) {
					append_binary_hash_data(fp_hash, sketches);
					fclose(fp_hash);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		if (my_rank == 0) cerr << "-----Successfully saved all kssd sketches into: " << folder_path << endl;
	}
}

void clust_from_sketches_fast_MPI(int my_rank, int comm_sz, int /* half_k */, int drlevel, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const double threshold, bool noSave, int threads) {
	bool isSave = !noSave;
	vector<KssdSketchInfo> sketches;
	KssdParameters info;

	double t_load_start = get_sec();
	sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
	std::sort(sketches.begin(), sketches.end(), [](const KssdSketchInfo& a, const KssdSketchInfo& b) { return a.hash32_arr.size() < b.hash32_arr.size(); });
	size_t total = sketches.size();
	const size_t chunkSize = 10000;
	size_t numChunks = (total + chunkSize - 1) / chunkSize;
	vector<vector<KssdSketchInfo>> buckets(numChunks);
	for (size_t i = 0; i < total; i++) buckets[i % numChunks].push_back(sketches[i]);
	sketches.clear();
	for (size_t i = 0; i < numChunks; i++)
		sketches.insert(sketches.end(), buckets[i].begin(), buckets[i].end());
	double t_load_end = get_sec();
	if (my_rank == 0) cerr << "========time of loading & sorting sketches is: " << t_load_end - t_load_start << " s========" << endl;

	double t_index_start = get_sec();
	char* sum_info_buffer = nullptr, * sum_hash_buffer = nullptr, * sum_index_buffer = nullptr, * sum_dict_buffer = nullptr;
	size_t sum_info_size = 0, sum_hash_size = 0, sum_index_size = 0, sum_dict_size = 0;
	transSketches_in_memory(sketches, info, threads, sum_info_buffer, sum_info_size, sum_hash_buffer, sum_hash_size, sum_index_buffer, sum_index_size, sum_dict_buffer, sum_dict_size);
	double t_index_end = get_sec();
	if (my_rank == 0) cerr << "========time of building inverted index is: " << t_index_end - t_index_start << " s========" << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	distribute_compute_clusters(my_rank, comm_sz, sketches, info, sketchByFile, outputFile, is_newick_tree, folder_path, threshold, isSave, threads, no_dense, isContainment, sum_index_buffer, sum_index_size, sum_dict_buffer, sum_dict_size);
	delete[] sum_info_buffer;
	delete[] sum_hash_buffer;
	delete[] sum_index_buffer;
	delete[] sum_dict_buffer;
}

void clust_from_genomes_fast_MPI(int my_rank, int comm_sz, const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads) {
	bool isSave = !noSave;
	vector<KssdSketchInfo> sketches;
	KssdParameters info;
	compute_kssd_sketches_mpi(my_rank, comm_sz, sketches, info, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
	int half_k = info.half_k;
	clust_from_sketches_fast_MPI(my_rank, comm_sz, half_k, drlevel, outputFile, folder_path, is_newick_tree, no_dense, sketchByFile, isContainment, threshold, noSave, threads);
}

// ============================================================================
// MinHash MPI functions
// ============================================================================

static void distribute_compute_minhash_clusters(int my_rank, int comm_sz, vector<SketchInfo>& sketches, bool sketchByFile, string output_file, bool is_newick_tree, double threshold, bool isSave, int threads, bool no_dense, bool isContainment, int kmerSize, MinHashInvertedIndex* inverted_index) {
	if (my_rank == 0) cerr << "-----MPI: computing local MinHash MST (OpenMP within each rank)..." << endl;
	int** denseArr = nullptr;
	uint64_t* aniArr = nullptr;
	int denseSpan = DENSE_SPAN;
	double t_mst_start = get_sec();
	vector<EdgeInfo> my_mst = compute_minhash_mst_mpi(my_rank, comm_sz, sketches, no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold, kmerSize, inverted_index);
	double t_mst_end = get_sec();
	if (my_rank == 0) cerr << "========time of MPI local MinHash MST computation is: " << t_mst_end - t_mst_start << " s========" << endl;
	if (my_rank == 0) cerr << "-----MPI: gathering edges to rank 0, running Kruskal..." << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank != 0) {
		size_t edge_count = my_mst.size();
		MPI_Send(&edge_count, 1, MPI_UINT64_T, 0, my_rank, MPI_COMM_WORLD);
		if (edge_count > 0)
			MPI_Send(my_mst.data(), (int)(edge_count * sizeof(EdgeInfo)), MPI_CHAR, 0, my_rank + comm_sz, MPI_COMM_WORLD);
		if (denseArr) {
			for (int i = 0; i < denseSpan; i++) delete[] denseArr[i];
			delete[] denseArr;
		}
		if (aniArr) delete[] aniArr;
		return;
	}

	double t_gather_start = get_sec();
	vector<EdgeInfo> sum_mst(my_mst.begin(), my_mst.end());
	for (int id = 1; id < comm_sz; id++) {
		size_t recv_edge_count;
		MPI_Recv(&recv_edge_count, 1, MPI_UINT64_T, id, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (recv_edge_count > 0) {
			vector<EdgeInfo> cur_MST(recv_edge_count);
			MPI_Recv(cur_MST.data(), (int)(recv_edge_count * sizeof(EdgeInfo)), MPI_CHAR, id, id + comm_sz, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum_mst.insert(sum_mst.end(), cur_MST.begin(), cur_MST.end());
		}
	}
	cerr << "========total edges gathered: " << sum_mst.size() << "========" << endl;
	sort(sum_mst.begin(), sum_mst.end(), cmpEdge);
	vector<EdgeInfo> final_mst = kruskalAlgorithm(sum_mst, (int)sketches.size());
	double t_gather_end = get_sec();
	cerr << "========time of gather + Kruskal is: " << t_gather_end - t_gather_start << " s========" << endl;
	if (denseArr) {
		for (int i = 0; i < denseSpan; i++) delete[] denseArr[i];
		delete[] denseArr;
	}
	if (aniArr) delete[] aniArr;

	if (is_newick_tree) {
		string output_newick_file = output_file + ".newick.tree";
		print_newick_tree(sketches, final_mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	vector<EdgeInfo> forest = generateForest(final_mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, (int)sketches.size());
	printResult(tmpClust, sketches, sketchByFile, output_file);
	cerr << "-----write the cluster result into: " << output_file << endl;
	cerr << "-----the cluster number of " << output_file << " is: " << tmpClust.size() << endl;
}

void compute_minhash_sketches_mpi(int my_rank, int comm_sz, vector<SketchInfo>& sketches, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int sketchSize, string sketchFunc, bool isContainment, int containCompress, int threads) {
	double t_sketch_start = get_sec();
	if (sketchByFile) {
		vector<string> file_list;
		ifstream ifs(inputFile);
		if (!ifs.is_open()) {
			cerr << "ERROR: cannot open input file list: " << inputFile << endl;
			exit(1);
		}
		string line;
		while (getline(ifs, line)) { if (!line.empty()) file_list.push_back(line); }
		ifs.close();
		if (!sketchFiles_mpi(file_list, my_rank, comm_sz, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)) {
			cerr << "ERROR: compute_minhash_sketches_mpi: sketch generation failed" << endl;
			exit(1);
		}
	} else {
		cerr << "ERROR: MPI path currently supports only -l/--list (sketch by file)" << endl;
		exit(1);
	}
	double t_sketch_end = get_sec();
	if (my_rank == 0) cerr << "========time of MPI MinHash sketch generation is: " << t_sketch_end - t_sketch_start << " s========" << endl;
	cerr << "-----rank " << my_rank << " MinHash sketches size: " << sketches.size() << endl;

	MPI_Barrier(MPI_COMM_WORLD);

	if (isSave) {
		if (my_rank == 0) {
			folder_path = currentDataTime();
			string command = "mkdir -p " + folder_path;
			if (system(command.c_str()) != 0) {
				cerr << "FATAL: Rank 0 failed to create directory: " << folder_path << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		int path_len;
		if (my_rank == 0) path_len = (int)folder_path.length();
		MPI_Bcast(&path_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

		char folder_path_buffer[path_len + 1];
		if (my_rank == 0) strcpy(folder_path_buffer, folder_path.c_str());
		MPI_Bcast(folder_path_buffer, path_len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
		if (my_rank != 0) folder_path = string(folder_path_buffer);

		if (my_rank == 0) cerr << "-----All MinHash sketches will be saved to directory: " << folder_path << endl;
		MPI_Barrier(MPI_COMM_WORLD);

		string info_file = folder_path + "/info.sketch";

		size_t local_sketch_count = sketches.size();
		size_t total_sketch_count = 0;
		MPI_Reduce(&local_sketch_count, &total_sketch_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_rank == 0) {
			FILE* fp_info = fopen(info_file.c_str(), "wb");
			if (fp_info) {
				fwrite(&sketchByFile, sizeof(bool), 1, fp_info);
				fwrite(&total_sketch_count, sizeof(size_t), 1, fp_info);
				fclose(fp_info);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int rank_to_write = 0; rank_to_write < comm_sz; ++rank_to_write) {
			if (my_rank == rank_to_write) {
				FILE* fp_info = fopen(info_file.c_str(), "ab");
				if (fp_info) {
					append_minhash_genome_info(fp_info, sketches, sketchByFile);
					fclose(fp_info);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		string hash_file = folder_path + "/hash.sketch";

		if (my_rank == 0) {
			FILE* fp_hash = fopen(hash_file.c_str(), "wb");
			if (fp_hash) {
				int sketch_func_id = 0;
				fwrite(&sketch_func_id, sizeof(int), 1, fp_hash);
				fwrite(&kmerSize, sizeof(int), 1, fp_hash);
				fwrite(&isContainment, sizeof(bool), 1, fp_hash);
				if (isContainment)
					fwrite(&containCompress, sizeof(int), 1, fp_hash);
				else
					fwrite(&sketchSize, sizeof(int), 1, fp_hash);
				fclose(fp_hash);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int rank_to_write = 0; rank_to_write < comm_sz; ++rank_to_write) {
			if (my_rank == rank_to_write) {
				FILE* fp_hash = fopen(hash_file.c_str(), "ab");
				if (fp_hash) {
					append_minhash_hash_data(fp_hash, sketches);
					fclose(fp_hash);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		if (my_rank == 0) cerr << "-----Successfully saved all MinHash sketches into: " << folder_path << endl;
	}
}

void clust_from_sketches_MPI(int my_rank, int comm_sz, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool isContainment, const double threshold, bool noSave, int threads) {
	bool isSave = !noSave;
	vector<SketchInfo> sketches;
	int sketch_func_id;

	double t_load_start = get_sec();
	bool sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
	if (sketch_func_id != 0) {
		cerr << "ERROR: clust_from_sketches_MPI: expected MinHash sketches (sketch_func_id=0), got " << sketch_func_id << endl;
		cerr << "For KSSD sketches, use --fast --mpi instead" << endl;
		exit(1);
	}
	sort(sketches.begin(), sketches.end(), cmpGenomeSize);
	double t_load_end = get_sec();
	if (my_rank == 0) cerr << "========time of loading & sorting MinHash sketches is: " << t_load_end - t_load_start << " s========" << endl;

	double t_index_start = get_sec();
	MinHashInvertedIndex inverted_index;
	vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> thread_local_indices(threads);
	#pragma omp parallel num_threads(threads)
	{
		int tid = omp_get_thread_num();
		auto& local_index = thread_local_indices[tid];
		#pragma omp for schedule(dynamic, 100)
		for (size_t i = 0; i < sketches.size(); i++) {
			if (sketches[i].minHash != nullptr) {
				sketches[i].id = i;
				vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
				for (uint64_t hash : hashes)
					local_index[hash].push_back(i);
			}
		}
	}
	for (int tid = 0; tid < threads; tid++) {
		for (auto& entry : thread_local_indices[tid]) {
			inverted_index.hash_map[entry.first].insert(
				inverted_index.hash_map[entry.first].end(),
				entry.second.begin(), entry.second.end());
		}
	}
	double t_index_end = get_sec();
	if (my_rank == 0) cerr << "========time of building MinHash inverted index is: " << t_index_end - t_index_start << " s========" << endl;
	if (my_rank == 0) cerr << "-----MinHash inverted index: " << inverted_index.hash_map.size() << " unique hashes" << endl;

	int kmerSize = 21;
	if (!sketches.empty() && sketches[0].minHash != nullptr) {
		kmerSize = sketches[0].minHash->getKmerSize();
		isContainment = sketches[0].isContainment;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	distribute_compute_minhash_clusters(my_rank, comm_sz, sketches, sketchByFile, outputFile, is_newick_tree, threshold, isSave, threads, no_dense, isContainment, kmerSize, &inverted_index);
}

void clust_from_genomes_MPI(int my_rank, int comm_sz, const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const int sketchSize, const double threshold, string sketchFunc, int containCompress, const int minLen, bool noSave, int threads) {
	bool isSave = !noSave;
	vector<SketchInfo> sketches;
	compute_minhash_sketches_mpi(my_rank, comm_sz, sketches, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, threads);
	clust_from_sketches_MPI(my_rank, comm_sz, outputFile, folder_path, is_newick_tree, no_dense, isContainment, threshold, noSave, threads);
}
#endif
