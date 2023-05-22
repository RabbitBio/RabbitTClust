#include "sub_command.h"
#include <assert.h>

using namespace std;

#ifdef GREEDY_CLUST
void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
	int sketch_func_id_0; 
	vector<SketchInfo> pre_sketches; 
	bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0); 
	if(pre_sketch_by_file != sketch_by_file){
		cerr << "ERROR: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
		//cerr << "the output cluster file may not have the genome file name" << endl;
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
	compute_sketches(append_sketches, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, false, threads);
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
}
#endif

#ifndef GREEDY_CLUST
void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
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
	loadDense(pre_dense_arr, folder_path, pre_dense_span, pre_genome_number);

	int ** dense_arr;
	int dense_span = DENSE_SPAN;
	uint64_t* ani_arr;
	vector<EdgeInfo> append_mst = modifyMST(final_sketches, pre_sketch_size, sketch_func_id_0, threads, dense_arr, dense_span, ani_arr);
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

	vector<EdgeInfo> forest = generateForest(final_mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, final_sketches.size());
	printResult(tmpClust, final_sketches, pre_sketch_by_file, output_file);
	cerr << "-----write the cluster result into: " << output_file << endl;
	cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
	
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
		saveMST(final_sketches, final_mst, new_folder_path, sketch_by_file);
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

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, double threshold, int threads){
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

	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printResult(tmpClust, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;
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
#endif

void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool sketchByFile, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads){
	bool isSave = !noSave;
	vector<SketchInfo> sketches;
	int sketch_func_id;
	if(sketchFunc == "MinHash")	sketch_func_id = 0;
	else if(sketchFunc == "KSSD") sketch_func_id = 1;

	compute_sketches(sketches, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, isSave, threads);

	compute_clusters(sketches, sketchByFile, outputFile, is_newick_tree, folder_path, sketch_func_id, threshold, isSave, threads);
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

void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, double threshold, int threads){
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
	cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
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
	vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, denseArr, denseSpan, aniArr);
	double time2 = get_sec();
	#ifdef Timer
	cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
	#endif
	if(is_newick_tree){
		string output_newick_file = outputFile + ".newick.tree";
		print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
		cerr << "-----write the newick tree into: " << output_newick_file << endl;
	}
	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printResult(tmpClust, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

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
	double time3 = get_sec();
	#ifdef Timer
	cerr << "========time of generator forest and cluster is: " << time3 - time2 << "========" << endl;
	#endif
//=======clust-mst======================================================================
#endif
}

void compute_sketches(vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress,  bool isSave, int threads){
	double t0 = get_sec();
	if(sketchByFile){
		if(!sketchFiles(inputFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){
			cerr << "ERROR: generate_sketches(), cannot finish the sketch generation by genome files" << endl;
			exit(1);
		}
	}//end sketch by sequence
	else{
		if(!sketchSequences(inputFile, kmerSize, sketchSize, minLen, sketchFunc, isContainment, containCompress, sketches, threads)){
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

void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads){
	vector<vector<int>> cluster;
	double t2 = get_sec();
#ifdef GREEDY_CLUST
//======clust-greedy====================================================================
	cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
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
	vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, denseArr, denseSpan, aniArr);
	double t3 = get_sec();
	#ifdef Timer
	cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
	#endif
	if(isSave){
		saveANI(folder_path, aniArr, sketch_func_id);
		saveDense(folder_path, denseArr, denseSpan, sketches.size());
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


//	for(int i = 0; i < denseSpan; i++){
//		for(int j = 0; j < sketches.size(); j++){
//			cout << denseArr[i][j] << endl;
//		}
//	}
	
	vector<EdgeInfo> forest = generateForest(mst, threshold);
	vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
	printResult(tmpClust, sketches, sketchByFile, outputFile);
	cerr << "-----write the cluster result into: " << outputFile << endl;
	cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;
	//tune cluster by noise cluster
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
//======clust-mst=======================================================================
#endif//endif GREEDY_CLUST
}



