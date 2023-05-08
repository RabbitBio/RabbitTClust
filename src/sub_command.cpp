#include "sub_command.h"

using namespace std;

void clust_from_mst(string folder_path, string outputFile, double threshold, int threads){
	vector<SketchInfo> sketches;
	vector<EdgeInfo> mst;
	vector<vector<int>> cluster;
	bool sketchByFile = loadMST(folder_path, sketches, mst);
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

void clust_from_genomes(string inputFile, string outputFile, bool sketchByFile, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads){
	bool isSave = !noSave;
	vector<SketchInfo> sketches;
	int sketch_func_id;
	if(sketchFunc == "MinHash")	sketch_func_id = 0;
	else if(sketchFunc == "KSSD") sketch_func_id = 1;

	compute_sketches(sketches, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, isSave, threads);

	compute_clusters(sketches, sketchByFile, outputFile, folder_path, sketch_func_id, threshold, isSave, threads);

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

void clust_from_sketches(string folder_path, string outputFile, double threshold, int threads){
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
	vector<EdgeInfo> mst = modifyMST(sketches, sketch_func_id, threads, denseArr, denseSpan, aniArr, outputFile, threshold);
	double time2 = get_sec();
	#ifdef Timer
	cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
	#endif
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

void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads){
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
	vector<EdgeInfo> mst = modifyMST(sketches, sketch_func_id, threads, denseArr, denseSpan, aniArr, outputFile, threshold);
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



