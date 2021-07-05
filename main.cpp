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
 * section 3: compute the distance matrix and create the Complete Distance Graph(CDG).
 * section 4: generate the Minimum Spanning Tree(MST) with the CDG by kruskal algrithm.
 * section 5: generate the clusters with the MST using different distance threshold.
 *
 * Author: Xiaoming Xu
 * Mar 5, 2021
 *
 */

#include <iostream>
#include "SketchInfo.h"
#include "Sketch.h"// need to add the include path in Makefile.
#include <sys/time.h>
#include <zlib.h>
//#include "kseq.h"
#include <stdio.h>
#include "MST.h"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include "UnionFind.h"
#include <algorithm>
#include "parameter.h"

//#include "getSketch.h"

//KSEQ_INIT(gzFile, gzread);

using namespace std;

bool cmpSketch(SketchInfo s1, SketchInfo s2){
	return s1.index < s2.index;
}

double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

void printUsage(void){
	fprintf(stdout, "usage: clust [-h] [-l] [-t] <int> [-d] <double> -F <string> [-f] genomeFile\n");
	fprintf(stdout, "usage: clust [-h] [-f] [-d] <double> genomeInfo mstInfo\n");
	fprintf(stdout, "	-h		: this help message\n");
	fprintf(stdout, "	-l		: genome clustering, inputFile is the path list of the genome files\n");
	fprintf(stdout, "	-t		: genome clustering, set the thread number\n");
	fprintf(stdout, "	-d		: set the threshold of the clusters from the Minimum Spanning Tree\n");
	fprintf(stdout, "	-f		: input files are genomeInfo and MST content\n");
	fprintf(stdout, "	-F		: sketch function, includes MinHash, WMH, OMH, HLL\n"); 
	
}

int main(int argc, char * argv[]){

	//section 1: init parameters
	int argIndex = 1;
	string inputFile = "genome.fna";
	string inputFile1 = "genome.info";
	string sketchFunc = "MinHash";
	int threads = 1;
	bool sketchByFile = false;
	bool useMST = false;
	double threshold = 1.0;
	while(argIndex < argc){
		if(useMST && argIndex + 2 == argc){
			inputFile1 = argv[argIndex];
		}
		else if(argIndex + 1 == argc){
			inputFile = argv[argIndex];
		}
		else{
			switch(argv[argIndex][1]){
				case 'h':
					printUsage();
					break;
				case 't':
					threads = atoi(argv[++argIndex]);
					if(threads < 1 || threads > 128){
						fprintf(stderr, "Invalid thread number %d\n", threads);
						return 1;
					}
					break;
				case 'l':
					sketchByFile = true;
					fprintf(stderr, "sketch by file: \n");
					break;
				case 'd':
					threshold = stod(argv[++argIndex]);
					fprintf(stderr, "set the threshold is: %lf \n", threshold);
					break;
				case 'f':
					useMST = true;
					fprintf(stderr, "input as Munimum Spanning Tree \n");
					break;
				case 'F':
					sketchFunc = argv[++argIndex];
					fprintf(stderr, "the sketch function is: %s \n", sketchFunc.c_str());
					break;
				default:
					fprintf(stderr, "Invalid option %s\n", argv[argIndex]);
					printUsage();
					return 1;
			}
		}
		++argIndex;
	}//end while file;

	if(useMST){
		cout << "input files is MST information and genome informations" << endl;
		fstream fs(inputFile);
		fstream fs1(inputFile1);

		if(!fs){
			fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
			printUsage();
			return 1;
		}
		if(!fs1){
			fprintf(stderr, "error open the inputFile: %s\n", inputFile1.c_str());
			printUsage();
			return 1;
		}
		vector<EdgeInfo> mst;
		string line;
		while(getline(fs, line)){
			stringstream ss;
			ss << line;
			int preNode, sufNode;
			double distance;
			ss >> preNode >> sufNode >> distance;

			EdgeInfo tmpEdge;
			tmpEdge.preNode = preNode;
			tmpEdge.sufNode = sufNode;
			tmpEdge.dist = distance;
			//printf("<%d, %d, %lf>\n", tmpEdge.preNode, tmpEdge.sufNode, tmpEdge.dist); 
			mst.push_back(tmpEdge);
		}//end while

		vector<SimilarityInfo> similarityInfos;
		while(1){
			if(!getline(fs1, line)) break;
			stringstream ss;
			ss << line;
			int id, length;
			string name, comment;
			ss >> id >> name >> length;
			
			getline(fs1, comment);
			

			SimilarityInfo tmpSimilarityInfo;
			tmpSimilarityInfo.id = id;
			tmpSimilarityInfo.name = name;
			tmpSimilarityInfo.comment = comment;
			tmpSimilarityInfo.length = length;
			similarityInfos.push_back(tmpSimilarityInfo);
		}

		vector<EdgeInfo> forest = generateForest(mst, threshold);

		vector<vector<int> > cluster = generateCluster(forest, mst.size()+1);
		cerr << "start the output " << endl;
		
		for(int i = 0; i < cluster.size(); i++){
			printf("the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++){
				//if(cluster[i][j] > similarityInfos.size()) continue;
				cout << j << '\t' << cluster[i][j] << '\t' << similarityInfos[cluster[i][j]].name << endl;
				string tmpComment("");
				for(int k = 1; k < similarityInfos[cluster[i][j]].comment.length(); k++){
					if(similarityInfos[cluster[i][j]].comment[k] != '>'){
						tmpComment += similarityInfos[cluster[i][j]].comment[k];	
					}
					else{
						cout << "\t\t\t" << tmpComment << endl;
						tmpComment = "";
					}
				}
				if(tmpComment.length() != 0){
					cout << "\t\t\t" << tmpComment << endl;
					tmpComment = "";
				}

			}
			cout << endl;
		}

		return 0;//end main 

	}

	if(sketchByFile) cout << "sketch by file!" << endl;
	else cout << "sketch by sequence!" << endl;
	cout << "the thread number is: " << threads << endl;
	cout << "the threshold is: " << threshold << endl;
	



	//section 2: read the files and create sketches.
	//vector<Sketch::MinHash*> minHashes;
	vector<SketchInfo> sketches;
	vector<SimilarityInfo> similarityInfos;

	//for weighted MinHash

	double t0 = get_sec();
	int sketchSize_ = MINHASH_SKETCH_SIZE;
	cerr << "the sketchSize is: " << sketchSize_ << endl;

	if(!sketchByFile){
		if(!sketchSequences(inputFile, sketchFunc, similarityInfos, sketches, threads)){
			printUsage();
			return 1;
		}
	
	}
	else{//sketch by file

		if(!sketchFiles(inputFile, sketchFunc, similarityInfos, sketches, threads)){
			printUsage();
			return 1;
		}

	}//end sketch by file

	//make the minHashes[i].index == i;
	sort(sketches.begin(), sketches.end(), cmpSketch);

	//make the similarityInfos[i].id == i;
	sort(similarityInfos.begin(), similarityInfos.end(), cmpInfo);

	double t1 = get_sec();
	cerr << "time of format the sequence and create sketch is: " << t1-t0 << endl;

	//section 3: compute the distance matrix and create the graph.
	//no save the whole graph, create subMST when construct the subGraph for reducing memory footprint.

	//vector< vector<EdgeInfo> > graphArr;
	vector <vector<EdgeInfo> > mstArr;
	mstArr.resize(threads);
	int subSize = 8;
	cerr << "the size of sketches is: " << sketches.size() << endl;
	
	//int index = 0;
	int id = 0;
	int tailNum = sketches.size() % subSize;
	#pragma omp parallel for num_threads(threads) schedule (dynamic)
	for(id = 0; id < sketches.size() - tailNum; id+=subSize){
		int thread_id = omp_get_thread_num();
		for(int i = id; i < id+subSize; i++){
			//EdgeInfo tmpEdge;
			for(int j = i+1; j < sketches.size(); j++){
				double tmpDist;
				if(sketchFunc == "MinHash")
				{
					tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
					//tmpDist = 1.0 - sketches[i].minHash->jaccard(sketches[j].minHash);
				}
				else if(sketchFunc == "WMH"){
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				}
				else if(sketchFunc == "HLL"){
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				}
				else if(sketchFunc == "OMH"){
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				}
				else	
					break;
					
				EdgeInfo tmpE;
				tmpE.preNode = i;
				tmpE.sufNode = j;
				tmpE.dist = tmpDist;
				mstArr[thread_id].push_back(tmpE);
			}

		}
		sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
		vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
		mstArr[thread_id].swap(tmpMst);
	
		vector<EdgeInfo>().swap(tmpMst);
	}
	
	if(tailNum != 0){
		for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
			for(int j = i+1; j < sketches.size(); j++){
				//double tmpDist = 1.0 - minHashes[i].minHash->jaccard(minHashes[j].minHash);
				double tmpDist;
				if(sketchFunc == "MinHash"){
					tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
					//tmpDist = sketches[i].minHash->jaccard(sketches[j].minHash);
				}
				else if(sketchFunc == "WMH")
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				else if(sketchFunc == "HLL")
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				else if(sketchFunc == "OMH")
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				else	
					break;

				EdgeInfo tmpE;
				tmpE.preNode = i;
				tmpE.sufNode = j;
				tmpE.dist = tmpDist;
				mstArr[0].push_back(tmpE);
			}
		}
		if(mstArr[0].size() != 0){
			sort(mstArr[0].begin(), mstArr[0].end(), cmpEdge);
			vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[0], sketches.size());
			mstArr[0].swap(tmpMst);
		}

	}
	double t2 = get_sec();

	cerr << "time of computing distance and creating graph including sort the graph is: " << t2-t1<< endl;


	//section 4: generate the MST
	
//	//vector <vector<EdgeInfo> > mstArr;
//	for(int i = 0; i < graphArr.size(); i++){
//		vector<EdgeInfo> tmpMst = kruskalAlgorithm(graphArr[i], minHashes.size());
//		mstArr.push_back(tmpMst);
//	}

	double t3 = get_sec();
	cerr << "time of generating tmpMSTs is: " << t3-t2 << endl;

	vector<EdgeInfo> finalGraph;
	for(int i = 0; i < mstArr.size(); i++){
		finalGraph.insert(finalGraph.end(), mstArr[i].begin(), mstArr[i].end());
		vector<EdgeInfo>().swap(mstArr[i]);
		//mstArr[i].clear();
	}

	double t3_1 = get_sec();

	sort(finalGraph.begin(), finalGraph.end(), cmpEdge);

	double t3_2 = get_sec();

	vector<EdgeInfo> mst = kruskalAlgorithm(finalGraph, sketches.size());
	vector<EdgeInfo>().swap(finalGraph);

	double t4 = get_sec();
	cerr << "time of merge subGraph(subMST) and sort finalGraph and generate finalMST is: " << t4-t3 << endl;
	cerr << "\t time: merge subGraph(subMST) is: " << t3_1-t3 << endl;
	cerr << "\t time: sort merged graph is: " << t3_2-t3_1 << endl;
	cerr << "\t time: generage finalMST is: " << t4-t3_2 << endl;
	
	
	//save the matching of graph id and genomeInfo 
	ofstream ofile;
	ofile.open(inputFile+sketchFunc+"GenomeInfo");
	for(int i = 0; i < similarityInfos.size(); i++){
		ofile << similarityInfos[i].id << ' ' << similarityInfos[i].name << ' ' << ' ' << similarityInfos[i].length << endl;
		ofile << similarityInfos[i].comment << endl;
	}
	ofile.close();

	//save the mst
	ofstream ofile1;
	ofile1.open(inputFile+sketchFunc+"MSTInfo");
	for(int i = 0; i < mst.size(); i++){
		printf("<%d, %d, %lf>\t%s\t%s\n", mst[i].preNode, mst[i].sufNode, mst[i].dist, similarityInfos[mst[i].preNode].name.c_str(), similarityInfos[mst[i].sufNode].name.c_str());
		//printf("<%d, %d, %lf>\n", mst[i].preNode, mst[i].sufNode, mst[i].dist);
		ofile1 << mst[i].preNode << ' ' << mst[i].sufNode << ' ' << mst[i].dist << endl;
	}
	cout << endl;
	ofile1.close();

	double t5 = get_sec();

	//section 5; generate the clustering 
	
	vector<EdgeInfo> forest;
	forest = generateForest(mst, threshold);

	vector<vector<int> >cluster = generateCluster(forest, sketches.size());

	double t6 = get_sec();
	cerr << "time of generating forest and cluster is: " << t6 - t5 << endl;

	for(int i = 0; i < cluster.size(); i++){
		printf("the cluster %d is: \n", i);
		for(int j = 0; j < cluster[i].size(); j++){
			//if(cluster[i][j] > similarityInfos.size()) continue;
			//cout << j << '\t' << cluster[i][j] << '\t' << similarityInfos[cluster[i][j]].id << '\t' << similarityInfos[cluster[i][j]].name << endl;;
			cout << j << '\t' << cluster[i][j] << '\t' << similarityInfos[cluster[i][j]].name << endl;
			string comment("");
			for(int k = 1; k < similarityInfos[cluster[i][j]].comment.length(); k++){
				if(similarityInfos[cluster[i][j]].comment[k] != '>'){
					comment += similarityInfos[cluster[i][j]].comment[k];
				}
				else{
					cout << "\t\t\t" << comment << endl;
					comment = "";
				}
			}
			if(comment.length() > 0){
				cout << "\t\t\t" << comment << endl;
				comment = "";
			}
				
			//cout << "\t\t" << similarityInfos[cluster[i][j]].comment << endl;
		}
		cout << endl;
	}


	return 0;
}//end main





