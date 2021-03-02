#include <iostream>
#include "SketchInfo.h"
#include "Sketch.h"// need to add the include path in Makefile.
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <stdio.h>
#include "MST.h"
#include <omp.h>
#include <fstream>
#include <string>

KSEQ_INIT(gzFile, gzread);

using namespace std;


double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

void printUsage(void){
	fprintf(stdout, "usage clust [-h] [-l] [-t] inputFile\n");
	fprintf(stdout, "	-h		: this help message\n");
	fprintf(stdout, "	-l		: genome clustering, inputFile is the path list of the genome files\n");
	fprintf(stdout, "	-t		: genome clustering, set the thread number\n");
	
}


int main(int argc, char * argv[]){


	int argIndex = 1;
	string inputFile = "genome.fna";
	int threads = 1;
	bool sketchByFile = false;
	while(argIndex < argc){
		if(argIndex + 1 == argc){
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
				default:
					fprintf(stderr, "Invalid option %s\n", argv[argIndex]);
					printUsage();
					return 1;
			}
		}
		++argIndex;
	}//end while file;


	//vector<SketchInfo> sketchInfos;
	vector<Sketch::MinHash*> minHashes;
	vector<SimilarityInfo> similarityInfos;

	double t0 = get_sec();

	if(!sketchByFile){
		fprintf(stderr, "input one file, sketch by sequence: \n");

		gzFile fp1;
		kseq_t* ks1;
	
		fp1 = gzopen(inputFile.c_str(), "r");
		if(fp1 == NULL){
			fprintf(stderr, "cannot open the genome file\n");
			return 1;
		}
	
		ks1 = kseq_init(fp1);
	
	
		int index = 0;
		while(1){
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
	
			//SketchInfo tmpSketchInfo;
			SimilarityInfo tmpSimilarityInfo;
		
			tmpSimilarityInfo.id = index;
			tmpSimilarityInfo.name = ks1->name.s;
			tmpSimilarityInfo.comment = ks1->comment.s;
			tmpSimilarityInfo.length = length;
			//tmpSimilarityInfo.seq = ks1->seq.s;//do not store the seq content;
	
			similarityInfos.push_back(tmpSimilarityInfo); 


			Sketch::MinHash *mh1 = new Sketch::MinHash();//the sketchSize(number of hashes) need to be passed properly.
			//mh1->update((char *)similarityInfos[i].seq.c_str());
			mh1->update(ks1->seq.s);
			minHashes.push_back(mh1);
	
			index++;
			//if(index == 5) break;
	
		}//end while
		cerr << "the number of the sequence is: " << index << endl;
	
	}
	else{
		fprintf(stderr, "input fileList, sketch by file\n");
		fstream fs(inputFile);
		if(!fs){
			fprintf(stderr, "error open the inputFile of fileList\n");
			return 1;
		}
		vector<string> fileList;
		string fileName;
		while(getline(fs, fileName)){
			fileList.push_back(fileName);
		} 
		
		for(int i = 0; i < fileList.size(); i++){
			//cout << fileList[i] << endl;
			gzFile fp1;
			kseq_t* ks1;
			
			fp1 = gzopen(fileList[i].c_str(), "r");
			if(fp1 == NULL){
				fprintf(stderr, "cannot open the genome file\n");
				return 1;
			}
	
			ks1 = kseq_init(fp1);

			Sketch::MinHash *mh1 = new Sketch::MinHash(21, 10000);

			int totalLength = 0;
			

			while(1){
				int length = kseq_read(ks1);
				totalLength += length;
				if(length < 0){
					break;
				}

				mh1->update(ks1->seq.s);




			}
			//cerr << "finish the file: 	" << fileList[i] << endl;
			minHashes.push_back(mh1);

			SimilarityInfo tmpSimilarityInfo;
			tmpSimilarityInfo.id = i;
			tmpSimilarityInfo.name = fileList[i];
			tmpSimilarityInfo.length = totalLength;
			similarityInfos.push_back(tmpSimilarityInfo);
		}







	}
	double t1 = get_sec();
	cerr << "the time of format the sequence is: " << t1-t0 << endl;

	//minHashes.reserve(similarityInfos.size());
	//#pragma omp parallel for num_threads(threads) schedule(dynamic)
	//for(int i = 0; i < similarityInfos.size(); i++){
	//	Sketch::MinHash *mh1 = new Sketch::MinHash();//the sketchSize(number of hashes) need to be passed properly.
	//	//Sketch::MinHash *mh1 = new Sketch::MinHash(21, 100000, 42, true);//the sketchSize(number of hashes) need to be passed properly.
	//	//the memcopy from Sketch::MinHash::reference to SketchInfo needs optimizing.
	//	//mh1->update(ks1->seq.s);
	//	//cout << similarityInfos[i].seq << endl;
	//	mh1->update((char *)similarityInfos[i].seq.c_str());
	//	
	//	//#pragma omp critical
	//	{
	//	//minHashes.push_back(mh1);
	//	minHashes[similarityInfos[i].id] = mh1;
	//	}
	//}
		
	double t2 = get_sec();

	cerr << "the time of create minHash sketches MultiThread is: " << t2-t1 << endl;

//	for(int i = 0; i < sketchInfo.size(); i++){
//		cout << sketchInfo[i].name << endl;
//		cout << sketchInfo[i].comment << endl;
//		cout << sketchInfo[i].length << endl;
//		for(int j = 0; j < sketchInfo[i].size; j++){
//			cout << sketchInfo[i].hashes[j] << endl;
//		}
//
//	}


	//create the graphs.
	vector<Graph> graphs;
	//graphs.reserve(minHashes.capacity());
	//graphs.reserve(similarityInfos.size());

	//Graph tmpG[threads];
	cerr << "the size of minHashes is: " << minHashes.size() << endl;
	cerr << "the capacity of minHashes is: " << minHashes.capacity() << endl;
	//return 0;

	//#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < minHashes.size(); i++){
	//for(int i = 0; i < minHashes.capacity(); i++){
		Graph tmpG;
		//int thread_id = omp_get_thread_num();
		//#pragma omp critical
		//{
		//cerr << "the thread id is: " << thread_id << endl;
		//}
		tmpG.node = i;
		tmpG.curNeighbor = 0;
		//printf("print the minHashes of the sequence %d\n", i);
		//minHashes[i]->printMinHashes();
		for(int j = i+1; j < minHashes.size(); j++){
		//for(int j = i+1; j < minHashes.capacity(); j++){
			//printf("print the minHashes of the sequence %d\n", j);
			//minHashes[j]->printMinHashes();
			double tmpDist = 1.0 - minHashes[i]->jaccard(minHashes[j]);
			//printf("<%d, %d, %lf>\n", i, j, tmpDist);
			tmpG.neighbor.push_back(NeighborNode(j, tmpDist));
		}
		//sort the suffixNodes by weight of the edge in an ascending order.
		std::sort(tmpG.neighbor.begin(), tmpG.neighbor.end(), cmpNeighbor);
		//cerr << "end the sort : " << endl;

		//#pragma omp critical
		//{
		graphs.push_back(tmpG);
		//graphs[i] = tmpG[thread_id];
		tmpG.neighbor.clear();
		//cerr << "end the sequence: " << i << endl;
		//exit(0);
		//}

	}

	double t3 = get_sec();

	cerr << "the size of the graphs is: " << graphs.size() << endl;
	cerr << "the capacity of the graphs is: " << graphs.capacity() << endl;
	cerr << "the time of create the graph is: " << t3-t2<< endl;

	MST mst;
	//creatMST(mst, graphs, index);
	creatMST(mst, graphs,minHashes.size());
	cout << endl;
	cout << endl;
	cout << endl;

	double t4 = get_sec();
	cout << "print the srcMST ================================================" << endl;
	cerr << "the time of createMST is: " << t4-t3 << endl;
	printMST(mst);

	double t5 = get_sec();
	cerr << "the time of printMST is: " << t5-t4 << endl;

	MST resForest;
	cout << "the size of mst.node is: " << mst.nodes.size() << endl;
	creatForest(mst, resForest, 1.0);
	double t6 = get_sec();
	cerr << "the time of creatForest is: " << t6-t5 << endl;
	cout << "the size of resForest.node is: " << resForest.nodes.size() << endl;
	cout << "print the resForest =============================================" << endl;
	printMST(resForest);
	double t7 = get_sec();
	cerr << "the time of printMST is: " << t7-t6 << endl;


	vector<Graph> resGraph;
	mst2Graph(resForest, resGraph);
	double t8 = get_sec();
	cerr << "the time of mst2Graph is: " << t8-t7 << endl;

	cout << "print the resGraph =============================================" << endl;
	printGraph(resGraph);
	
	double t9 = get_sec();
	cerr << "the time of printGraph is: " << t9-t8 << endl;

	//vector< vector<int> > clusters;
	vector< unordered_set<int> > clusters;
	creatClust(resGraph, clusters);
	
	double t10 = get_sec();
	cerr << "the time of creatClust is: " << t10-t9 << endl;

	cout << "print the result Cluster ========================================" << endl;
	for(int i = 0; i < clusters.size(); i++){
		printf("the cluster %d is: \n", i);
		//for(int j = 0; j < clusters[i].size(); j++){
		for(const int &x : clusters[i]){
			//printf("%d\t", clusters[i][j]);
			printf("%d\t", x);
		}
		printf("\n");
	}
	double t11 = get_sec();
	cerr << "the time of print the resultClust is: " << t11-t10 << endl;



	


	

	



	
		



	return 0;
}//end main






























