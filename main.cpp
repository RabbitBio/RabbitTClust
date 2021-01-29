#include <iostream>
#include "SketchInfo.h"
#include "Sketch.h"// need to add the include path in Makefile.
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <stdio.h>
#include "MST.h"
#include <omp.h>

KSEQ_INIT(gzFile, gzread);

using namespace std;


double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}


int main(int argc, char * argv[]){

	gzFile fp1;
	kseq_t* ks1;

	fp1 = gzopen(argv[1], "r");
	if(fp1 == NULL){
		fprintf(stderr, "cannot open the input file\n");
		return 1;
	}

	ks1 = kseq_init(fp1);
	

	//vector<SketchInfo> sketchInfos;
	vector<Sketch::MinHash*> minHashes;
	vector<SimilarityInfo> similarityInfos;

	double t0 = get_sec();

	int index = 0;
	while(1){
		int length = kseq_read(ks1);
		if(length < 0){
			break;
		}

		//SketchInfo tmpSketchInfo;
		SimilarityInfo tmpSimilarityInfo;

		//tmpSketchInfo.id = index;
		//tmpSketchInfo.name = ks1->name.s;
		//tmpSketchInfo.comment = ks1->comment.s;
		//tmpSketchInfo.length = length;
		////not consider the fastq strand property.
		
	
		tmpSimilarityInfo.id = index;
		tmpSimilarityInfo.name = ks1->name.s;
		tmpSimilarityInfo.comment = ks1->comment.s;
		tmpSimilarityInfo.length = length;
		tmpSimilarityInfo.seq = ks1->seq.s;

		

		similarityInfos.push_back(tmpSimilarityInfo); 

		index++;
		//if(index == 5) break;

			

	}//end while

	cerr << "the number of the sequence is: " << index << endl;
	double t1 = get_sec();
	cerr << "the time of format the sequence is: " << t1-t0 << endl;

	int threads = 32;

	minHashes.reserve(similarityInfos.size());

	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < similarityInfos.size(); i++){

		Sketch::MinHash *mh1 = new Sketch::MinHash();//the sketchSize(number of hashes) need to be passed properly.
		//Sketch::MinHash *mh1 = new Sketch::MinHash(21, 100000, 42, true);//the sketchSize(number of hashes) need to be passed properly.
		//the memcopy from Sketch::MinHash::reference to SketchInfo needs optimizing.
		//mh1->update(ks1->seq.s);
		//cout << similarityInfos[i].seq << endl;
		mh1->update((char *)similarityInfos[i].seq.c_str());
		
		//#pragma omp critical
		{
		//minHashes.push_back(mh1);
		minHashes[similarityInfos[i].id] = mh1;
		}
	}
		
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
	graphs.reserve(similarityInfos.size());

	Graph tmpG[threads];
	cerr << "the size of minHashes is: " << minHashes.size() << endl;
	cerr << "the capacity of minHashes is: " << minHashes.capacity() << endl;
	//return 0;

	#pragma omp parallel for num_threads(1) schedule(dynamic)
	//for(int i = 0; i < minHashes.size(); i++){
	for(int i = 0; i < minHashes.capacity(); i++){
		//Graph tmpG;
		int thread_id = omp_get_thread_num();
		//#pragma omp critical
		//{
		//cerr << "the thread id is: " << thread_id << endl;
		//}
		tmpG[thread_id].node = i;
		printf("print the minHashes of the sequence %d\n", i);
		minHashes[i]->printMinHashes();
		//for(int j = i+1; j < minHashes.size(); j++){
		for(int j = i+1; j < minHashes.capacity(); j++){
			printf("print the minHashes of the sequence %d\n", j);
			minHashes[j]->printMinHashes();
			double tmpDist = 1.0 - minHashes[i]->jaccard(minHashes[j]);
			printf("<%d, %d, %lf>\n", i, j, tmpDist);
			tmpG[thread_id].neighbor.push_back(NeighborNode(j, tmpDist));
		}
		//sort the suffixNodes by weight of the edge in an ascending order.
		std::sort(tmpG[thread_id].neighbor.begin(), tmpG[thread_id].neighbor.end(), cmpNeighbor);
		//cerr << "end the sort : " << endl;

		#pragma omp critical
		{
		graphs.push_back(tmpG[thread_id]);
		//graphs[i] = tmpG[thread_id];
		tmpG[thread_id].neighbor.clear();
		cerr << "end the sequence: " << i << endl;
		exit(0);
		}

	}

	double t3 = get_sec();

	cerr << "the size of the graphs is: " << graphs.size() << endl;
	cerr << "the capacity of the graphs is: " << graphs.capacity() << endl;
	cerr << "the time of create the graph is: " << t3-t2<< endl;

	MST mst;
	creatMST(mst, graphs, index);
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






























