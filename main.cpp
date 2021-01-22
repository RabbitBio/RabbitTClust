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

		Sketch::MinHash *mh1 = new Sketch::MinHash();//the sketchSize(number of hashes) need to be passed properly.
		//the memcopy from Sketch::MinHash::reference to SketchInfo needs optimizing.
		mh1->update(ks1->seq.s);
		//store the hashes to the SketchInfo.
		
		//tmpSketchInfo.size = mh1->getSketchSize();
		//TODO: change the RabbitSketch src file of Sketch::MinHash::printSketch.
		//vector<uint64_t> hashes = mh1->printMinHashes();
		//tmpSketchInfo.hashes = hashes;

		//sketchInfos.push_back(tmpSketchInfo);
		minHashes.push_back(mh1);

		similarityInfos.push_back(tmpSimilarityInfo); 

		index++;
		//if(index == 5) break;

			

	}//end while

	cerr << "the number of the sequence is: " << index << endl;

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
	int threads = 32;

	Graph tmpG[threads];

	#pragma omp parallel for num_threads(threads)
	for(int i = 0; i < minHashes.size(); i++){
		//Graph tmpG;
		int thread_id = omp_get_thread_num();
		#pragma omp critical
		{
		cerr << "the thread id is: " << thread_id << endl;
		}
		tmpG[thread_id].node = i;
		for(int j = i+1; j < minHashes.size(); j++){
			double tmpDist = 1.0 - minHashes[i]->jaccard(minHashes[j]);
			//printf("<%d, %d, %lf>\n", i, j, tmpDist);
			tmpG[thread_id].neighbor.push_back(NeighborNode(j, tmpDist));
		}
		//sort the suffixNodes by weight of the edge in an ascending order.
		std::sort(tmpG[thread_id].neighbor.begin(), tmpG[thread_id].neighbor.end(), cmpNeighbor);

		#pragma omp critical
		{
		graphs.push_back(tmpG[thread_id]);
		tmpG[thread_id].neighbor.clear();
		//cerr << "end the sequence: " << i << endl;
		}

	}

	cerr << "the size of the graphs is: " << graphs.size() << endl;

	MST mst;
	creatMST(mst, graphs, index);
	cout << endl;
	cout << endl;
	cout << endl;

	cout << "print the srcMST ================================================" << endl;
	printMST(mst);


	MST resForest;
	cout << "the size of mst.node is: " << mst.nodes.size() << endl;
	creatForest(mst, resForest, 0.9);
	cout << "the size of resForest.node is: " << resForest.nodes.size() << endl;
	cout << "print the resForest =============================================" << endl;
	printMST(resForest);


	vector<Graph> resGraph;
	mst2Graph(resForest, resGraph);

	cout << "print the resGraph =============================================" << endl;
	printGraph(resGraph);
	

	//vector< vector<int> > clusters;
	vector< unordered_set<int> > clusters;
	creatClust(resGraph, clusters);
	

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



	


	

	



	
		



	return 0;
}//end main






























