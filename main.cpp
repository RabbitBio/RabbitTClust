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
 * section 5: generate the clusters with the MST using different distance threshhold.
 *
 * Author: Xiaoming Xu
 * Mar 5, 2021
 *
 */

int main(int argc, char * argv[]){

	//section 1: init parameters
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



	//section 2: read the files and create sketches.
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


			//Sketch::MinHash *mh1 = new Sketch::MinHash();//the sketchSize(number of hashes) need to be passed properly.
			Sketch::MinHash *mh1 = new Sketch::MinHash(21, 10000);
			//mh1->update((char *)similarityInfos[i].seq.c_str());
			mh1->update(ks1->seq.s);
			minHashes.push_back(mh1);
	
			index++;
			//if(index == 5) break;
	
		}//end while
		cerr << "the number of the sequence is: " << index << endl;
	
	}
	else{//sketch by file
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

	}//end sketch by file
	double t1 = get_sec();
	cerr << "the time of format the sequence and create sketch is: " << t1-t0 << endl;
		

	//section 3: compute the distance matrix and create the graph.

	vector<EdgeInfo> graph;
	cerr << "the size of minHashes is: " << minHashes.size() << endl;
	cerr << "the capacity of minHashes is: " << minHashes.capacity() << endl;

	for(int i = 0; i < minHashes.size(); i++){
		//EdgeInfo tmpEdge;
		for(int j = i+1; j < minHashes.size(); j++){
			double tmpDist = 1.0 - minHashes[i]->jaccard(minHashes[j]);
			graph.push_back(EdgeInfo(i, j, tmpDist));
		}

	}
	double t2 = get_sec();
	std::sort(graph.begin(), graph.end(), cmpEdge);

	double t3 = get_sec();

	cerr << "the size of the graphs is: " << graph.size() << endl;
	cerr << "the time of computing distance and creating graph is: " << t2-t1<< endl;
	cerr << "the time of sort graph is: " << t3-t2<< endl;

	//section 4: generate the MST
	//MST mst;
	vector<EdgeInfo> mst = kruskalAlgorithm(graph, minHashes.size());

	double t4 = get_sec();
	cerr << "the time of createMST is: " << t4-t3 << endl;
	printMST(mst);

	//double t5 = get_sec();


	return 0;
}//end main






























