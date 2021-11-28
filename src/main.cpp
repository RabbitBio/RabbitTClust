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
 * section 3: compute the distance matrix and generate the Minimum Spanning Tree(MST).
 * section 4: generate the clusters with the MST using different distance threshold.
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
#include "MST.h"
#include <omp.h>
#include "UnionFind.h"
#include <algorithm>
#include "parameter.h"
#include "MST_IO.h"
#include <math.h>

using namespace std;


double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

int main(int argc, char * argv[]){

	//section 1: init parameters
	int argIndex = 1;
	string inputFile = "genome.fna";
	string inputFile1 = "genome.info";
	string sketchFunc = "MinHash";
	string outputFile = "result.out";
	int threads = 1;
	bool sketchByFile = false;
	bool isContainment = false;
	bool useMST = false;
	double threshold = 0.05;
	int kmerSize = 21;
	int sketchSize = 1000;
	while(argIndex < argc){
		//if(useMST && argIndex + 2 == argc){
		//	inputFile1 = argv[argIndex];
		//}
		//else if(argIndex + 1 == argc){
		//	if(argv[argIndex][1] == 'h'){
		//		printUsage();
		//		return 1;
		//	}
		//	inputFile = argv[argIndex];
		//}
		//else{
		switch(argv[argIndex][1]){
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
			case 'c':
				isContainment = true;
				//threshold = 1 - exp(-threshold * KMER_SIZE) / (2 - exp(-threshold * KMER_SIZE));//generator from mash distance
				threshold = 0.30;
				fprintf(stderr, "compute containment\n");
				break;
			case 'k':
				kmerSize = stoi(argv[++argIndex]);
				fprintf(stderr, "set kmerSize: %d\n", kmerSize);
				break;
			case 's':
				sketchSize = stoi(argv[++argIndex]);
				fprintf(stderr, "set sketchSize:  %d\n", sketchSize);
				break;
			case 'd':
				threshold = stod(argv[++argIndex]);
				fprintf(stderr, "set the threshold: %lf \n", threshold);
				break;
			case 'f':
				useMST = true;
				fprintf(stderr, "input as Munimum Spanning Tree \n");
				break;
			case 'F':
				sketchFunc = argv[++argIndex];
				fprintf(stderr, "the sketch function is: %s \n", sketchFunc.c_str());
				break;
			case 'o':
				outputFile = argv[++argIndex];
				fprintf(stderr, "set output file: %s \n", outputFile.c_str());
				break;
			case 'i':
			{
				if(!useMST)
				{
					inputFile = argv[++argIndex];
					fprintf(stderr, "the inputFile is: %s \n", inputFile.c_str());
				}
				else
				{
					inputFile1 = argv[++argIndex];
					if(argIndex == argc)
					{
						fprintf(stderr, "error input File with -f options, exit\n");
						printUsage();
						return 1;
					}
					inputFile = argv[++argIndex];
					fprintf(stderr, "the genomeInfo and MSTInfo are: %s and %s\n", inputFile1.c_str(), inputFile.c_str());
				}
				break;
			}
			default:
				fprintf(stderr, "Invalid option %s\n", argv[argIndex]);
				printUsage();
				return 1;
		}
		//}
		++argIndex;
	}//end while argument parse;

	//input as GenomeInfo and MSTInfo
	if(useMST){
		MST2Cluster(inputFile, inputFile1, outputFile, threshold);
		return 0;//end main 
	}

	if(sketchByFile) cout << "sketch by file!" << endl;
	else cout << "sketch by sequence!" << endl;

	
	int sketchSize_;
	if(isContainment)
		sketchSize_ = sketchByFile ? SKETCH_COMPRESS_GENOME : SKETCH_COMPRESS_SEQUENCE;
	else 
		sketchSize_ = sketchSize;
	//int kmerSize_ = KMER_SIZE;
#ifdef DEBUG
	cerr << "the kmerSize is: " << kmerSize << endl;
	cerr << "the thread number is: " << threads << endl;
	cerr << "the threshold is: " << threshold << endl;
	if(isContainment)
		cerr << "the sketchSize is in proportion with 1/" << sketchSize_ << endl;
	else
		cerr << "the sketchSize is: " << sketchSize_ << endl;
#endif
	
	//section 2: read the files and create sketches.
	vector<SketchInfo> sketches;

#ifdef Timer
	double t0 = get_sec();
#endif

	if(!sketchByFile){
		if(!sketchSequences(inputFile, kmerSize, sketchSize, sketchFunc, isContainment, sketches, threads)){
			printUsage();
			return 1;
		}
	
	}//end sketch by sequence
	else{
		if(!sketchFiles(inputFile, kmerSize, sketchSize, sketchFunc, isContainment, sketches, threads)){
			printUsage();
			return 1;
		}
	}//end sketch by file

	cerr << "the size of sketches(number of genomes or sequences) is: " << sketches.size() << endl;
#ifdef Timer
	double t1 = get_sec();
	cerr << "========time of sketch is: " << t1 - t0 << "========" << endl;
#endif

	//section 3: compute the distance matrix and generate MST
	vector<EdgeInfo> mst = generateMST(sketches, sketchFunc, threads);
#ifdef Timer
	double t2 = get_sec();
	cerr << "========time of generateMST is: " << t2 - t1 << "========" << endl;
#endif

	//save the matching of graph id and genomeInfo 
	saveMST(inputFile, sketchFunc, isContainment, sketches, mst, sketchByFile, sketchSize_, kmerSize);

#ifdef Timer
	double t3 = get_sec();
	cerr << "========time of saveMST is: " << t3 - t2 << "========" << endl;
#endif

	//section 4: generate the clustering 
	
	vector<EdgeInfo> forest = generateForest(mst, threshold);

	vector<vector<int> >cluster = generateCluster(forest, sketches.size());

#ifdef Timer
	double t4 = get_sec();
	cerr << "========time of generator forest and cluster is: " << t4 - t3 << "========" << endl;
#endif

	printResult(cluster, sketches, sketchByFile, outputFile);
#ifdef Timer
	double t5 = get_sec();
	cerr << "========time of save result is: " << t5 - t4 << "========" << endl;
#endif

	return 0;
}//end main





