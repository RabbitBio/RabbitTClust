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
#include <zlib.h>
#include "MST.h"
#include <omp.h>
#include "UnionFind.h"
#include <algorithm>
#include "parameter.h"
#include "MST_IO.h"
#include <math.h>

#ifdef GREEDY_CLUST
#include "Sketch_IO.h"
#include "greedy.h"
#endif

using namespace std;

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
	bool useFile = false;
	double threshold = 0.05;
	int kmerSize = 21;
	int sketchSize = 1000;
	int containCompress = 10000;
	while(argIndex < argc){
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
				containCompress = stoi(argv[++argIndex]);
				threshold = 0.10;
				fprintf(stderr, "compute containment, The sketchSize is in proportion with 1/%d \n", containCompress);
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
				useFile = true;
				#ifdef GREEDY_CLUST
				fprintf(stderr, "input as sketches informations for greedy incremental cluster\n");
				#else
				fprintf(stderr, "input as MST for Minimum Spanning Tree cluster\n");
				#endif
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
				if(!useFile)
				{
					inputFile = argv[++argIndex];
					fprintf(stderr, "the inputFile is: %s \n", inputFile.c_str());
				}
				else
				{
					inputFile = argv[++argIndex];
					if(argIndex == argc)
					{
						fprintf(stderr, "error input File with -f options, exit\n");
						printUsage();
						return 1;
					}
					inputFile1 = argv[++argIndex];
					#ifdef GREEDY_CLUST
					fprintf(stderr, "the genomeInfo and SketchInfo are: %s and %s\n", inputFile.c_str(), inputFile1.c_str());
					#else
					fprintf(stderr, "the genomeInfo and MSTInfo are: %s and %s\n", inputFile.c_str(), inputFile1.c_str());
					#endif
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


#ifdef GREEDY_CLUST
	cerr << "use the Greedy cluster" << endl;
#else
	cerr << "use the MST cluster" << endl;
#endif


	//input as GenomeInfo and MSTInfo or sketch informations
	if(useFile){
#ifdef GREEDY_CLUST
		Sketch2Clust(inputFile, inputFile1, outputFile, threshold, threads);
#else
		MST2Cluster(inputFile, inputFile1, outputFile, threshold);
#endif
		return 0;//end main 
	}

	if(sketchByFile) cout << "sketch by file!" << endl;
	else cout << "sketch by sequence!" << endl;

	
#ifdef DEBUG
	cerr << "the kmerSize is: " << kmerSize << endl;
	cerr << "the thread number is: " << threads << endl;
	cerr << "the threshold is: " << threshold << endl;
	if(isContainment)
		cerr << "the sketchSize is in proportion with 1/" << containCompress << endl;
	else
		cerr << "the sketchSize is: " << sketchSize << endl;
#endif
	
	//section 2: read the files and create sketches.
	vector<SketchInfo> sketches;

#ifdef Timer
	double t0 = get_sec();
#endif

	if(!sketchByFile){
		if(!sketchSequences(inputFile, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){
			printUsage();
			return 1;
		}
	
	}//end sketch by sequence
	else{
		if(!sketchFiles(inputFile, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){
			printUsage();
			return 1;
		}
	}//end sketch by file

	cerr << "the size of sketches(number of genomes or sequences) is: " << sketches.size() << endl;
#ifdef Timer
	double t1 = get_sec();
	cerr << "========time of sketch is: " << t1 - t0 << "========" << endl;
#endif


#ifdef GREEDY_CLUST
	saveSketches(sketches, inputFile, sketchFunc, isContainment, containCompress, sketchByFile, sketchSize, kmerSize);
#ifdef Timer
	double t2 = get_sec();
	cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
	vector<vector<int> >cluster = greedyCluster(sketches, sketchFunc, threshold, threads);
#ifdef Timer
	double t3 = get_sec();
	cerr << "========time of greedyCluster is: " << t3 - t2 << "========" << endl;
#endif
	printResult(cluster, sketches, sketchByFile, outputFile);

#else
	//section 3: compute the distance matrix and generate MST
	vector<EdgeInfo> mst = generateMST(sketches, sketchFunc, threads);
#ifdef Timer
	double t2 = get_sec();
	cerr << "========time of generateMST is: " << t2 - t1 << "========" << endl;
#endif

	//save the matching of graph id and genomeInfo 
	saveMST(inputFile, sketchFunc, isContainment, containCompress, sketches, mst, sketchByFile, sketchSize, kmerSize);

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
#endif//endif GREEDY_CLUST

	return 0;
}//end main





