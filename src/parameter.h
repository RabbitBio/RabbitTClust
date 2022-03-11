#ifndef H_PARAMETER
#define H_PARAMETER

/*
 * The parameter.h is for the basic parameter of sketch, such as sketchSize, kmerSize.
 * 	KMER_SIZE is the kmer size for slide-window of all sketch functions.
 * 	MINHASH_SKETCH_SIZE is the fixed sketch size for minHash resemblance computing.
 * 	SKETCH_COMPRESS_SEQUENCE is the proportion sketch size with genome sequences for containment computing.
 * 	SKETCH_COMPRESS_GENOME is the proportion sketch size with genome size for containment computing.
 * 	WMH_SKETCH_SIZE is the sketch size for WeightedMinHash.
 * 	WINDOW_SIZE is the window size for minimizer in WeightedMinHash.
 * 	HLL_SKETCH_BIT is the bit number to define the sketch size for HyperLogLog.
 *
 */

#include <cstdio>

//#define KMER_SIZE 21
//#define MINHASH_SKETCH_SIZE 10000
//#define SKETCH_COMPRESS_SEQUENCE 1000
//#define SKETCH_COMPRESS_GENOME 1000
#define WMH_SKETCH_SIZE 50
#define WINDOW_SIZE 20
#define HLL_SKETCH_BIT 10

#include <sys/time.h>
inline double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

inline void printUsage(){
	fprintf(stdout, "usage: clust [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> -i <string> \n");
	fprintf(stdout, "usage: clust [-h] [-f] [-d] <double> -i <string> <string> -o <string>\n");
	fprintf(stdout, "	-h\t\t: this help message\n");
	fprintf(stdout, "	-k <int>\t: set kmer size, default 21\n");
	fprintf(stdout, "	-s <int>\t: set sketch size, default 1000\n");
	fprintf(stdout, "	-l\t\t: cluster for genomes(not sequences), list input. Lines in each <input> specify paths to genome files, one per line.\n");
	fprintf(stdout, "	-c <int>\t: compute the containment of genomes, set proportion sketchSize = genomeSize/compress, ATTENTION with MinHash function\n");
	fprintf(stdout, "	-t <int>\t: set the thread number, default 1\n");
	fprintf(stdout, "	-d <double>\t: set the threshold of the clusters from the Minimum Spanning Tree\n");
	fprintf(stdout, "	-f\t\t: input files are genomeInfo and MST contents\n");
	fprintf(stdout, "	-F <string>\t: sketch function, includes MinHash, WMH, OMH, HLL, default MinHash\n"); 
	fprintf(stdout, "	-o <string>\t: path of output file \n"); 
	fprintf(stdout, "	-i <string>\t: path of input file, ATTENTION with -f option \n"); 

	fprintf(stdout, "example as follows:\n");
	fprintf(stdout, "\tfor genomes clustering, input is a genome file list:\n");
	fprintf(stdout, "\t ./clust -l -d 0.05 -t 48 -F MinHash -i data/bacteriaList -o bacteria.clust\n");

	fprintf(stdout, "\tfor genomes clustering, input is a single genome file:\n");
	fprintf(stdout, "\t ./clust -d 0.05 -t 48 -F MinHash -i data/bacteria.fna -o bacteria.clust\n");

	fprintf(stdout, "\tfor redundancy detection with containment, input is a genome file list:\n");
	fprintf(stdout, "\t ./clust -l -c 10000 -d 0.10 -t 48 -F MinHash -i data/bacteriaList -o bacteria.out\n");

	fprintf(stdout, "\tfor redundancy detection with containment, input is a single genome file:\n");
	fprintf(stdout, "\t ./clust -c 10000 -d 0.10 -t 48 -F MinHash -i data/bacteria.fna -o bacteria.out\n");

	fprintf(stdout, "\tfor generator cluster from exist MST:\n");
	fprintf(stdout, "\t ./clust -f -d 0.05 -i bacteriaListMinHashGenomeInfo bacteriaListMinHashMSTInfo -o bacteria.clust\n");
	fprintf(stdout, "\t ATTENTION: the -f must in front of the -i option\n");

}




#endif
