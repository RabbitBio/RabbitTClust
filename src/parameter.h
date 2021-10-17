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
#define SKETCH_COMPRESS_SEQUENCE 1000
#define SKETCH_COMPRESS_GENOME 1000
#define WMH_SKETCH_SIZE 50
#define WINDOW_SIZE 20
#define HLL_SKETCH_BIT 10

inline void printUsage(){
	fprintf(stdout, "usage: clust [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> genomeFile\n");
	fprintf(stdout, "usage: clust [-h] [-f] [-d] <double> genomeInfo mstInfo\n");
	fprintf(stdout, "	-h		: this help message\n");
	fprintf(stdout, "	-k		: set kmer size, default 21\n");
	fprintf(stdout, "	-s		: set sketch size, default 1000\n");
	fprintf(stdout, "	-l		: genome clustering, inputFile is the path list of the genome files\n");
	fprintf(stdout, "	-c		: compute the containment of sequences(genomes), need sketch function MinHash\n");
	fprintf(stdout, "	-t		: set the thread number\n");
	fprintf(stdout, "	-d		: set the threshold of the clusters from the Minimum Spanning Tree\n");
	fprintf(stdout, "	-f		: input files are genomeInfo and MST content\n");
	fprintf(stdout, "	-F		: sketch function, includes MinHash, WMH, OMH, HLL\n"); 
	fprintf(stdout, "	-o		: output file \n"); 
	fprintf(stdout, "example as follows:\n");
	fprintf(stdout, "\tfor genomes clustering, input are genome file list:\n");
	fprintf(stdout, "\t ./clust -l -d 0.05 -t 48 -F MinHash -o result.out viralList\n");
	fprintf(stdout, "\tfor generator cluster from exist MST:\n");
	fprintf(stdout, "\t ./clust -f -d 0.05 viralListMinHashGenomeInfo viralListMinHashMSTInfo -o result.out\n");
}




#endif
