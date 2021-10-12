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
#define SKETCH_COMPRESS_SEQUENCE 100
#define SKETCH_COMPRESS_GENOME 1000
#define WMH_SKETCH_SIZE 50
#define WINDOW_SIZE 20
#define HLL_SKETCH_BIT 10


void printUsage();





#endif
