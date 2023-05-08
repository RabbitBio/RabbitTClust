#ifndef HPP_COMMON
#define HPP_COMMON

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
#include <iostream>

//#define KMER_SIZE 21
//#define MINHASH_SKETCH_SIZE 10000
//#define SKETCH_COMPRESS_SEQUENCE 1000
//#define SKETCH_COMPRESS_GENOME 1000
#define WMH_SKETCH_SIZE 50
#define WINDOW_SIZE 20
#define HLL_SKETCH_BIT 10
#define DENSE_SPAN 100;

#include <sys/time.h>
inline double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

#include <ctime>
inline const string currentDataTime(){
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y_%m_%d_%H-%M-%S", &tstruct);

	return buf;
}

#endif
