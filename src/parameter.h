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

#include <ctime>
inline const string currentDataTime(){
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y_%m_%d_%H-%M-%S", &tstruct);

	return buf;
}

inline void printUsage(){
	fprintf(stdout, "usage: clust-mst [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> -i <string> \n");
	fprintf(stdout, "usage: clust-mst [-h] [-f] [-E] [-d] <double> -i <string> <string> -o <string>\n");
	fprintf(stdout, "usage: clust-greedy [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> -i <string> \n");
	fprintf(stdout, "usage: clust-greedy [-h] [-f] [-d] <double> -i <string> <string> -o <string>\n");
	fprintf(stdout, "	-h\t\t: this help message\n");
	fprintf(stdout, "	-k <int>\t: set kmer size, default 21, for both clust-greedy and clust-mst\n");
	fprintf(stdout, "	-s <int>\t: set sketch size, default 1000, for both clust-greedy and clust-mst\n");
	fprintf(stdout, "	-l\t\t: cluster for genomes(not sequences), list input. Lines in each <input> specify paths to genome files, one per line. for both clust-greedy and clust-mst\n");
	fprintf(stdout, "	-c <int>\t: compute the containment of genomes, set proportion sketchSize = genomeSize/compress, ATTENTION with MinHash function, for both clust-greedy and clust-mst\n");
	fprintf(stdout, "	-t <int>\t: set the thread number, default 1, for both clust-greedy and clust-mst\n");
	fprintf(stdout, "	-d <double>\t: set the threshold of the clusters from the Minimum Spanning Tree and greedy cluster threshold, default 0.05(0.01) for clust-mst(clust-greedy)\n");
	fprintf(stdout, "	-f\t\t: input files are genomeInfo and MST contents(sketch contents) for clust-mst(clust-greedy)\n");
	fprintf(stdout, "	-E\t\t: input files are genomeInfo and sketch contents contents for clust-mst\n");
	fprintf(stdout, "	-F <string>\t: sketch function, includes MinHash and KSSD default MinHash, for both clust-greedy and clust-mst\n"); 
	fprintf(stdout, "	-o <string>\t: path of output file, for both clust-greedy and clust-mst\n"); 
	fprintf(stdout, "	-i <string>\t: path of input file, ATTENTION with -f and -E option \n"); 

	fprintf(stdout, "example as follows:\n");
	fprintf(stdout, "\tfor genomes clustering, input is a genome file list:\n");
	fprintf(stdout, "\t ./clust-mst -l -t 48 -i bacteriaList -o bacteria.clust\n");
	fprintf(stdout, "\t ./clust-greedy -l -t 48 -i bacteriaList -o bacteria.clust\n");
	fprintf(stdout,"\n");

	fprintf(stdout, "\tfor genomes clustering, input is a single genome file:\n");
	fprintf(stdout, "\t ./clust-mst -d 0.05 -t 48 -i bacteria.fna -o bacteria.clust\n");
	fprintf(stdout, "\t ./clust-greedy -d 0.05 -t 48 -i bacteria.fna -o bacteria.clust\n");
	fprintf(stdout,"\n");

	fprintf(stdout, "\tfor redundancy detection with containment, input is a genome file list:\n");
	fprintf(stdout, "\t ./clust-mst -l -c 10000 -t 48 -i bacteriaList -o bacteria.out\n");
	fprintf(stdout, "\t ./clust-greedy -l -c 10000 -t 48 -i bacteriaList -o bacteria.out\n");
	fprintf(stdout,"\n");

	fprintf(stdout, "\tfor redundancy detection with containment, input is a single genome file:\n");
	fprintf(stdout, "\t ./clust-mst -c 10000 -t 48 -i bacteria.fna -o bacteria.out\n");
	fprintf(stdout, "\t ./clust-greedy -c 10000 -t 48 -i bacteria.fna -o bacteria.out\n");
	fprintf(stdout,"\n");

	fprintf(stdout, "\tfor generator cluster from exist MST:\n");
	fprintf(stdout, "\t ./clust-mst -f -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashMSTInfo -o bacteria.clust\n");
	fprintf(stdout, "\t ATTENTION: the -f must in front of the -i option\n");
	fprintf(stdout,"\n");

	fprintf(stdout, "\tfor generator cluster from exist sketches:\n");
	fprintf(stdout, "\t ./clust-mst -E -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashSketchInfo -o bacteria.clust\n");
	fprintf(stdout, "\t ./clust-greedy -f -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashSketchInfo -o bacteria.clust\n");
	fprintf(stdout, "\t ATTENTION: the -E and -f must in front of the -i option\n");
	fprintf(stdout,"\n");

}




#endif
