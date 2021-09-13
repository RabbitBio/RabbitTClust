#include "parameter.h"

void printUsage(){
	fprintf(stdout, "usage: clust [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> genomeFile\n");
	fprintf(stdout, "usage: clust [-h] [-f] [-d] <double> genomeInfo mstInfo\n");
	fprintf(stdout, "	-h		: this help message\n");
	fprintf(stdout, "	-l		: genome clustering, inputFile is the path list of the genome files\n");
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
