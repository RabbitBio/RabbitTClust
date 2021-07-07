#include "parameter.h"
void printUsage(){
	fprintf(stdout, "usage: clust [-h] [-l] [-t] <int> [-d] <double> -F <string> [-f] genomeFile\n");
	fprintf(stdout, "usage: clust [-h] [-f] [-d] <double> genomeInfo mstInfo\n");
	fprintf(stdout, "	-h		: this help message\n");
	fprintf(stdout, "	-l		: genome clustering, inputFile is the path list of the genome files\n");
	fprintf(stdout, "	-t		: genome clustering, set the thread number\n");
	fprintf(stdout, "	-d		: set the threshold of the clusters from the Minimum Spanning Tree\n");
	fprintf(stdout, "	-f		: input files are genomeInfo and MST content\n");
	fprintf(stdout, "	-F		: sketch function, includes MinHash, WMH, OMH, HLL\n"); 
	
}
