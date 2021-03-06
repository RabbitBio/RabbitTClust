#ifndef H_MST_IO
#define H_MST_IO

#include <fstream>
#include <sstream>
#include "SketchInfo.h"//include <iostream> <string>
#include "MST.h" //include <vector>
#include "parameter.h"

struct ClusterInfo{
	int id;
	uint64_t length;
};

void MST2Cluster(std::string inputFile, std::string inputFile1, std::string outputFile, double threshold);

void printResult(std::vector< std::vector<int> > clusterOrigin, std::vector<SketchInfo> sketches, bool sketchByFile, string outputFile);

void saveMST(string inputFile, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo> sketches, vector<EdgeInfo> mst, bool sketchByFile, int sketchSize, int kmerSize);


#endif

