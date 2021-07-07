#ifndef H_MST_IO
#define H_MST_IO

#include <fstream>
#include <sstream>
#include "SketchInfo.h"//include <iostream> <string>
#include "MST.h" //include <vector>
#include "parameter.h"


void MST2Cluster(std::string inputFile, std::string inputFile1, double threshold);

void printResult(std::vector< std::vector<int> > cluster, std::vector<SimilarityInfo> similarityInfos, bool sketchByFile);


void saveMST(string inputFile, string sketchFunc, vector<SimilarityInfo> similarityInfos, vector<EdgeInfo> mst, bool sketchByFile);



#endif

