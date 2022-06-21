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
void loadDense(int** &denseArr, string inputFile, int denseSpan, vector<SketchInfo> sketches);

bool loadMSTs(string inputFile, string inputFile1, vector<SketchInfo>& sketches, vector<EdgeInfo>& mst);

void printResult(std::vector< std::vector<int> > clusterOrigin, std::vector<SketchInfo> sketches, bool sketchByFile, string outputFile);

void saveMST(string folderPath, string inputFile, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo> sketches, vector<EdgeInfo> mst, bool sketchByFile, int sketchSize, int kmerSize);

void saveDense(string folderPath, string prefixName, int** denseArr, int denseSpan, vector<SketchInfo> sketches);

void saveANI(string folderPath, string prefixName, uint64_t* aniArr, string sketchFunc);

#endif

