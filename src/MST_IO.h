#ifndef H_MST_IO
#define H_MST_IO

#include <fstream>
#include <sstream>
#include "SketchInfo.h"//include <iostream> <string>
#include "MST.h" //include <vector>
#include "common.hpp"

struct ClusterInfo{
	int id;
	uint64_t length;
};

void print_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output);
void printResult(std::vector<std::vector<int>>& clusterOrigin, std::vector<SketchInfo>& sketches, bool sketchByFile, string outputFile);

void loadMST(string folderPath, vector<EdgeInfo>& mst);
void loadDense(int** &denseArr, string folderPath, int& denseSpan, int& genome_number);
void loadANI(string folderPath, uint64_t* &aniArr, int sketch_func_id);

void saveMST(vector<SketchInfo>& sketches, vector<EdgeInfo>& mst, string folderPath, bool sketchByFile);
void saveDense(string folderPath, int** denseArr, int denseSpan, int genome_number);
void saveANI(string folderPath, uint64_t* aniArr, int sketch_func_id);

#endif

