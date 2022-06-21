#ifndef H_SKETCH_IO
#define H_SKETCH_IO
#include "SketchInfo.h"
#include "parameter.h"

void saveSketches(vector<SketchInfo> sketches, string folderPath, string inputFile, string sketchFunc, bool isContainment, int containCompress, bool sketchByFile, int sketchSize, int kmerSize);

bool loadSketches(string inputFile0, string inputFile1, int threads, vector<SketchInfo>& sketches);
#endif
