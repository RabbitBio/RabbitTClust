#ifndef H_SKETCH_IO
#define H_SKETCH_IO
#include "SketchInfo.h"
#include "common.hpp"

void read_sketch_parameters(string folder_path, int& sketch_func_id, int& kmer_size, bool& is_containment, int& contain_compress, int& sketch_size, int& half_k, int& half_subk, int& drlevel);
void save_genome_info(vector<SketchInfo>& sketches, string folderPath, string type, bool sketchByFile);
void saveSketches(vector<SketchInfo>& sketches, string folderPath, bool sketchByFile, string sketchFunc, bool isContainment, int containCompress, int sketchSize, int kmerSize);

bool loadSketches(string folderPath, int threads, vector<SketchInfo>& sketches, int& sketch_func_id);
bool load_genome_info(string folderPath, string type, vector<SketchInfo>& sketches);
#endif
