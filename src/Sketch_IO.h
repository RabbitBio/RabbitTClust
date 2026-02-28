#ifndef H_SKETCH_IO
#define H_SKETCH_IO
#include "SketchInfo.h"
#include "common.hpp"

void read_sketch_parameters(string folder_path, int& sketch_func_id, int& kmer_size, bool& is_containment, int& contain_compress, int& sketch_size, int& half_k, int& half_subk, int& drlevel);
void save_genome_info(vector<SketchInfo>& sketches, string folderPath, string type, bool sketchByFile);
void save_kssd_genome_info(const vector<KssdSketchInfo>& sketches, const string folderPath, const string type, bool sketchByFile);
void saveSketches(vector<SketchInfo>& sketches, string folderPath, bool sketchByFile, string sketchFunc, bool isContainment, int containCompress, int sketchSize, int kmerSize);
void saveKssdSketches(const vector<KssdSketchInfo>& sketches, const KssdParameters info, const string folderPath, bool sketchByFile);

bool loadSketches(string folderPath, int threads, vector<SketchInfo>& sketches, int& sketch_func_id);
bool load_genome_info(string folderPath, string type, vector<SketchInfo>& sketches);
bool load_kssd_genome_info(string folderPath, string type, vector<KssdSketchInfo>& sketches);
bool loadKssdSketches(string folderPath, int threads, vector<KssdSketchInfo>& sketches, KssdParameters& info);

#ifdef USE_MPI
void append_binary_genome_info(FILE* fp_info, const vector<KssdSketchInfo>& sketches, bool sketchByFile);
void append_binary_hash_data(FILE* fp_hash, const vector<KssdSketchInfo>& sketches);
void append_minhash_genome_info(FILE* fp, const vector<SketchInfo>& sketches, bool sketchByFile);
void append_minhash_hash_data(FILE* fp, const vector<SketchInfo>& sketches);
#endif
#endif
