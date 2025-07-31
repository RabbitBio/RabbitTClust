#ifndef H_SKETCH_INFO
#define H_SKETCH_INFO

#include "Sketch.h"
#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>
#include <random>

using namespace std;

//for sequence information
struct SequenceInfo{
	string name;
	string comment;
	int strand;
	int length;
};

typedef vector<SequenceInfo> Vec_SeqInfo;
struct SketchInfo{
	int id;
	string fileName;//for sketch files;
	uint64_t totalSeqLength;
	Vec_SeqInfo fileSeqs;//for sketch files;
	SequenceInfo seqInfo;//for sketch sequence;
	bool isContainment = false;
	
	Sketch::MinHash* minHash;
	Sketch::KSSD* KSSD;
	Sketch::WMinHash* WMinHash;
	Sketch::HyperLogLog* HLL;
	Sketch::OrderMinHash * OMH;
};

struct KssdSketchInfo{
	int id;
	string fileName;
	uint64_t totalSeqLength;
	Vec_SeqInfo fileSeqs;
	SequenceInfo seqInfo;
	bool use64;
	vector<uint32_t> hash32_arr;
	vector<uint64_t> hash64_arr;
  uint32_t sketchsize;
};

struct KssdParameters{
	int id;
	int half_k;
	int half_subk;
	int drlevel;
	int genomeNumber;
};


bool cmpGenomeSize(SketchInfo s1, SketchInfo s2);
bool cmpSeqSize(SketchInfo s1, SketchInfo s2);

void calSize(bool sketchByFile, string inputFile, int threads, uint64_t minLen, uint64_t &maxSize, uint64_t& minSize, uint64_t& averageSize);
bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool cmpSketch(SketchInfo s1, SketchInfo s2);
//bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, int threads);
bool sketchFileWithKssd(const vector<string> inputFile, const uint64_t minLen, int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads);
bool sketchSequencesWithKssd(const string inputFile, const int minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads);
void transSketches(const vector<KssdSketchInfo>& sketches, const KssdParameters& info, const string folder_path, int numThreads);

void transSketches_in_memory(
    const std::vector<KssdSketchInfo>& sketches, 
    const KssdParameters& info, 
    int numThreads,
        char*& sum_info_buffer_out, size_t& sum_info_size_out,
            char*& sum_hash_buffer_out, size_t& sum_hash_size_out,  
    char*& sum_index_buffer_out, size_t& sum_index_size_out,
    char*& sum_dict_buffer_out, size_t& sum_dict_size_out
    
    ); 


#endif //H_SKETCH_INFO
