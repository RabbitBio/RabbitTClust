#ifndef H_SKETCH_INFO
#define H_SKETCH_INFO

#include "Sketch.h"
#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>
#include <random>

using namespace std;

//struct SketchInfo{
//	int id;
//	std::string name;
//	std::string comment;
//	std::string strand;//for fastq files.
//	uint64_t length;
//	int size;
//	std::vector<uint64_t> hashes;//the "uint64_t" need to be a parameter for the sketch.
//
//};

struct SketchInfo{
	Sketch::MinHash* minHash;
	Sketch::WMinHash* WMinHash;
	Sketch::HyperLogLog* HLL;
	Sketch::OrderMinHash * OMH;
	int index;
};

struct SimilarityInfo{
	//for sequence information
	int id;
	string name;
	string comment;
	string strand;
	string seq;
	uint64_t length;
	int size;
	
};


bool sketchSequences(string inputFile, string sketchFunc, vector<SimilarityInfo>& similarityInfos, vector<SketchInfo>& sketches, int threads);
bool sketchFiles(string inputFile, string sketchFunc, vector<SimilarityInfo>& similarityInfos, vector<SketchInfo>& sketches, int threads);
























#endif //H_SKETCH_INFO
