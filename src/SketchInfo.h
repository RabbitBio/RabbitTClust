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
	SequenceInfo seqInfo;//for sketch squence;
	bool isContainment = false;
	
	Sketch::MinHash* minHash;
	Sketch::KSSD* KSSD;
	Sketch::WMinHash* WMinHash;
	Sketch::HyperLogLog* HLL;
	Sketch::OrderMinHash * OMH;
};


bool cmpGenomeSize(SketchInfo s1, SketchInfo s2);
bool cmpSeqSize(SketchInfo s1, SketchInfo s2);

void calSize(bool sketchByFile, string inputFile, int threads, uint64_t minLen, uint64_t &maxSize, uint64_t& minSize, uint64_t& averageSize);
bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool cmpSketch(SketchInfo s1, SketchInfo s2);
























#endif //H_SKETCH_INFO
