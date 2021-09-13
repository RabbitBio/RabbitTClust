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
	
	Sketch::MinHash* minHash;
	Sketch::WMinHash* WMinHash;
	Sketch::HyperLogLog* HLL;
	Sketch::OrderMinHash * OMH;
};



bool sketchSequences(string inputFile, string sketchFunc, vector<SketchInfo>& sketches, int threads);
bool sketchFiles(string inputFile, string sketchFunc, vector<SketchInfo>& sketches, int threads);
bool cmpSketch(SketchInfo s1, SketchInfo s2);
























#endif //H_SKETCH_INFO
