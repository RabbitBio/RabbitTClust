#ifndef H_SKETCH_INFO
#define H_SKETCH_INFO

#include "Sketch.h"
#include "phmap.h"
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

// Inverted index structure for incremental building during sketch generation
struct KssdInvertedIndex{
	bool use64;
	// For use64 case: hash_map from hash value to list of genome IDs
	phmap::flat_hash_map<uint64_t, vector<uint32_t>> hash_map_64;
	// For non-use64 case: vector of vectors indexed by hash value
	vector<vector<uint32_t>> hash_map_32;
	size_t hashSize;  // Only used for non-use64 case
	
	KssdInvertedIndex() : use64(false), hashSize(0) {}
	
	void init(bool use64_val, size_t hashSize_val = 0) {
		use64 = use64_val;
		if(use64) {
			hash_map_64.reserve(1000000);  // Pre-allocate
		} else {
			hashSize = hashSize_val;
			hash_map_32.resize(hashSize);
		}
	}
	
	// Thread-safe insert with critical section
	void insert_hash_safe(uint64_t hash_val, uint32_t genome_id) {
		if(use64) {
			#pragma omp critical
			{
				hash_map_64[hash_val].push_back(genome_id);
			}
		} else {
			uint32_t hash = (uint32_t)hash_val;
			if(hash < hash_map_32.size()) {
				#pragma omp critical
				{
					hash_map_32[hash].push_back(genome_id);
				}
			}
		}
	}
};

bool cmpGenomeSize(SketchInfo s1, SketchInfo s2);
bool cmpSeqSize(SketchInfo s1, SketchInfo s2);

void calSize(bool sketchByFile, string inputFile, int threads, uint64_t minLen, uint64_t &maxSize, uint64_t& minSize, uint64_t& averageSize);
bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads);
bool cmpSketch(SketchInfo s1, SketchInfo s2);
//bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, int threads);
bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index = nullptr);
bool sketchSequencesWithKssd(const string inputFile, const int minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index = nullptr);
void transSketches(const vector<KssdSketchInfo>& sketches, const KssdParameters& info, const string folder_path, int numThreads);
void transSketchesFromIndex(const KssdInvertedIndex& inverted_index, const KssdParameters& info, const string folder_path);



#endif //H_SKETCH_INFO
