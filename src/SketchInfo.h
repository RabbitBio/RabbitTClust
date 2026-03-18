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
	// For non-use64 case: phmap (sparse, same as use64) to save memory
	phmap::flat_hash_map<uint32_t, vector<uint32_t>> hash_map_32;
	
	KssdInvertedIndex() : use64(false) {}
	
	void init(bool use64_val, size_t /* hashSize_val unused for phmap */ = 0) {
		use64 = use64_val;
		if(use64) {
			hash_map_64.reserve(1000000);  // Pre-allocate
		} else {
			hash_map_32.reserve(1000000);  // Pre-allocate, sparse
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
			#pragma omp critical
			{
				hash_map_32[hash].push_back(genome_id);
			}
		}
	}
};

// Inverted index structure for MinHash (similar to KssdInvertedIndex)
struct MinHashInvertedIndex{
	// Hash map from hash value to list of genome IDs
	phmap::flat_hash_map<uint64_t, vector<uint32_t>> hash_map;
	
	MinHashInvertedIndex() {
		hash_map.reserve(1000000);  // Pre-allocate
	}
	
	// Thread-safe insert with critical section
	void insert_hash_safe(uint64_t hash_val, uint32_t genome_id) {
		#pragma omp critical
		{
			hash_map[hash_val].push_back(genome_id);
		}
	}

	static string index_path(const string& folder) {
		return folder + "/minhash.sketch.index";
	}

	bool save_to_file(const string& folder) const {
		string path = index_path(folder);
		FILE* fp = fopen(path.c_str(), "wb");
		if (!fp) { cerr << "WARNING: cannot save MinHash index to: " << path << endl; return false; }
		const char magic[9] = "MHIDX001";
		fwrite(magic, 1, 8, fp);
		size_t n = hash_map.size();
		fwrite(&n, sizeof(size_t), 1, fp);
		for (const auto& kv : hash_map) {
			uint64_t h = kv.first;
			uint32_t m = (uint32_t)kv.second.size();
			fwrite(&h, sizeof(uint64_t), 1, fp);
			fwrite(&m, sizeof(uint32_t), 1, fp);
			fwrite(kv.second.data(), sizeof(uint32_t), m, fp);
		}
		fclose(fp);
		cerr << "-----MinHash inverted index saved: " << path
		     << " (" << n << " unique hashes)" << endl;
		return true;
	}

	bool load_from_file(const string& folder) {
		string path = index_path(folder);
		FILE* fp = fopen(path.c_str(), "rb");
		if (!fp) return false;
		char magic[9] = {};
		if (fread(magic, 1, 8, fp) != 8 || string(magic, 8) != "MHIDX001") {
			fclose(fp); return false;
		}
		size_t n;
		fread(&n, sizeof(size_t), 1, fp);
		hash_map.clear();
		hash_map.reserve(n);
		for (size_t i = 0; i < n; i++) {
			uint64_t h; uint32_t m;
			fread(&h, sizeof(uint64_t), 1, fp);
			fread(&m, sizeof(uint32_t), 1, fp);
			auto& vec = hash_map[h];
			vec.resize(m);
			fread(vec.data(), sizeof(uint32_t), m, fp);
		}
		fclose(fp);
		cerr << "-----MinHash inverted index loaded: " << path
		     << " (" << n << " unique hashes)" << endl;
		return true;
	}
};

bool cmpGenomeSize(SketchInfo s1, SketchInfo s2);
bool cmpSeqSize(SketchInfo s1, SketchInfo s2);

void calSize(bool sketchByFile, string inputFile, int threads, uint64_t minLen, uint64_t &maxSize, uint64_t& minSize, uint64_t& averageSize);
bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads, MinHashInvertedIndex* inverted_index = nullptr);
bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads, MinHashInvertedIndex* inverted_index = nullptr);
bool cmpSketch(SketchInfo s1, SketchInfo s2);
//bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, int threads);
bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index = nullptr);
bool sketchSequencesWithKssd(const string inputFile, const int minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index = nullptr);
void transSketches(const vector<KssdSketchInfo>& sketches, const KssdParameters& info, const string folder_path, int numThreads);
void transSketchesFromIndex(const KssdInvertedIndex& inverted_index, const KssdParameters& info, const string folder_path);
void saveMinHashIndex(const MinHashInvertedIndex& inverted_index, const string folder_path);
bool loadMinHashIndex(const string folder_path, MinHashInvertedIndex& inverted_index);



#endif //H_SKETCH_INFO
