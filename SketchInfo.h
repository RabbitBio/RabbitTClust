#ifndef H_SKETCH_INFO
#define H_SKETCH_INFO

#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>

struct SketchInfo{
	int id;
	std::string name;
	std::string comment;
	std::string strand;//for fastq files.
	uint64_t length;
	int size;
	std::vector<uint64_t> hashes;//the "uint64_t" need to be a parameter for the sketch.

};

struct SimilarityInfo{
	//for sequence information
	int id;
	std::string name;
	std::string comment;
	std::string strand;
	std::string seq;
	uint64_t length;
	int size;


	//for similarity information
	//vector<int> suffixIds;
	//vector<double> similaritys;
	
};

























#endif //H_SKETCH_INFO
