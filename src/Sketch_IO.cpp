#include "Sketch_IO.h"
#include <fstream>
#include <sstream>
#include "greedy.h"
#include "MST_IO.h"
#include <assert.h>

bool cmpIndex(SketchInfo s1, SketchInfo s2){
	return s1.id < s2.id;
}

void read_sketch_parameters(string folder_path, int& sketch_func_id, int& kmer_size, bool& is_containment, int& contain_compress, int& sketch_size, int& half_k, int& half_subk, int& drlevel){
	string hash_file = folder_path + '/' + "hash.sketch";
	FILE * fp_hash = fopen(hash_file.c_str(), "r");
	if(!fp_hash){
		cerr << "ERROR: read_sketch_parameters(), cannot open the file: " << hash_file << endl;
		exit(1);
	}
	fread(&sketch_func_id, sizeof(int), 1, fp_hash);
	if(sketch_func_id == 0){
		fread(&kmer_size, sizeof(int), 1, fp_hash);
		fread(&is_containment, sizeof(bool), 1, fp_hash);
		if(is_containment)
			fread(&contain_compress, sizeof(int), 1, fp_hash);
		else
			fread(&sketch_size, sizeof(int), 1, fp_hash);
	}
	else if(sketch_func_id == 1){
		fread(&half_k, sizeof(int), 1, fp_hash);
		fread(&half_subk, sizeof(int), 1, fp_hash);
		fread(&drlevel, sizeof(int), 1, fp_hash);
	}
	fclose(fp_hash); 
}

void save_genome_info(vector<SketchInfo>& sketches, string folderPath, string type, bool sketchByFile){
	assert(type == "sketch" || type == "mst");
	string info_file = folderPath + '/' + "info." + type;
	FILE * fp_info = fopen(info_file.c_str(), "w+");
	if(!fp_info){
		cerr << "ERROR: save_genome_info(), cannot open the file: " << info_file << endl;
		exit(1);
	}
	fwrite(&sketchByFile, sizeof(bool), 1, fp_info);
	size_t sketch_number = sketches.size();
	fwrite(&sketch_number, sizeof(size_t), 1, fp_info);

	if(sketchByFile)
	{
		for(int i = 0; i < sketches.size(); i++)
		{
			Vec_SeqInfo curFileSeqs = sketches[i].fileSeqs;
			int file_name_length = sketches[i].fileName.length();
			int seq0_name_length = curFileSeqs[0].name.length();
			int seq0_comment_length = curFileSeqs[0].comment.length();
			fwrite(&file_name_length, sizeof(int), 1, fp_info);
			fwrite(&seq0_name_length, sizeof(int), 1, fp_info);
			fwrite(&seq0_comment_length, sizeof(int), 1, fp_info);
			fwrite(&curFileSeqs[0].strand, sizeof(int), 1, fp_info);
			fwrite(&sketches[i].totalSeqLength, sizeof(uint64_t), 1, fp_info);
			fwrite(sketches[i].fileName.c_str(), sizeof(char), file_name_length, fp_info);
			fwrite(curFileSeqs[0].name.c_str(), sizeof(char), seq0_name_length, fp_info);
			fwrite(curFileSeqs[0].comment.c_str(), sizeof(char), seq0_comment_length, fp_info);
		}
	}
	else//sketchBySequence
	{
		for(int i = 0; i < sketches.size(); i++)
		{
			SequenceInfo curSeq = sketches[i].seqInfo;
			int seq_name_length = curSeq.name.length();
			int seq_comment_length = curSeq.comment.length();
			fwrite(&seq_name_length, sizeof(int), 1, fp_info);
			fwrite(&seq_comment_length, sizeof(int), 1, fp_info);
			fwrite(&curSeq.strand, sizeof(int), 1, fp_info);
			fwrite(&curSeq.length, sizeof(int), 1, fp_info);
			fwrite(curSeq.name.c_str(), sizeof(char), seq_name_length, fp_info);
			fwrite(curSeq.comment.c_str(), sizeof(char), seq_comment_length, fp_info);
		}
	}
	fclose(fp_info);
}


void saveSketches(vector<SketchInfo>& sketches, string folderPath, bool sketchByFile, string sketchFunc, bool isContainment, int containCompress, int sketchSize, int kmerSize)
{
	//-----save the info.sketch
	save_genome_info(sketches, folderPath, "sketch", sketchByFile);

	//-----save the hash.sketch
	string hash_file = folderPath + '/' + "hash.sketch";
	FILE * fp_hash = fopen(hash_file.c_str(), "w+");
	if(!fp_hash){
		cerr << "ERROR: saveSketch(), cannot open the file: " << hash_file << endl;
		exit(1);
	}
	int sketch_func_id = 0;
	if(sketchFunc == "MinHash"){
		sketch_func_id = 0;
	}
	else if(sketchFunc == "KSSD"){
		sketch_func_id = 1;
	}
	else{
		cerr << "ERROR: saveSketches(), noexistent hash function: " << sketchFunc << endl;
		exit(1);
	}
	fwrite(&sketch_func_id, sizeof(int), 1, fp_hash);
	if(sketch_func_id == 0){
		fwrite(&kmerSize, sizeof(int), 1, fp_hash);
		fwrite(&isContainment, sizeof(bool), 1, fp_hash);
		if(isContainment)
			fwrite(&containCompress, sizeof(int), 1, fp_hash);
		else
			fwrite(&sketchSize, sizeof(int), 1, fp_hash);
		//-----saving hash number and hash values for each sketch
		for(int i = 0; i < sketches.size(); i++){
			vector<uint64_t> hashArr = sketches[i].minHash->storeMinHashes();
			size_t cur_sketch_size = hashArr.size();
			fwrite(&cur_sketch_size, sizeof(size_t), 1, fp_hash);
			fwrite(hashArr.data(), sizeof(uint64_t), cur_sketch_size, fp_hash);
		}
	}
	else if(sketch_func_id == 1){
		int half_k, half_subk, drlevel;
		half_k = sketches[0].KSSD->get_half_k();
		half_subk = sketches[0].KSSD->get_half_subk();
		drlevel = sketches[0].KSSD->get_drlevel();
		fwrite(&half_k, sizeof(int), 1, fp_hash);
		fwrite(&half_subk, sizeof(int), 1, fp_hash);
		fwrite(&drlevel, sizeof(int), 1, fp_hash);
		for(int i = 0; i < sketches.size(); i++){
			vector<uint64_t> hashArr = sketches[i].KSSD->storeHashes();
			size_t cur_sketch_size = hashArr.size();
			fwrite(&cur_sketch_size, sizeof(size_t), 1, fp_hash);
			fwrite(hashArr.data(), sizeof(uint64_t), cur_sketch_size, fp_hash);
		}
	}
	fclose(fp_hash); 
	cerr << "-----save the sketches into: " << folderPath << endl;
	
}

bool loadSketches(string folderPath, int threads, vector<SketchInfo>& sketches, int& sketch_func_id){
	string hash_file = folderPath + '/' + "hash.sketch";
	FILE * fp_hash = fopen(hash_file.c_str(), "r");
	if(!fp_hash){
		cerr << "ERROR: loadSketches(), cannot open the file: " << hash_file << endl;
		exit(1);
	}
	fread(&sketch_func_id, sizeof(int), 1, fp_hash);
	int kmer_size, contain_compress, sketch_size;
	int half_k(10), half_subk(6), drlevel(3);
	bool is_containment;
	if(sketch_func_id == 0){ //MinHash
		fread(&kmer_size, sizeof(int), 1, fp_hash);
		fread(&is_containment, sizeof(bool), 1, fp_hash);
		if(is_containment)
			fread(&contain_compress, sizeof(int), 1, fp_hash);
		else
			fread(&sketch_size, sizeof(int), 1, fp_hash);
	}
	else if(sketch_func_id == 1){//KSSD
		fread(&half_k, sizeof(int), 1, fp_hash);
		fread(&half_subk, sizeof(int), 1, fp_hash);
		fread(&drlevel, sizeof(int), 1, fp_hash);
	}

	Sketch::KSSDParameters kssdPara(half_k, half_subk, drlevel);

	bool sketch_by_file = load_genome_info(folderPath, "sketch", sketches);
	int max_hash_number = 1 << 20;
	uint64_t * buffer_hash_arr = new uint64_t [max_hash_number];
	for(size_t i = 0; i < sketches.size(); i++){
		Sketch::MinHash * mh1;
		Sketch::KSSD * kssd;
		size_t cur_sketch_size;
		fread(&cur_sketch_size, sizeof(size_t), 1, fp_hash);
		if(cur_sketch_size > max_hash_number){
			max_hash_number = cur_sketch_size;
			buffer_hash_arr = new uint64_t [max_hash_number];
		}
		int cur_hash_number = fread(buffer_hash_arr, sizeof(uint64_t), cur_sketch_size, fp_hash);
		assert(cur_hash_number == cur_sketch_size);
		vector<uint64_t> hash_arr(buffer_hash_arr, buffer_hash_arr+cur_hash_number);
		if(sketch_func_id == 0){//MinHash
			if(is_containment){
				mh1 = new Sketch::MinHash(kmer_size, contain_compress);
				sketches[i].isContainment = true;
			}
			else{
				mh1 = new Sketch::MinHash(kmer_size, sketch_size);
			}
			mh1->loadMinHashes(hash_arr);
			sketches[i].minHash = mh1;
		}
		else if(sketch_func_id == 1){
			kssd = new Sketch::KSSD(kssdPara);
			kssd->loadHashes(hash_arr);
			sketches[i].KSSD = kssd;
		}
		sketches[i].id = i;
	}

	//std::sort(sketches.begin(), sketches.end(), cmpIndex);
	return sketch_by_file;
}

bool load_genome_info(string folderPath, string type, vector<SketchInfo>& sketches){
	assert(type == "sketch" || type == "mst");
	string file_genome_info = folderPath + '/' + "info." + type;
	FILE* fp_info = fopen(file_genome_info.c_str(), "r");
	if(!fp_info){
		cerr << "ERROR: loadMST(), cannot open file: " << file_genome_info << endl;
		exit(1);
	}
	bool sketch_by_file;
	fread(&sketch_by_file, sizeof(bool), 1, fp_info);
	size_t sketch_number;
	fread(&sketch_number, sizeof(size_t), 1, fp_info);
	int max_file_name_length = 1 << 12;
	int max_seq_name_length = 1 << 12;
	int max_seq_comment_length = 1 << 16;
	char * buffer_file_name = new char[max_file_name_length+1];
	char * buffer_seq_name = new char[max_seq_name_length+1];
	char * buffer_seq_comment = new char[max_seq_comment_length+1];
	if(sketch_by_file){
		int file_name_length, seq0_name_length, seq0_comment_length, strand;
		uint64_t total_seq_length;
		for(int i = 0; i < sketch_number; i++){
			SketchInfo cur_sketch_info;
			Sketch::MinHash * mh1;
			Sketch::KSSD * kssd;
			fread(&file_name_length, sizeof(int), 1, fp_info);
			fread(&seq0_name_length, sizeof(int), 1, fp_info);
			fread(&seq0_comment_length, sizeof(int), 1, fp_info);
			fread(&strand, sizeof(int), 1, fp_info);
			fread(&total_seq_length, sizeof(uint64_t), 1, fp_info);
			if(file_name_length > max_file_name_length){
				max_file_name_length = file_name_length;
				buffer_file_name = new char [max_file_name_length+1];
			}
			if(seq0_name_length > max_seq_name_length){
				max_seq_name_length = seq0_name_length;
				buffer_seq_name = new char [max_seq_name_length+1];
			}
			if(seq0_comment_length > max_seq_comment_length){
				max_seq_comment_length = seq0_comment_length;
				buffer_seq_comment = new char [max_seq_comment_length+1];
			}
			int cur_name_length = fread(buffer_file_name, sizeof(char), file_name_length, fp_info);
			assert(cur_name_length == file_name_length);
			int cur_seq_name_length = fread(buffer_seq_name, sizeof(char), seq0_name_length, fp_info);
			assert(cur_seq_name_length = seq0_name_length);
			int cur_seq_comment_length = fread(buffer_seq_comment, sizeof(char), seq0_comment_length, fp_info);
			assert(cur_seq_comment_length == seq0_comment_length);
			string file_name, seq0_name, seq0_comment;
			file_name.assign(buffer_file_name, buffer_file_name + cur_name_length);
			seq0_name.assign(buffer_seq_name, buffer_seq_name + cur_seq_name_length);
			seq0_comment.assign(buffer_seq_comment, buffer_seq_comment + cur_seq_comment_length);
			SequenceInfo tmp_seq{seq0_name, seq0_comment, strand, 0};
			Vec_SeqInfo cur_file_seqs;
			cur_file_seqs.push_back(tmp_seq);
			cur_sketch_info.fileName = file_name;
			cur_sketch_info.totalSeqLength = total_seq_length;
			cur_sketch_info.fileSeqs = cur_file_seqs;
			cur_sketch_info.id = i;
			sketches.push_back(cur_sketch_info);
		}
	}
	else{
		int seq_name_length, seq_comment_length, strand, length;
		for(int i = 0; i < sketch_number; i++){
			SketchInfo cur_sketch_info;
			Sketch::MinHash * mh1;
			Sketch::KSSD * kssd;
			fread(&seq_name_length, sizeof(int), 1, fp_info);
			fread(&seq_comment_length, sizeof(int), 1, fp_info);
			fread(&strand, sizeof(int), 1, fp_info);
			fread(&length, sizeof(int), 1, fp_info);
			if(seq_name_length > max_seq_name_length){
				max_seq_name_length = seq_name_length;
				buffer_seq_name = new char [max_seq_name_length+1];
			}
			if(seq_comment_length > max_seq_comment_length){
				max_seq_comment_length = seq_comment_length;
				buffer_seq_comment = new char [max_seq_comment_length+1];
			}
			int cur_seq_name_length = fread(buffer_seq_name, sizeof(char), seq_name_length, fp_info);
			assert(cur_seq_name_length == seq_name_length);
			int cur_seq_comment_length = fread(buffer_seq_comment, sizeof(char), seq_comment_length, fp_info);
			assert(cur_seq_comment_length == seq_comment_length);
			string seq_name, seq_comment;
			seq_name.assign(buffer_seq_name, buffer_seq_name + cur_seq_name_length);
			seq_comment.assign(buffer_seq_comment, buffer_seq_comment + cur_seq_comment_length);
			SequenceInfo cur_seq{seq_name, seq_comment, strand, length};
			cur_sketch_info.seqInfo = cur_seq;
			cur_sketch_info.id = i;
			sketches.push_back(cur_sketch_info);
		}
	}
	fclose(fp_info);
	return sketch_by_file;
}





