#include "SketchInfo.h"
#include "kseq.h"
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <sys/stat.h>

#ifdef THREADPOOL_MINHASH
#include "ThreadPool.h"
#endif

#ifdef RABBIT_FX
#include "FastxStream.h"
#include "FastxChunk.h"
#include "Formater.h"
#include "DataQueue.h"
#include <thread>
#include <cstdint>
#include <unordered_map>
#endif

#include "common.hpp"

#define _64MASK 0xffffffffffffffffLLU

KSEQ_INIT(gzFile, gzread);
using namespace std;

bool cmpSketch(SketchInfo s1, SketchInfo s2){
	return s1.id < s2.id;
}
bool cmpGenomeSize(SketchInfo s1, SketchInfo s2){
	if(s1.totalSeqLength > s2.totalSeqLength)	return true;
	else if(s1.totalSeqLength == s2.totalSeqLength)	return s1.id < s2.id;
	else	return false;
}


bool KssdcmpGenomeSize(KssdSketchInfo s1, KssdSketchInfo s2){
	if(s1.totalSeqLength > s2.totalSeqLength)	return true;
	else if(s1.totalSeqLength == s2.totalSeqLength)	return s1.id < s2.id;
	else	return false;
}

bool KssdcmpSketchSize(KssdSketchInfo s1, KssdSketchInfo s2){
	if(s1.sketchsize > s2.sketchsize)	return true;
	else if(s1.sketchsize == s2.sketchsize)	return s1.id < s2.id;
	else	return false;
}

bool cmpSeqSize(SketchInfo s1, SketchInfo s2){
	if(s1.seqInfo.length > s2.seqInfo.length)	return true;
	else if(s1.seqInfo.length == s2.seqInfo.length)	return s1.id < s2.id;
	else return false;
}

int * shuffle(int arr[], int length, uint64_t seed)

{
	if(length > RAND_MAX){
		fprintf(stderr, "shuffling array length %d must be less than RAND_MAX: %d", length, RAND_MAX);
		exit(1);
	}
	srand(seed);
	int j, tmp;
	for(int i = length-1; i > 0; i--)
	{
		j = rand() % (i + 1);
		tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}
	
	return arr;
}

int * shuffleN(int n, int base)
{
	int * arr;
	arr = (int* ) malloc(n * sizeof(int));
	for(int i = 0; i < n; i++){
		arr[i] = i + base;
	}
	
	return shuffle(arr, n, 23);
}

int * generate_shuffle_dim(int half_subk){
	int dim_size = 1 << 4 * half_subk;
	int * shuffled_dim = shuffleN(dim_size, 0);
	shuffled_dim = shuffle(shuffled_dim, dim_size, 348842630);

	//printf("print the shuffle_dim : \n");
	//for(int i = 0; i < dim_size; i++)
	//	printf("%lx\n", shuffled_dim[i]);
	//exit(0);

	return shuffled_dim;
}

#ifdef THREADPOOL_MINHASH
struct SketchInput{
	SketchInput(SketchInfo sketchInfoNew, string sketchFuncNew, char * seqNew, int indexNew): sketchInfo(sketchInfoNew), sketchFunc(sketchFuncNew), seq(seqNew), index(indexNew){}

	SketchInfo sketchInfo;
	string sketchFunc;
	char * seq;
	int index;

};

struct SketchOutput{
	SketchInfo sketchInfo;
	//int index;
};

SketchOutput* sketchBySequence(SketchInput* input){
	SketchOutput * output = new SketchOutput();
	SketchInfo sketchInfo = input->sketchInfo;
	string sketchFunc = input->sketchFunc;
	if(sketchFunc == "MinHash"){
		sketchInfo.minHash->update(input->seq);
	}
	else if(sketchFunc == "KSSD"){
		sketchInfo.KSSD->update(input->seq);
	}
	else if(sketchFunc == "WMH"){
		//wmh->update(ks1->seq.s);
		//wmh->computeHistoSketch();
		sketchInfo.WMinHash->update(input->seq);
		sketchInfo.WMinHash->computeHistoSketch();
	}
	else if(sketchFunc == "HLL"){
		sketchInfo.HLL->update(input->seq);
	}
	else if(sketchFunc == "OMH"){
		sketchInfo.OMH->buildSketch(input->seq);
	}
	free(input->seq);
	sketchInfo.id = input->index;
	output->sketchInfo = sketchInfo;
	return output;
}


void useThreadOutput(SketchOutput * output, vector<SketchInfo> &sketches){
	//SketchInfo tmpSimilarityInf
	sketches.push_back(output->sketchInfo);
}

#endif

#ifdef RABBIT_FX
typedef rabbit::core::TDataQueue<rabbit::fa::FastaChunk> FaChunkQueue;
int producer_fasta_task(std::string file, rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq){
	//std::cerr << "filename: " << file << std::endl;
	rabbit::fa::FastaFileReader* faFileReader;
	faFileReader = new rabbit::fa::FastaFileReader(file, *fastaPool, false);
	int n_chunks = 0;
	while(1){
		rabbit::fa::FastaChunk* faChunk = new rabbit::fa::FastaChunk;
		faChunk = faFileReader->readNextChunkList();
		if(faChunk == NULL) break;
		n_chunks++;
		dq.Push(n_chunks, faChunk);
		//cerr << "reading : " << n_chunks << endl;

	}
	dq.SetCompleted();

	return 0;
}

void consumer_fasta_seqSize(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, uint64_t minLen, uint64_t* maxSize, uint64_t* minSize, uint64_t* totalSize,	uint64_t* number, uint64_t* badNumber){
	rabbit::int64 id = 0;
	rabbit::fa::FastaChunk *faChunk;
	while(dq.Pop(id, faChunk)){
		std::vector<Reference> data;
		int ref_num = rabbit::fa::chunkListFormat(*faChunk, data);
		for(Reference &r: data){
			uint64_t length = (uint64_t)r.length;
			if(length < minLen){
				*badNumber = *badNumber + 1;
				continue;
			}
			*maxSize = std::max(*maxSize, length);
			*minSize = std::min(*minSize, length);
			*totalSize += length;
			*number = *number + 1;
			//sizeArr->push_back(length);
		}
		rabbit::fa::FastaDataChunk *tmp = faChunk->chunk;
		do{
			if(tmp != NULL){
			fastaPool->Release(tmp);
			tmp = tmp->next;
			}
		}while(tmp!=NULL);
	}
}

void consumer_fasta_task(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, Sketch::KSSDParameters *kssdPara, Sketch::WMHParameters * parameters, vector<SketchInfo> *sketches){
	int line_num = 0;
	rabbit::int64 id = 0;

	rabbit::fa::FastaChunk *faChunk;
	while(dq.Pop(id, faChunk)){
		std::vector<Reference> data;
		int ref_num = rabbit::fa::chunkListFormat(*faChunk, data);
		for(Reference &r: data){
			string name = r.name;
			string comment = r.comment;
			int length = r.length;
			if(length < minLen) continue;
			SequenceInfo curSeq{name, comment, 0, length};

			SketchInfo tmpSketchInfo; 
			tmpSketchInfo.seqInfo = curSeq;
			if(sketchFunc == "MinHash"){
				Sketch::MinHash * mh1;
				if(isContainment)
				{
					int curSketchSize = std::max(length/containCompress, 100);
					mh1 = new Sketch::MinHash(kmerSize, curSketchSize);
					tmpSketchInfo.isContainment = true;
				}
				else
					mh1 = new Sketch::MinHash(kmerSize, sketchSize);
				mh1->update((char*)r.seq.c_str());
				tmpSketchInfo.minHash = mh1;
			}
			else if(sketchFunc == "KSSD"){
				Sketch::KSSD * kssd = new Sketch::KSSD(*kssdPara);
				kssd->update((char*)r.seq.c_str());
				tmpSketchInfo.KSSD = kssd;
			}
			else if(sketchFunc == "WMH"){
				Sketch::WMinHash * wmh = new Sketch::WMinHash(*parameters);
				wmh->update((char*)r.seq.c_str());
				wmh->computeHistoSketch();
				tmpSketchInfo.WMinHash = wmh;
			}
			else if(sketchFunc == "HLL"){
				Sketch::HyperLogLog * hll = new Sketch::HyperLogLog(HLL_SKETCH_BIT);
				hll->update((char*)r.seq.c_str());
				tmpSketchInfo.HLL = hll;
			}
			else if(sketchFunc == "OMH"){
				Sketch::OrderMinHash* omh = new Sketch::OrderMinHash();
				omh->buildSketch((char*)r.seq.c_str());
				tmpSketchInfo.OMH = omh;
			}

			//tmpSketchInfo.minHash = mh1;
			tmpSketchInfo.id = r.gid;
			//if(length >= 10000)
			sketches->push_back(tmpSketchInfo);

		}
		rabbit::fa::FastaDataChunk *tmp = faChunk->chunk;
		do{
			if(tmp != NULL){
			fastaPool->Release(tmp);
			tmp = tmp->next;
			}
		}while(tmp!=NULL);
	}

}

void consumer_fasta_task_with_kssd(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, int minLen, int kmerSize, int drlevel, robin_hood::unordered_map<uint32_t, int> *shuffled_map, vector<KssdSketchInfo> *sketches){
	static const int BaseMap[128] = 
	{
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	// generate the shuffle file.
	int half_k = (kmerSize + 1) / 2;
	bool use64 = half_k - drlevel > 8 ? true : false;
	int half_subk = 6 - drlevel >= 2 ? 6 : drlevel + 2;
	int dim_size = 1 << 4 * half_subk;
	int dim_start = 0;
	int dim_end = 1 << 4 * (half_subk - drlevel);

	int comp_bittl = 64 - 4 * half_k;
	int half_outctx_len = half_k - half_subk;
	int rev_add_move = 4 * half_k - 2;
	//cout << "the comp_bittl is: " << comp_bittl << endl;

	uint64_t tupmask = _64MASK >> comp_bittl;
	uint64_t domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
	uint64_t undomask = (tupmask ^ domask) & tupmask;
	uint64_t undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
	uint64_t undomask0 = undomask ^ undomask1;


	rabbit::int64 id = 0;
	rabbit::fa::FastaChunk *faChunk;
	uint64_t totalLength = 0;
	while(dq.Pop(id, faChunk)){
		std::vector<Reference> data;
		int ref_num = rabbit::fa::chunkListFormat(*faChunk, data);
		for(Reference &r: data){
			string name = r.name;
			string comment = r.comment;
			int length = r.length;
			if(length < minLen) continue;
			SequenceInfo curSeq{name, comment, 0, length};

			//SketchInfo tmpSketchInfo; 
			KssdSketchInfo tmpKssdSketchInfo;
			tmpKssdSketchInfo.seqInfo = curSeq;

			// parse the sequences 
			unordered_set<uint32_t> hashValueSet;
			unordered_set<uint64_t> hashValueSet64;

			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int base = 1;
			//slide window to generate k-mers and get hashes
			for(int i = 0; i < length; i++)
			{
				//char ch = sequence[i];
				char ch = r.seq[i];
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;

				}
				if(base > kmerSize)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);

					//pfilter = shuffled_dim[dim_id];
					//if((pfilter >= dim_end) || (pfilter < dim_start)){
					//	continue;
					//}

					if(shuffled_map->count(dim_id) == 0){
						continue;
					}
					pfilter = (*shuffled_map)[dim_id];

					pfilter -= dim_start;
					//dr_tuple = (((uni_tuple & undomask0) + ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) + pfilter; 
					////only when the dim_end is 4096(the pfilter is 12bit in binary and the lowerst 12bit is all 0 
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmerSize * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 

					if(use64)
						hashValueSet64.insert(dr_tuple);
					else
						hashValueSet.insert(dr_tuple);
				}//end if, i > kmer_size
			}//end for, of a sequence 

			vector<uint32_t> hashArr32;
			vector<uint64_t> hashArr64;
			if(use64){
				for(auto x : hashValueSet64){
					hashArr64.push_back(x);
				}
			}
			else{
				for(auto x : hashValueSet){
					hashArr32.push_back(x);
				}
			}

			tmpKssdSketchInfo.id = r.gid;
			tmpKssdSketchInfo.totalSeqLength = length;
			tmpKssdSketchInfo.use64 = use64;
			tmpKssdSketchInfo.hash32_arr = hashArr32;
			tmpKssdSketchInfo.hash64_arr = hashArr64;
			tmpKssdSketchInfo.sketchsize = hashArr32.size();
      sketches->push_back(tmpKssdSketchInfo);

		}
		rabbit::fa::FastaDataChunk *tmp = faChunk->chunk;
		do{
			if(tmp != NULL){
			fastaPool->Release(tmp);
			tmp = tmp->next;
			}
		}while(tmp!=NULL);
	}//end while loop for this file
}

#endif


void calSize(bool sketchByFile, string inputFile, int threads, uint64_t minLen, uint64_t &maxSize, uint64_t& minSize, uint64_t& averageSize){
	maxSize = 0;
	minSize = 1 << 31;
	averageSize = 0;
	uint64_t totalSize = 0;
	int number = 0;
	int badNumber = 0;
	if(sketchByFile){
		ifstream ifs(inputFile);
		if(!ifs){
			cerr << "ERROR: calSize(), cannot open the inputFile: " << inputFile << endl;
			exit(1);
		}
		string line;
		while(getline(ifs, line)){
			uint64_t curSize;
			string fileSuffix = line.substr(line.length()-2);
			if(fileSuffix == "gz")//gz file
			{ 
				FILE *fp = fopen(line.c_str(), "r");
				fseek(fp, -4, SEEK_END);
				int nUnCompress = 0;
				fread(&nUnCompress, sizeof(int), 1, fp);
				curSize = nUnCompress;
				fclose(fp);
				//cerr << "compressed fileLength is: " << fileLength << endl;
			}
			else
			{
				struct stat statbuf;
				stat(line.c_str(), &statbuf);
				curSize = statbuf.st_size;
			}
			if(curSize < minLen){
				badNumber++;
				continue;
			}
			maxSize = std::max(maxSize, curSize);
			minSize = std::min(minSize, curSize);
			totalSize += curSize;
			number++;
		}
		ifs.close();
	}
	else{//sketch by sequence
		#ifdef RABBIT_FX
		int th = std::max(threads-1, 1);
		uint64_t maxSizeArr[th];
		uint64_t minSizeArr[th];
		uint64_t totalSizeArr[th];
		uint64_t numArr[th];
		uint64_t badNumArr[th];
		for(int i = 0; i < th; i++){
			maxSizeArr[i] = 1;
			minSizeArr[i] = 1 << 31;
			totalSizeArr[i] = 0;
			numArr[i] = 0;
			badNumArr[i] = 0;
		}

		rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1<< 24);
		FaChunkQueue queue1(128, 1);
		std::thread producer(producer_fasta_task, inputFile, fastaPool, std::ref(queue1));
		std::thread **threadArr = new std::thread* [th];
		for(int t = 0; t < th; t++){
			threadArr[t] = new std::thread(std::bind(consumer_fasta_seqSize, fastaPool, std::ref(queue1), minLen, &maxSizeArr[t], &minSizeArr[t], &totalSizeArr[t], &numArr[t], &badNumArr[t]));
		}
		producer.join();

		for(int t = 0; t < th; t++){
			threadArr[t]->join();
		}
		for(int i = 0; i < th; i++){
			maxSize = std::max(maxSize, maxSizeArr[i]);
			minSize = std::min(minSize, minSizeArr[i]);
			totalSize += totalSizeArr[i];
			number += numArr[i];
			badNumber += badNumArr[i];
		}
		#else 
			gzFile fp1 = gzopen(inputFile.c_str(), "r");
			if(!fp1){
				fprintf(stderr, "cannot open the genome file, %s\n", inputFile.c_str());
				exit(1);
			}
			kseq_t * ks1 = kseq_init(fp1);
			while(1){
				int length = kseq_read(ks1);
				if(length < 0) break;
				if(length < minLen){
					badNumber++;
					continue;
				}
				maxSize = std::max(maxSize, (uint64_t)length);
				minSize = std::min(minSize, (uint64_t)length);
				totalSize += (uint64_t)length;
				number++;
			}
		#endif
	}
	averageSize = totalSize / number;
	int totalNumber = number + badNumber;
	cerr << "\t===the genome number for clustering is: " << number << endl;
	cerr << "\t===the genome number below the minimum genome length threshold is: " << badNumber << endl;
	cerr << "\t===the total genome number is: " << totalNumber << endl;
	
	if((double)badNumber / totalNumber >= 0.2)
	{
		fprintf(stderr, "Warning: there are %d poor quality (length < %ld) genome assemblies in the total %d genome assemblied.\n", badNumber, minLen, totalNumber);
	}
	cerr << "\t===the totalSize is: " << totalSize << endl;
	cerr << "\t===the maxSize is: " << maxSize << endl;
	cerr << "\t===the minSize is: " << minSize << endl;
	cerr << "\t===the averageSize is: " << averageSize << endl;
}

bool sketchSequencesWithKssd(const string inputFile, const int minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads){
	int sufIndex = inputFile.find_last_of('.');
	string sufName = inputFile.substr(sufIndex+1);
	if(sufName != "fasta" && sufName != "fna" && sufName != "fa")
	{
		cerr << "error input format file: " << inputFile << endl;
		cerr << "Only support FASTA files" << endl;
		exit(0);
	}
#ifdef RABBIT_FX
	int th = std::max(threads - 1, 1);//consumer threads number;
	vector<KssdSketchInfo> sketchesArr[th];

	int half_subk = 6 - drlevel >= 2 ? 6 : drlevel + 2;
	int dim_size = 1 << 4 * half_subk;
	int dim_start = 0;
	int dim_end = 1 << 4 * (half_subk - drlevel);
	int* shuffled_dim = generate_shuffle_dim(half_subk);
	robin_hood::unordered_map<uint32_t, int> shuffled_map;
	for(int t = 0; t < dim_size; t++){
		if(shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start){
			shuffled_map.insert({t, shuffled_dim[t]});
		}
	}
	
	rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1<< 24);
	FaChunkQueue queue1(128, 1);
	//cout << "--------------------" << endl;
	std::thread producer(producer_fasta_task, inputFile, fastaPool, std::ref(queue1));
	std::thread **threadArr = new std::thread* [th];

	for(int t = 0; t < th; t++){
		threadArr[t] = new std::thread(std::bind(consumer_fasta_task_with_kssd, fastaPool, std::ref(queue1), minLen, kmerSize, drlevel, &shuffled_map, &sketchesArr[t]));
	}
	producer.join();
	for(int t = 0; t < th; t++){
		threadArr[t]->join();
	}

	for(int i = 0; i < th; i++){
		for(int j = 0; j < sketchesArr[i].size(); j++){
			sketches.push_back(sketchesArr[i][j]);
			//similarityInfos.push_back(similarityInfosArr[i][j]);
		}
	}
	int half_k = (kmerSize + 1) / 2;
	//int half_subk = 6 - drlevel >= 2 ? 6 : drlevel + 2;
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	info.id = (half_k << 8) + (half_subk << 4) + drlevel;
	info.genomeNumber = sketches.size();
	
	return true;
#else
	cerr << "need the RabbitFX to parse the genome sequences" << endl;
	return false;
#endif
}


bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads){
	//cerr << "input File is: " << inputFile << endl;
	//cerr << "in sketchSequences(), the minLen is: " << minLen << endl;
	int sufIndex = inputFile.find_last_of('.');
	string sufName = inputFile.substr(sufIndex+1);
	if(sufName != "fasta" && sufName != "fna" && sufName != "fa")
	{
		cerr << "error input format file: " << inputFile << endl;
		cerr << "Only support FASTA files" << endl;
		exit(0);
	}
	gzFile fp1;
	kseq_t * ks1;
	fp1 = gzopen(inputFile.c_str(), "r");
	if(fp1 == NULL){
		fprintf(stderr, "ERROR: sketchSequences(), cannot open the genome file, %s\n", inputFile.c_str());
		//printfUsage();
		return false;
	}
	ks1 = kseq_init(fp1);

	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	Sketch::KSSDParameters kssdPara(half_k, half_subk, drlevel);
	
	Sketch::WMHParameters parameters;
	if(sketchFunc == "WMH"){
		parameters.kmerSize = kmerSize;
		parameters.sketchSize = WMH_SKETCH_SIZE;
		parameters.windowSize = WINDOW_SIZE;
		parameters.r = (double *)malloc(parameters.sketchSize * pow(parameters.kmerSize, 4) * sizeof(double));
		parameters.c = (double *)malloc(parameters.sketchSize * pow(parameters.kmerSize, 4) * sizeof(double));
		parameters.b = (double *)malloc(parameters.sketchSize * pow(parameters.kmerSize, 4) * sizeof(double));
		getCWS(parameters.r, parameters.c, parameters.b, parameters.sketchSize, pow(parameters.kmerSize, 4));
	}

	int index = 0;
	#ifdef THREADPOOL_MINHASH
	ThreadPool<SketchInput, SketchOutput> threadPool(0, threads);
	while(1){
		int length = kseq_read(ks1);
		if(length < 0) break;
		if(length < kmerSize) continue;

		string name("noName");
		string comment("noComment");
		if(ks1->name.s != NULL)
			name = ks1->name.s;
		if(ks1->comment.s != NULL)
			comment = ks1->comment.s;
		SequenceInfo curSeq{name, comment, 0, length};
		
		char * seqCopy = (char*) malloc((length+1) * sizeof(char));
		memcpy(seqCopy, ks1->seq.s, length+1);
		SketchInfo sketchInfo;
		sketchInfo.seqInfo = curSeq;
		if(sketchFunc == "MinHash"){
			Sketch::MinHash * mh1;
			if(isContainment)
			{
				int curSketchSize = std::max(length/containCompress, 100);
				mh1 = new Sketch::MinHash(kmerSize, curSketchSize);
				sketchInfo.isContainment = true;
			}
			else
				mh1 = new Sketch::MinHash(kmerSize, sketchSize);
			sketchInfo.minHash = mh1;
		}
		else if(sketchFunc == "KSSD"){
			Sketch::KSSD * kssd = new Sketch::KSSD(kssdPara);
			sketchInfo.KSSD = kssd;
		}
		else if(sketchFunc == "WMH"){
			Sketch::WMinHash * wmh = new Sketch::WMinHash(parameters);
			sketchInfo.WMinHash = wmh;
		}
		else if(sketchFunc == "HLL"){
			Sketch::HyperLogLog *hll = new Sketch::HyperLogLog(HLL_SKETCH_BIT);
			sketchInfo.HLL = hll;
		}
		else if(sketchFunc == "OMH"){
			Sketch::OrderMinHash *omh = new Sketch::OrderMinHash();
			sketchInfo.OMH = omh;
		}

		threadPool.runWhenThreadAvailable(new SketchInput(sketchInfo, sketchFunc, seqCopy, index), sketchBySequence);

		while(threadPool.outputAvailable()){
			useThreadOutput(threadPool.popOutputWhenAvailable(), sketches);
		}
		index++;

	}//end while
	while(threadPool.running()){
		useThreadOutput(threadPool.popOutputWhenAvailable(), sketches);
	}
	#else 
	#ifdef RABBIT_FX
	int th = std::max(threads - 1, 1);//consumer threads number;
	vector<SketchInfo>  sketchesArr[th];
	//vector<SimilarityInfo>  similarityInfosArr[th];
	
	rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1<< 24);
	FaChunkQueue queue1(128, 1);
	//cout << "--------------------" << endl;
	std::thread producer(producer_fasta_task, inputFile, fastaPool, std::ref(queue1));
	std::thread **threadArr = new std::thread* [th];

	for(int t = 0; t < th; t++){
		threadArr[t] = new std::thread(std::bind(consumer_fasta_task, fastaPool, std::ref(queue1), kmerSize, sketchSize, minLen, sketchFunc, isContainment, containCompress, &kssdPara, &parameters, &sketchesArr[t]));
	}
	producer.join();
	for(int t = 0; t < th; t++){
		threadArr[t]->join();
	}

	for(int i = 0; i < th; i++){
		for(int j = 0; j < sketchesArr[i].size(); j++){
			sketches.push_back(sketchesArr[i][j]);
			//similarityInfos.push_back(similarityInfosArr[i][j]);
		}
	}

	//cerr << "-----the size of sketches is: " << sketches.size() << endl;
	//cerr << "the size of similarityInfos is: " << similarityInfos.size() << endl;
//	exit(0);

	#else 
	//for single thread sketch
	while(1){
		int length = kseq_read(ks1);
		if(length < 0) break;
		if(length < kmerSize) continue;

		string name("noName");
		string comment("noComment");
		if(ks1->name.s != NULL)
			name = ks1->name.s;
		if(ks1->comment.s != NULL)
			comment = ks1->comment.s;
		SequenceInfo curSeq{name, comment, 0, length};

		SketchInfo tmpSketchInfo;
		tmpSketchInfo.seqInfo = curSeq;
		
		if(sketchFunc == "MinHash"){
			Sketch::MinHash* mh1;
			if(isContainment)
			{
				int curSketchSize = std::max(length/containCompress, 100);
				mh1 = new Sketch::MinHash(kmerSize, curSketchSize);
				tmpSketchInfo.isContainment = true;
			}
			else
				mh1 = new Sketch::MinHash(kmerSize, sketchSize);
			mh1->update(ks1->seq.s);
			tmpSketchInfo.minHash = mh1;
		}
		else if(sketchFunc == "KSSD"){
			Sketch::KSSD * kssd = new Sketch::KSSD(kssdPara);
			kssd->update(ks1->seq.s);
			tmpSketchInfo.KSSD = kssd;
		}
		else if(sketchFunc == "WMH"){
			Sketch::WMinHash *wmh = new Sketch::WMinHash(parameters);
			wmh->update(ks1->seq.s);
			wmh->computeHistoSketch();
			tmpSketchInfo.WMinHash = wmh;
		}
		else if(sketchFunc == "HLL"){
			Sketch::HyperLogLog* hll = new Sketch::HyperLogLog(HLL_SKETCH_BIT);
			hll->update(ks1->seq.s);
			tmpSketchInfo.HLL = hll;
		}
		else if(sketchFunc == "OMH"){
			Sketch::OrderMinHash* omh = new Sketch::OrderMinHash();
			omh->buildSketch(ks1->seq.s);
			tmpSketchInfo.OMH = omh;
		}

		tmpSketchInfo.id = index;
		sketches.push_back(tmpSketchInfo);
		index++;
	}//end while
	#endif
	#endif
	//cerr << "the number of sequence is: " << index << endl;
	gzclose(fp1);
	kseq_destroy(ks1);

	//sort(sketches.begin(), sketches.end(), cmpSketch);
	sort(sketches.begin(), sketches.end(), cmpSeqSize);

	return true;
}

bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads){
	fprintf(stderr, "-----input fileList, sketch by file\n");
	fstream fs(inputFile);
	if(!fs){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
		return false;
	}
	vector<string> fileList;
	string fileName;
	while(getline(fs, fileName)){
		fileList.push_back(fileName);
	}

	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	Sketch::KSSDParameters kssdPara(half_k, half_subk, drlevel);

	Sketch::WMHParameters parameter;
	if(sketchFunc == "WMH"){
		parameter.kmerSize = kmerSize;
		parameter.sketchSize = WMH_SKETCH_SIZE;
		parameter.windowSize = WINDOW_SIZE;
		parameter.r = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
		parameter.c = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
		parameter.b = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
		getCWS(parameter.r, parameter.c, parameter.b, parameter.sketchSize, pow(parameter.kmerSize, 4));
	}

	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < fileList.size(); i++){
		//cerr << "start the file: " << fileList[i] << endl;
		gzFile fp1;
		kseq_t* ks1;
		fp1 = gzopen(fileList[i].c_str(), "r");
		if(fp1 == NULL){
			fprintf(stderr, "cannot open the genome file: %s\n", fileList[i].c_str());
			exit(1);
			//return false;
		}
		ks1 = kseq_init(fp1);
		uint64_t totalLength = 0;
		int fileLength = 0;

		Sketch::MinHash * mh1;
		Sketch::KSSD * kssd;
		Sketch::WMinHash * wmh1;
		Sketch::HyperLogLog * hll;
		Sketch::OrderMinHash * omh;

		//get the fileSize(for the sketchSize init)
		if(isContainment)
		{
			FILE * fp = fopen(fileList[i].c_str(), "r");
			if(!fp){
				fprintf(stderr, "cannot open the genome file: %s and get the size of genome file\n", fileList[i].c_str());
				exit(1);
			}

			string fileSuffix = fileList[i].substr(fileList[i].length()-2);
			if(fileSuffix == "gz")//gz file
			{
				fseek(fp, -4, SEEK_END);
				int nUnCompress = 0;
				fread(&nUnCompress, sizeof(int), 1, fp);
				fileLength = nUnCompress;
				//cerr << "compressed fileLength is: " << fileLength << endl;
			}
			else//fna file
			{
				fseek(fp, 0, SEEK_END);
				int nUnCompress = ftell(fp);
				fileLength =  nUnCompress;
				//cerr << "uncompressed fileLength is: " << fileLength << endl;
			}
			fclose(fp);
		}

		if(sketchFunc == "MinHash"){
			if(isContainment){
				int curSketchSize = std::max(fileLength / containCompress, 100);
				mh1 = new Sketch::MinHash(kmerSize, curSketchSize);
			}
			else
				mh1 = new Sketch::MinHash(kmerSize, sketchSize);
		}
		else if(sketchFunc == "KSSD"){
			kssd = new Sketch::KSSD(kssdPara);
		}
		else if(sketchFunc == "WMH"){
			wmh1 = new Sketch::WMinHash(parameter);
		}
		else if(sketchFunc == "HLL"){
			hll = new Sketch::HyperLogLog(HLL_SKETCH_BIT);
		}
		else if(sketchFunc == "OMH"){
			omh = new Sketch::OrderMinHash();
		}
		else{
			fprintf(stderr, "Invaild sketch function: %s\n", sketchFunc.c_str());
			exit(1);
			//return false;
		}


		Vec_SeqInfo curFileSeqs;
		
		while(1){
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;
			string name("noName");
			string comment("noName");
			if(ks1->name.s != NULL)
				name = ks1->name.s;
			if(ks1->comment.s != NULL)
				comment = ks1->comment.s;
			SequenceInfo tmpSeq{name, comment, 0, length};
			
			if(sketchFunc == "MinHash"){
				mh1->update(ks1->seq.s);
			}
			else if(sketchFunc == "KSSD"){
				kssd->update(ks1->seq.s);
			}
			else if(sketchFunc == "WMH"){
				wmh1->update(ks1->seq.s);
			}
			else if(sketchFunc == "HLL"){
				hll->update(ks1->seq.s);
			}
			else if(sketchFunc == "OMH"){
				omh->buildSketch(ks1->seq.s);
			}

			//only save the info of the first sequence for reducing the memory footprint
			//and the overhead of sorting the sketches vector array.
			if(curFileSeqs.size() == 0)
				curFileSeqs.push_back(tmpSeq);
		}//end while, end sketch current file.
		if(sketchFunc == "WMH"){
			wmh1->computeHistoSketch();
		}

		#pragma omp critical
		{
			SketchInfo tmpSketchInfo;
			if(sketchFunc == "MinHash"){
				if(isContainment)
					tmpSketchInfo.isContainment = true;
				tmpSketchInfo.minHash = mh1;
			}
			else if(sketchFunc == "KSSD"){
				tmpSketchInfo.KSSD = kssd;
			}
			else if(sketchFunc == "WMH"){
				tmpSketchInfo.WMinHash = wmh1;
			}
			else if(sketchFunc == "HLL"){
				tmpSketchInfo.HLL = hll;
			}
			else if(sketchFunc == "OMH"){
				tmpSketchInfo.OMH = omh;
			}

			tmpSketchInfo.id = i;
			tmpSketchInfo.fileName = fileList[i];
			tmpSketchInfo.totalSeqLength = totalLength;
			tmpSketchInfo.fileSeqs = curFileSeqs;
			if(totalLength >= minLen)//filter the poor quality genome assemblies whose length less than minLen(fastANI paper)
				sketches.push_back(tmpSketchInfo);
			if(i % 10000 == 0)	cerr << "---finished sketching: " << i << " genomes" << endl;
		}

		gzclose(fp1);
		kseq_destroy(ks1);
	}//end for

	//sort(sketches.begin(), sketches.end(), cmpSketch);
	sort(sketches.begin(), sketches.end(), cmpGenomeSize);

	return true;
}

bool sketchFileWithKssd(const vector<string> fileList, const uint64_t minLen, int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads){
	//fprintf(stderr, "-----input fileList, sketch by file\n");
	//ifstream ifs(inputFile);
	//if(!ifs){
	//	fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
	//	return false;
	//}
	//vector<string> fileList;
	//string fileName;
	//while(getline(ifs, fileName)){
	//	fileList.push_back(fileName);
	//}

	static const int BaseMap[128] = 
	{
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
	-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	// generate the shuffle file.
	int half_k = (kmerSize + 1) / 2;
	kmerSize = half_k * 2;
	bool use64 = half_k - drlevel > 8 ? true : false;
	int half_subk = 6 - drlevel >= 2 ? 6 : drlevel + 2;
	int dim_size = 1 << 4 * half_subk;
	int dim_start = 0;
	int dim_end = 1 << 4 * (half_subk - drlevel);
	int* shuffled_dim = generate_shuffle_dim(half_subk);
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	info.id = (half_k << 8) + (half_subk << 4) + drlevel;
	info.genomeNumber = fileList.size();

	int comp_bittl = 64 - 4 * half_k;
	int half_outctx_len = half_k - half_subk;
	int rev_add_move = 4 * half_k - 2;
	//cout << "the comp_bittl is: " << comp_bittl << endl;

	uint64_t tupmask = _64MASK >> comp_bittl;
	uint64_t domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
	uint64_t undomask = (tupmask ^ domask) & tupmask;
	uint64_t undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
	uint64_t undomask0 = undomask ^ undomask1;

	robin_hood::unordered_map<uint32_t, int> shuffled_map;
	for(int t = 0; t < dim_size; t++){
		if(shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start){
			shuffled_map.insert({t, shuffled_dim[t]});
		}
	}

	// parse the sequences 
	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for (int i = 0; i < fileList.size(); i++){
		//cerr << "start the file: " << fileList[i] << endl;
		gzFile fp1;
		kseq_t* ks1;
		fp1 = gzopen(fileList[i].c_str(), "r");
		if(fp1 == NULL){
			fprintf(stderr, "cannot open the genome file: %s\n", fileList[i].c_str());
			exit(1);
			//return false;
		}
		ks1 = kseq_init(fp1);
		uint64_t totalLength = 0;

		Vec_SeqInfo curFileSeqs;
		unordered_set<uint32_t> hashValueSet;
		unordered_set<uint64_t> hashValueSet64;
		
		while(1){
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;
			string name("noName");
			string comment("noName");
			if(ks1->name.s != NULL)
				name = ks1->name.s;
			if(ks1->comment.s != NULL)
				comment = ks1->comment.s;
			SequenceInfo tmpSeq{name, comment, 0, length};
			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int base = 1;
			//cerr << "the length is: " << length << endl;
			//slide window to generate k-mers and get hashes
			for(int i = 0; i < length; i++)
			{
				//char ch = sequence[i];
				char ch = ks1->seq.s[i];
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;

				}
				if(base > kmerSize)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);

					//pfilter = shuffled_dim[dim_id];
					//if((pfilter >= dim_end) || (pfilter < dim_start)){
					//	continue;
					//}
					
					if(shuffled_map.count(dim_id) == 0){
						continue;
					}
					pfilter = shuffled_map[dim_id];

					pfilter -= dim_start;
					//dr_tuple = (((uni_tuple & undomask0) + ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) + pfilter; 
					////only when the dim_end is 4096(the pfilter is 12bit in binary and the lowerst 12bit is all 0 
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmerSize * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 
					
					if(use64)
						hashValueSet64.insert(dr_tuple);
					else
						hashValueSet.insert(dr_tuple);
				}//end if, i > kmer_size
			}//end for, of a sequence 
			//exit(0);
			
			//only save the info of the first sequence for reducing the memory footprint
			//and the overhead of sorting the sketches vector array.
			if(curFileSeqs.size() == 0)
				curFileSeqs.push_back(tmpSeq);
		}//end while, end sketch current file.

		vector<uint32_t> hashArr32;
		vector<uint64_t> hashArr64;
		
    if(use64){
			for(auto x : hashValueSet64){
				hashArr64.push_back(x);
			}
      std::sort(hashArr64.begin(), hashArr64.end());
		}
		else{
			for(auto x : hashValueSet){
				hashArr32.push_back(x);
			}
      std::sort(hashArr32.begin(), hashArr32.end());
		}

		#pragma omp critical
		{
			KssdSketchInfo tmpKssdSketchInfo;

			tmpKssdSketchInfo.id = i;
			tmpKssdSketchInfo.fileName = fileList[i];
			tmpKssdSketchInfo.totalSeqLength = totalLength;
			tmpKssdSketchInfo.fileSeqs = curFileSeqs;
			tmpKssdSketchInfo.use64 = use64;
			tmpKssdSketchInfo.hash32_arr = hashArr32;
			tmpKssdSketchInfo.hash64_arr = hashArr64;
			tmpKssdSketchInfo.sketchsize = hashArr32.size();
      if(totalLength >= minLen)//filter the poor quality genome assemblies whose length less than minLen(fastANI paper)
				sketches.push_back(tmpKssdSketchInfo);
			if(i % 10000 == 0)	cerr << "---finished sketching: " << i << " genomes" << endl;
		}

		gzclose(fp1);
		kseq_destroy(ks1);
	}//end for
	//sort(sketches.begin(), sketches.end(), KssdcmpSketchSize);
	return true;
}

void transSketches(const vector<KssdSketchInfo>& sketches, const KssdParameters& info, const string folder_path, int numThreads){
	//cerr << "the folder path in transSKetch is: " << folder_path << endl;
	double t0 = get_sec();
	int half_k = info.half_k;
	int drlevel = info.drlevel;
	bool use64 = half_k - drlevel > 8 ? true : false;
	if(use64)
		cerr << "transSketches: use64" << endl;
	else
		cerr << "transSketches: not use64" << endl;

	if(use64){
		double t0 = get_sec();
		robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
		for(size_t i = 0; i < sketches.size(); i++){
			//#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
			for(size_t j = 0; j < sketches[i].hash64_arr.size(); j++){
				uint64_t cur_hash = sketches[i].hash64_arr[j];
				hash_map_arr.insert({cur_hash, vector<uint32_t>()});
				hash_map_arr[cur_hash].push_back(i);
			}
		}
		double t1 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of generate the bloom dictionary and hash_map_arr is: " << t1 - t0 << endl;
		#endif
		size_t hash_number = hash_map_arr.size();
		cerr << "the hash_number is: " << hash_number << endl;
		size_t total_size = 0;
		uint64_t* hash_arr = (uint64_t*)malloc(hash_number * sizeof(uint64_t));
		uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
		FILE* fp_dict = fopen((cur_dict_file).c_str(), "w+");
		if(!fp_dict){
			cerr << "ERROR: transSketches, cannot open dictFile: " << cur_dict_file << endl;
			exit(1);
		}
		size_t cur_id = 0;
		for(auto x : hash_map_arr){
			hash_arr[cur_id] = x.first;
			fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
			hash_size_arr[cur_id] = x.second.size();
			total_size += x.second.size();
			cur_id++;
		}
		cerr << "the total size is: " << total_size << endl;
		fclose(fp_dict);
		double t2 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of writing dictFile is: " << t2 - t1 << endl;
		#endif

		string cur_index_file = folder_path + '/' + "kssd.sketch.index";
		FILE* fp_index = fopen(cur_index_file.c_str(), "w+");
		if(!fp_index){
			cerr << "ERROR: transSketches, cannot open indexFile: " << cur_index_file << endl;
			exit(1);
		}
		fwrite(&hash_number, sizeof(size_t), 1, fp_index);
		fwrite(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		double t3 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of writing indexFile is: " << t3 - t2 << endl;
		#endif
	}
	else{
		size_t hashSize = 1LLU << (4 * (half_k - drlevel));
		vector<vector<uint32_t>> hashMapId;
		for(size_t i = 0; i < hashSize; i++){
			hashMapId.push_back(vector<uint32_t>());
		}
		uint32_t* offsetArr = (uint32_t*)calloc(hashSize, sizeof(uint32_t));

		cerr << "the hashSize is: " << hashSize << endl;
		for(size_t i = 0; i < sketches.size(); i++){
			#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
			for(size_t j = 0; j < sketches[i].hash32_arr.size(); j++){
				uint32_t hash = sketches[i].hash32_arr[j];
				hashMapId[hash].push_back(i);
			}
		}
		double tt0 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of generate the idx by multiple threads are: " << tt0 - t0 << endl;
		#endif

		string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
		FILE * fp0 = fopen(cur_dict_file.c_str(), "w+");
		uint64_t totalIndex = 0;
		for(size_t hash = 0; hash < hashSize; hash++){
			offsetArr[hash] = 0;
			if(hashMapId[hash].size() != 0){
				fwrite(hashMapId[hash].data(), sizeof(uint32_t), hashMapId[hash].size(), fp0);
				totalIndex += hashMapId[hash].size();
				offsetArr[hash] = hashMapId[hash].size();
			}
		}
		fclose(fp0);
		
		double t1 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of merge multiple idx into final hashMap is: " << t1 - tt0 << endl;
		#endif

		string cur_index_file = folder_path + '/' + "kssd.sketch.index";
		FILE * fp1 = fopen(cur_index_file.c_str(), "w+");
		fwrite(&hashSize, sizeof(size_t), 1, fp1);
		fwrite(&totalIndex, sizeof(uint64_t), 1, fp1);
		fwrite(offsetArr, sizeof(uint32_t), hashSize, fp1);
		double t2 = get_sec();
		fclose(fp1);
		#ifdef Timer_inner
		cerr << "the time of write output file is: " << t2 - t1 << endl;
		#endif
	}

	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "the totalIndex is: " << totalIndex << endl;
}





