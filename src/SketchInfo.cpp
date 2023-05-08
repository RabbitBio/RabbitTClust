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
bool cmpSeqSize(SketchInfo s1, SketchInfo s2){
	if(s1.seqInfo.length > s2.seqInfo.length)	return true;
	else if(s1.seqInfo.length == s2.seqInfo.length)	return s1.id < s2.id;
	else return false;
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
	cerr << "\t===the number is: " << number << endl;
	cerr << "\t===the badNumber is: " << badNumber << endl;
	cerr << "\t===the totalNumber is: " << totalNumber << endl;
	
	if((double)badNumber / totalNumber >= 0.2)
	{
		fprintf(stderr, "Warning: there are %d poor quality (length < %ld) genome assemblies in the total %d genome assemblied.\n", badNumber, minLen, totalNumber);
	}
	cerr << "\t===the totalSize is: " << totalSize << endl;
	cerr << "\t===the maxSize is: " << maxSize << endl;
	cerr << "\t===the minSize is: " << minSize << endl;
	cerr << "\t===the averageSize is: " << averageSize << endl;
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






