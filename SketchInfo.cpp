#include "SketchInfo.h"
#include "kseq.h"
#include <fstream>
#include <sstream>
#include <zlib.h>

#ifdef THREADPOOL_MINHASH
#include "ThreadPool.h"
#endif

#ifdef RABBIT_IO
#include "FastxStream.h"
#include "FastxChunk.h"
#include "Formater.h"
#include "DataQueue.h"
#include <thread>
#include <cstdint>
#include <unordered_map>
#endif

KSEQ_INIT(gzFile, gzread);
using namespace std;

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
	sketchInfo.index = input->index;

	output->sketchInfo = sketchInfo;
	//output->index = input->index;
	return output;
}


void useThreadOutput(SketchOutput * output, vector<SketchInfo> &sketches){
	//SketchInfo tmpSimilarityInf
	sketches.push_back(output->sketchInfo);
}

#endif

#ifdef RABBIT_IO
typedef rabbit::core::TDataQueue<rabbit::fa::FastaChunk> FaChunkQueue;
int producer_fasta_task(std::string file, rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq){
	std::cout << "filename" << file << std::endl;
	rabbit::fa::FastaFileReader* faFileReader;
	faFileReader = new rabbit::fa::FastaFileReader(file, *fastaPool, false);
	int n_chunks = 0;
	while(1){
		rabbit::fa::FastaChunk* faChunk = new rabbit::fa::FastaChunk;
		faChunk = faFileReader->readNextChunk();
		if(faChunk == NULL) break;
		n_chunks++;
		dq.Push(n_chunks, faChunk);

	}
	dq.SetCompleted();

	return 0;
}

void consumer_fasta_task(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, string sketchFunc, Sketch::WMHParameters * parameters, vector<SketchInfo> *sketches, vector<SimilarityInfo> *similarityInfos){
	int line_num = 0;
	rabbit::int64 id = 0;

	rabbit::fa::FastaChunk *faChunk;
	while(dq.Pop(id, faChunk)){
		std::vector<Reference> data;
		int ref_num = rabbit::fa::chunkFormat(*faChunk, data);
		for(Reference &r: data){
			SimilarityInfo tmpSimilarityInfo;
			tmpSimilarityInfo.id = r.gid;
			tmpSimilarityInfo.name = r.name;
			tmpSimilarityInfo.comment = r.comment;
			tmpSimilarityInfo.length = r.length;
			similarityInfos->push_back(tmpSimilarityInfo);

			SketchInfo tmpSketchInfo; 
			if(sketchFunc == "MinHash"){
				Sketch::MinHash * mh1 = new Sketch::MinHash(21, 10000);
				mh1->update((char*)r.seq.c_str());
				tmpSketchInfo.minHash = mh1;
			}
			else if(sketchFunc == "WMH"){
				Sketch::WMinHash * wmh = new Sketch::WMinHash(*parameters);
				wmh->update((char*)r.seq.c_str());
				wmh->computeHistoSketch();
				tmpSketchInfo.WMinHash = wmh;
			}
			else if(sketchFunc == "HLL"){
				Sketch::HyperLogLog * hll = new Sketch::HyperLogLog(20);
				hll->update((char*)r.seq.c_str());
				tmpSketchInfo.HLL = hll;
			}
			else if(sketchFunc == "OMH"){
				Sketch::OrderMinHash* omh = new Sketch::OrderMinHash();
				omh->buildSketch((char*)r.seq.c_str());
				tmpSketchInfo.OMH = omh;
			}

			//tmpSketchInfo.minHash = mh1;
			tmpSketchInfo.index = r.gid;
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

void getCWS(double * r, double * c, double * b, int sketchSize, int dimension){
	const int DISTRIBUTION_SEED = 1;
	default_random_engine generator(DISTRIBUTION_SEED);
	gamma_distribution<double> gamma(2.0, 1.0);
	uniform_real_distribution<double> uniform(0.0, 1.0);

	for(int i = 0; i < sketchSize * dimension; ++i){
		r[i] = gamma(generator);
		c[i] = log(gamma(generator));
		b[i] = uniform(generator) * r[i];
	}
}

bool sketchSequences(string inputFile, string sketchFunc, vector<SimilarityInfo>& similarityInfos, vector<SketchInfo>& sketches, int threads){
	gzFile fp1;
	kseq_t * ks1;
	fp1 = gzopen(inputFile.c_str(), "r");
	if(fp1 == NULL){
		fprintf(stderr, "cannot open the genome file, %s\n", inputFile.c_str());
		//printfUsage();
		return false;
	}
	
	ks1 = kseq_init(fp1);
	Sketch::WMHParameters parameters;
	if(sketchFunc == "WMH"){
		//Sketch::WMHParameters parameters;
		parameters.kmerSize = 21;
		parameters.sketchSize = 50;
		parameters.windowSize = 20;
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
		if(length < 21) continue;

		SimilarityInfo tmpSimilarityInfoInfo;
		tmpSimilarityInfoInfo.id = index;
		tmpSimilarityInfoInfo.name = ks1->name.s;
		tmpSimilarityInfoInfo.comment = ks1->comment.s;
		tmpSimilarityInfoInfo.length = length;

		similarityInfos.push_back(tmpSimilarityInfoInfo);
		
		char * seqCopy = (char*) malloc((length+1) * sizeof(char));
		memcpy(seqCopy, ks1->seq.s, length+1);
		SketchInfo sketchInfo;
		if(sketchFunc == "MinHash"){
			Sketch::MinHash * mh1 = new Sketch::MinHash(21, 10000);
			sketchInfo.minHash = mh1;
		}
		else if(sketchFunc == "WMH"){
			Sketch::WMinHash * wmh = new Sketch::WMinHash(parameters);
			sketchInfo.WMinHash = wmh;
		}
		else if(sketchFunc == "HLL"){
			Sketch::HyperLogLog *hll = new Sketch::HyperLogLog(20);
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
	#ifdef RABBIT_IO
	int th = threads - 1;//consumer threads number;
	vector<SketchInfo>  sketchesArr[th];
	vector<SimilarityInfo>  similarityInfosArr[th];
	
	rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1<< 24);
	FaChunkQueue queue1(128, 1);
	//cout << "--------------------" << endl;
	std::thread producer(producer_fasta_task, inputFile, fastaPool, std::ref(queue1));
	std::thread **threadArr = new std::thread* [th];

	for(int t = 0; t < th; t++){
		threadArr[t] = new std::thread(std::bind(consumer_fasta_task, fastaPool, std::ref(queue1), sketchFunc, &parameters, &sketchesArr[t], &similarityInfosArr[t]));
	}
	producer.join();
	for(int t = 0; t < th; t++){
		threadArr[t]->join();
	}

	if(sketchFunc == "MinHash"){
		unordered_map<int, Sketch::MinHash* > sketchMap;
		for(int i = 0; i < th; i++){
			for(int j = 0; j < sketchesArr[i].size(); j++){
				if(!sketchMap.count(sketchesArr[i][j].index)){
					sketchMap.insert({sketchesArr[i][j].index, sketchesArr[i][j].minHash});
				}
				else{
					Sketch::MinHash * tmpMH = sketchMap.at(sketchesArr[i][j].index);
					tmpMH->merge(*sketchesArr[i][j].minHash);
					//sketches.push_back(sketchesArr[i][j]);
				}
			}
		}

		for(auto &x : sketchMap){
			SketchInfo tmpSketch;
			tmpSketch.index = x.first;
			tmpSketch.minHash = x.second;
			sketches.push_back(tmpSketch);
		}

		unordered_map<int, SimilarityInfo> similarityMap;
		for(int i = 0; i < th; i++){
			for(int j = 0; j < similarityInfosArr[i].size(); j++){
				if(!similarityMap.count(similarityInfosArr[i][j].id)){
					similarityMap.insert({similarityInfosArr[i][j].id, similarityInfosArr[i][j]});
				}
				else{
					SimilarityInfo &tmpSimilarityInfo = similarityMap.at(similarityInfosArr[i][j].id);
					tmpSimilarityInfo.name +=similarityInfosArr[i][j].name;
					tmpSimilarityInfo.comment +=similarityInfosArr[i][j].comment;
				}
			}
		}
		for(auto & x : similarityMap){
			similarityInfos.push_back(x.second);
		}
	}
	else{
		cerr << "RabbitIO does not support unmerged MinHash functions" << endl;
		return false;
	}

	cerr << "the size of sketches is: " << sketches.size() << endl;
	cerr << "the size of similarityInfos is: " << similarityInfos.size() << endl;

	
	#else 
	while(1){
		int length = kseq_read(ks1);
		if(length < 0){
			break;
		}
		SimilarityInfo tmpSimilarityInfo;
		tmpSimilarityInfo.id = index;
		tmpSimilarityInfo.name = ks1->name.s;
		tmpSimilarityInfo.comment = ks1->comment.s;
		tmpSimilarityInfo.length = length;

		similarityInfos.push_back(tmpSimilarityInfo);

		SketchInfo tmpSketchInfo;
		if(sketchFunc == "MinHash"){
			Sketch::MinHash *mh1 = new Sketch::MinHash(21, 10000);
			mh1->update(ks1->seq.s);
			tmpSketchInfo.minHash = mh1;
		}
		else if(sketchFunc == "WMH"){
			Sketch::WMinHash *wmh = new Sketch::WMinHash(parameters);
			wmh->update(ks1->seq.s);
			wmh->computeHistoSketch();
			tmpSketchInfo.WMinHash = wmh;
		}
		else if(sketchFunc == "HLL"){
			Sketch::HyperLogLog* hll = new Sketch::HyperLogLog(20);
			hll->update(ks1->seq.s);
			tmpSketchInfo.HLL = hll;

		}
		else if(sketchFunc == "OMH"){
			Sketch::OrderMinHash* omh = new Sketch::OrderMinHash();
			omh->buildSketch(ks1->seq.s);
			tmpSketchInfo.OMH = omh;
		}

		tmpSketchInfo.index = index;
		sketches.push_back(tmpSketchInfo);
		index++;

	}//end while
	#endif
	#endif
	cerr << "the number of sequence is: " << index << endl;
	gzclose(fp1);
	kseq_destroy(ks1);

	return true;
}

bool sketchFiles(string inputFile, string sketchFunc, vector<SimilarityInfo>& similarityInfos, vector<SketchInfo>& sketches, int threads){
	fprintf(stderr, "input fileList, sketch by file\n");
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

	Sketch::WMHParameters parameter;
	if(sketchFunc == "WMH"){
		parameter.kmerSize = 21;
		parameter.sketchSize = 50;
		parameter.windowSize = 20;
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

		Sketch::MinHash * mh1;
		Sketch::WMinHash * wmh1;
		Sketch::HyperLogLog * hll;
		Sketch::OrderMinHash * omh;

		if(sketchFunc == "MinHash"){
			mh1 = new Sketch::MinHash(21, 10000);
		}
		else if(sketchFunc == "WMH"){
			wmh1 = new Sketch::WMinHash(parameter);
		}
		else if(sketchFunc == "HLL"){
			hll = new Sketch::HyperLogLog(10);
		}
		else if(sketchFunc == "OMH"){
			omh = new Sketch::OrderMinHash();
		}
		else{
			fprintf(stderr, "Invaild sketch function: %s\n", sketchFunc.c_str());
			exit(0);
			//return false;
		}


		long long int totalLength = 0;
		string comment("");

		
		while(1){
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;

			comment += '>';
			comment += ks1->name.s;
			comment += ' ';
			comment += ks1->comment.s;
			
			if(sketchFunc == "MinHash"){
				mh1->update(ks1->seq.s);
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

		}
		if(sketchFunc == "WMH"){
			wmh1->computeHistoSketch();
		}
		#pragma omp critical
		{
			SketchInfo tmpSketchInfo;
			if(sketchFunc == "MinHash"){
				tmpSketchInfo.minHash = mh1;
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

			tmpSketchInfo.index = i;
			sketches.push_back(tmpSketchInfo);

			SimilarityInfo tmpSimilarityInfo;
			tmpSimilarityInfo.id = i;
			tmpSimilarityInfo.name = fileList[i];
			tmpSimilarityInfo.comment = comment;
			comment = "";
			tmpSimilarityInfo.length = totalLength;
			similarityInfos.push_back(tmpSimilarityInfo);
			//cerr << "end the file: " << fileList[i] << endl;

		}
		gzclose(fp1);
		kseq_destroy(ks1);
	}//end for

	return true;
}
























