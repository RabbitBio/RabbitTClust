#include "SketchInfo.h"
#include "kseq.h"
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>

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
#include "phmap.h"  // Google's Swiss Tables for faster hash maps

#include <omp.h>

#define _64MASK 0xffffffffffffffffLLU

KSEQ_INIT(gzFile, gzread);
using namespace std;

string canonicalSketchName(const string& name){
	string key;
	key.reserve(name.size());
	for(unsigned char ch : name){
		if(ch == '-' || ch == '_') continue;
		key.push_back(static_cast<char>(tolower(ch)));
	}
	if(key == "minhash" || key == "mh") return "MinHash";
	if(key == "kssd") return "KSSD";
	if(key == "fastkmv" || key == "kmv") return "FastKMV";
	if(key == "probminhash" || key == "probmh" || key == "pmh") return "ProbMinHash";
	if(key == "bindash") return "BinDash";
	if(key == "hll" || key == "hyperloglog") return "HLL";
	if(key == "setsketch" || key == "ss") return "SetSketch";
	return "";
}

int sketchFunctionId(const string& name){
	const string canonical = canonicalSketchName(name);
	if(canonical == "MinHash") return SKETCH_MINHASH;
	if(canonical == "KSSD") return SKETCH_KSSD;
	if(canonical == "HLL") return SKETCH_HLL;
	if(canonical == "FastKMV") return SKETCH_FASTKMV;
	if(canonical == "ProbMinHash") return SKETCH_PROBMINHASH;
	if(canonical == "BinDash") return SKETCH_BINDASH;
	if(canonical == "SetSketch") return SKETCH_SETSKETCH;
	return -1;
}

bool isGenericRabbitSketch(int sketch_func_id){
	return sketch_func_id == SKETCH_FASTKMV ||
	       sketch_func_id == SKETCH_PROBMINHASH ||
	       sketch_func_id == SKETCH_BINDASH ||
	       sketch_func_id == SKETCH_HLL ||
	       sketch_func_id == SKETCH_SETSKETCH;
}

static double normalizeMashDistance(double distance){
	if(!std::isfinite(distance) || distance >= 1.0) return 1.0;
	if(distance <= 0.0) return 0.0;
	return distance;
}

static double mashDistanceFromJaccard(double jaccard, int kmer_size){
	if(jaccard >= 1.0) return 0.0;
	if(jaccard <= 0.0 || kmer_size <= 0) return 1.0;
	return normalizeMashDistance(
		-log(2.0 * jaccard / (1.0 + jaccard)) /
		static_cast<double>(kmer_size));
}

double sketchDistance(const SketchInfo& lhs, const SketchInfo& rhs, int sketch_func_id){
	switch(sketch_func_id){
		case SKETCH_MINHASH:
			if(!lhs.minHash || !rhs.minHash) return 1.0;
			return normalizeMashDistance(lhs.isContainment
				? lhs.minHash->containDistance(rhs.minHash)
				: lhs.minHash->distance(rhs.minHash));
		case SKETCH_KSSD:
			return lhs.KSSD && rhs.KSSD
				? normalizeMashDistance(lhs.KSSD->distance(rhs.KSSD)) : 1.0;
		case SKETCH_WMH:
			return lhs.WMinHash && rhs.WMinHash
				? normalizeMashDistance(lhs.WMinHash->distance(rhs.WMinHash)) : 1.0;
		case SKETCH_HLL:
			return lhs.HLL && rhs.HLL
				? mashDistanceFromJaccard(lhs.HLL->jaccard_index(*rhs.HLL), lhs.kmerSize)
				: 1.0;
		case SKETCH_OMH:
			return lhs.OMH && rhs.OMH
				? normalizeMashDistance(lhs.OMH->distance(*rhs.OMH)) : 1.0;
		case SKETCH_FASTKMV:
			return lhs.fastKMV && rhs.fastKMV
				? normalizeMashDistance(lhs.fastKMV->distance(*rhs.fastKMV)) : 1.0;
		case SKETCH_PROBMINHASH:
			return lhs.probMinHash && rhs.probMinHash
				? mashDistanceFromJaccard(
					lhs.probMinHash->jaccard_weighted(*rhs.probMinHash), lhs.kmerSize)
				: 1.0;
		case SKETCH_BINDASH:
			return lhs.binDash && rhs.binDash
				? mashDistanceFromJaccard(lhs.binDash->jaccard(*rhs.binDash), lhs.kmerSize)
				: 1.0;
		case SKETCH_SETSKETCH:
			return lhs.setSketch && rhs.setSketch
				? mashDistanceFromJaccard(
					lhs.setSketch->jaccard_index(*rhs.setSketch), lhs.kmerSize)
				: 1.0;
		default:
			return 1.0;
	}
}

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

static int registerBitsForSketchSize(int sketch_size){
	// RabbitSketch's HLL/SetSketch estimators are calibrated for at least
	// 2^10 registers; smaller arrays become unstable on bacterial genomes.
	int bits = 10;
	uint64_t registers = 1ULL << bits;
	while(registers < static_cast<uint64_t>(std::max(sketch_size, 16)) && bits < 20){
		++bits;
		registers <<= 1;
	}
	return bits;
}

static void initializeGenericSketch(
	SketchInfo& info,
	const string& sketch_func,
	int kmer_size,
	int sketch_size,
	bool is_containment,
	int contain_compress,
	uint64_t length_hint)
{
	info.kmerSize = kmer_size;
	info.sketchSize = sketch_size;
	if(sketch_func == "MinHash"){
		const uint64_t desired_size = is_containment
			? std::max<uint64_t>(length_hint / static_cast<uint64_t>(std::max(1, contain_compress)), 100)
			: static_cast<uint64_t>(sketch_size);
		const int actual_size = static_cast<int>(std::min<uint64_t>(desired_size, INT_MAX));
		info.minHash = new Sketch::MinHash(kmer_size, actual_size);
		info.isContainment = is_containment;
	}else if(sketch_func == "FastKMV"){
		info.fastKMV = new Sketch::FastKMV(std::max(sketch_size, 2), kmer_size, 42);
	}else if(sketch_func == "ProbMinHash"){
		info.probMinHash = new Sketch::ProbMinHash4(std::max(sketch_size, 2), kmer_size, 42, 0);
	}else if(sketch_func == "BinDash"){
		const uint32_t groups = static_cast<uint32_t>((std::max(sketch_size, 64) + 63) / 64);
		info.binDash = new Sketch::BinDash(groups, kmer_size, 16, 42);
	}else if(sketch_func == "HLL"){
		info.HLL = new Sketch::HyperLogLog(registerBitsForSketchSize(sketch_size), kmer_size);
	}else if(sketch_func == "SetSketch"){
		info.setSketch = new Sketch::SetSketch(
			registerBitsForSketchSize(sketch_size), 2.0, 5.0, kmer_size, false);
	}
}

static inline uint8_t encodeFrequencyBase(unsigned char ch){
	switch(ch){
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		default: return 0xff;
	}
}

static void accumulateFrequencyKmers(
	const char* seq,
	uint64_t length,
	int kmer_size,
	vector<uint64_t>& sampled,
	uint64_t& total)
{
	static constexpr uint64_t sampling_divisor = 8;
	static constexpr uint64_t sampling_threshold = UINT64_MAX / sampling_divisor;
	static constexpr uint64_t sampling_mix = 0x9E3779B97F4A7C15ULL;
	if(kmer_size < 1 || kmer_size > 32 || length < static_cast<uint64_t>(kmer_size)) return;
	const uint64_t mask = kmer_size == 32 ? UINT64_MAX : ((1ULL << (2 * kmer_size)) - 1);
	const uint64_t high_shift = 2ULL * static_cast<uint64_t>(kmer_size - 1);
	uint64_t forward = 0;
	uint64_t reverse = 0;
	int valid = 0;
	for(uint64_t i = 0; i < length; ++i){
		const uint8_t base = encodeFrequencyBase(static_cast<unsigned char>(seq[i]));
		if(base == 0xff){
			valid = 0;
			forward = reverse = 0;
			continue;
		}
		forward = ((forward << 2) | base) & mask;
		reverse = (reverse >> 2) | (static_cast<uint64_t>(3 ^ base) << high_shift);
		if(++valid >= kmer_size){
			const uint64_t canonical = std::min(forward, reverse);
			++total;
			if(canonical * sampling_mix <= sampling_threshold)
				sampled.push_back(canonical);
		}
	}
}

static void updateGenericSketch(
	SketchInfo& info,
	const string& sketch_func,
	char* seq,
	uint64_t length,
	vector<uint64_t>& frequency_kmers,
	uint64_t& frequency_total)
{
	if(sketch_func == "MinHash") info.minHash->update(seq);
	else if(sketch_func == "FastKMV") info.fastKMV->update(seq, length);
	else if(sketch_func == "ProbMinHash")
		accumulateFrequencyKmers(seq, length, info.kmerSize, frequency_kmers, frequency_total);
	else if(sketch_func == "BinDash") info.binDash->update(seq, length);
	else if(sketch_func == "HLL") info.HLL->update(seq);
	else if(sketch_func == "SetSketch") info.setSketch->update(seq, static_cast<size_t>(length));
}

static void finalizeGenericSketch(
	SketchInfo& info,
	const string& sketch_func,
	vector<uint64_t>& frequency_kmers,
	uint64_t frequency_total)
{
	if(sketch_func == "MinHash"){
		info.minHash->finalize();
	}else if(sketch_func == "FastKMV"){
		info.fastKMV->finalize();
	}else if(sketch_func == "ProbMinHash"){
		sort(frequency_kmers.begin(), frequency_kmers.end());
		if(frequency_total > 0){
			const double inverse_total = 1.0 / static_cast<double>(frequency_total);
			for(size_t i = 0; i < frequency_kmers.size();){
				size_t j = i + 1;
				while(j < frequency_kmers.size() && frequency_kmers[j] == frequency_kmers[i]) ++j;
				info.probMinHash->addHash(
					frequency_kmers[i], static_cast<double>(j - i) * inverse_total);
				i = j;
			}
		}
		info.probMinHash->finalize();
	}else if(sketch_func == "BinDash"){
		info.binDash->finalize();
	}else if(sketch_func == "HLL"){
		// Populate the lazy cardinality cache before pairwise parallel reads.
		(void)info.HLL->cardinality();
	}else if(sketch_func == "SetSketch"){
		// SetSketch's joint MLE consults a lazy cardinality cache. Warm it while
		// this object is still thread-local to avoid concurrent first writes.
		(void)info.setSketch->cardinality();
	}
}

static void releaseSketchObjects(SketchInfo& info){
	delete info.minHash;
	delete info.KSSD;
	delete info.WMinHash;
	delete info.HLL;
	delete info.OMH;
	delete info.fastKMV;
	delete info.probMinHash;
	delete info.binDash;
	delete info.setSketch;
	info.minHash = nullptr;
	info.KSSD = nullptr;
	info.WMinHash = nullptr;
	info.HLL = nullptr;
	info.OMH = nullptr;
	info.fastKMV = nullptr;
	info.probMinHash = nullptr;
	info.binDash = nullptr;
	info.setSketch = nullptr;
}

static void initializeKssdMetadata(
	KssdParameters& info,
	int half_k,
	int half_subk,
	int drlevel,
	size_t genome_number)
{
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	info.id = (half_k << 8) + (half_subk << 4) + drlevel;
	info.genomeNumber = static_cast<int>(genome_number);
}

static void rebuildKssdInvertedIndex(
	vector<KssdSketchInfo>& sketches,
	bool use64,
	KssdInvertedIndex* inverted_index)
{
	if(inverted_index == nullptr) return;
	inverted_index->init(use64);
	for(size_t i = 0; i < sketches.size(); ++i){
		sketches[i].id = static_cast<int>(i);
		if(use64){
			for(uint64_t hash : sketches[i].hash64_arr)
				inverted_index->hash_map_64[hash].push_back(static_cast<uint32_t>(i));
		}else{
			for(uint32_t hash : sketches[i].hash32_arr)
				inverted_index->hash_map_32[hash].push_back(static_cast<uint32_t>(i));
		}
	}
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
	vector<uint64_t> frequency_kmers;
	uint64_t frequency_total = 0;
	if(sketchFunc == "ProbMinHash")
		frequency_kmers.reserve(static_cast<size_t>(sketchInfo.seqInfo.length / 8 + 1));
	if(sketchFunc == "MinHash" || sketchFunc == "FastKMV" ||
	   sketchFunc == "ProbMinHash" || sketchFunc == "BinDash" ||
	   sketchFunc == "HLL" || sketchFunc == "SetSketch"){
		updateGenericSketch(
			sketchInfo, sketchFunc, input->seq,
			static_cast<uint64_t>(sketchInfo.seqInfo.length),
			frequency_kmers, frequency_total);
		finalizeGenericSketch(sketchInfo, sketchFunc, frequency_kmers, frequency_total);
	}
	else if(sketchFunc == "WMH"){
		//wmh->update(ks1->seq.s);
		//wmh->computeHistoSketch();
		sketchInfo.WMinHash->update(input->seq);
		sketchInfo.WMinHash->computeHistoSketch();
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
	if(output == nullptr) return;
	sketches.push_back(std::move(output->sketchInfo));
	delete output;
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

// Retained only for the legacy non-thread-pool build. The current build uses
// sketchBySequence(), which supports every --sketch implementation and the
// current RabbitSketch Kssd API.
#ifndef THREADPOOL_MINHASH
void consumer_fasta_task(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, Sketch::WMHParameters * parameters, vector<SketchInfo> *sketches){
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

void consumer_fasta_task_with_kssd(rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq, int minLen, int kmerSize, int drlevel, const phmap::flat_hash_map<uint32_t, int>* shuffled_map, vector<KssdSketchInfo> *sketches){
	// Thread-Local Memory Pool: Reuse data structures across iterations
	// Theory: Avoid frequent malloc/free, improve cache locality, reduce memory fragmentation
	// Reference: Berger et al. "Hoard: A Scalable Memory Allocator for Multithreaded Applications"
	// Quantifiable: Reduces allocation calls from O(N*sequences) to O(N*threads)
	thread_local vector<uint32_t> reusable_hash32;
	thread_local vector<uint64_t> reusable_hash64;
	thread_local phmap::flat_hash_set<uint32_t> thread_hashSet32;
	thread_local phmap::flat_hash_set<uint64_t> thread_hashSet64;
	
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

			// Reuse thread-local hash sets (clear but keep capacity)
			thread_hashSet32.clear();
			thread_hashSet64.clear();
			
			// Adaptive Pre-allocation: Estimate unique k-mers based on sequence length and dimension reduction
			// Theory: Expected unique k-mers ≈ min(L - k + 1, |Σ|^k) × (dim_end / dim_size) × uniqueness_factor
			// where L is sequence length, k is k-mer size, |Σ| = 4 (DNA alphabet)
			int estimated_kmers = (length >= kmerSize) ? (length - kmerSize + 1) : 0;
			double dim_filter_ratio = (double)(dim_end - dim_start) / dim_size;
			double uniqueness_factor = 0.85;  // Empirical: ~85% k-mers are unique in genomic data
			int estimated_capacity = (int)(estimated_kmers * dim_filter_ratio * uniqueness_factor);
			
			// Reuse existing capacity or expand if needed
			if(use64) {
				if(thread_hashSet64.capacity() < estimated_capacity) {
					thread_hashSet64.reserve(estimated_capacity);
				}
			} else {
				if(thread_hashSet32.capacity() < estimated_capacity) {
					thread_hashSet32.reserve(estimated_capacity);
				}
			}

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
					uint32_t dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);

					// phmap lookup (fast hash map) - single find operation
					auto it = shuffled_map->find(dim_id);
					if(it == shuffled_map->end()) {
						continue;
					}
					pfilter = it->second - dim_start;
					//dr_tuple = (((uni_tuple & undomask0) + ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) + pfilter; 
					////only when the dim_end is 4096(the pfilter is 12bit in binary and the lowerst 12bit is all 0 
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmerSize * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 

					if(use64)
						thread_hashSet64.emplace(dr_tuple);
					else
						thread_hashSet32.emplace(dr_tuple);
				}//end if, i > kmer_size
			}//end for, of a sequence 

			// Memory Pool Reuse: Use thread-local buffers to reduce allocation overhead
			// Strategy: Reuse vector capacity across iterations via swap
			// Benefit: Amortized O(1) allocation after first iteration
			reusable_hash32.clear();
			reusable_hash64.clear();
			
			if(use64){
				// Reuse capacity from previous iteration
				if(reusable_hash64.capacity() < thread_hashSet64.size()) {
					reusable_hash64.reserve(thread_hashSet64.size());
				}
				reusable_hash64.assign(thread_hashSet64.begin(), thread_hashSet64.end());
			}
			else{
				if(reusable_hash32.capacity() < thread_hashSet32.size()) {
					reusable_hash32.reserve(thread_hashSet32.size());
				}
				reusable_hash32.assign(thread_hashSet32.begin(), thread_hashSet32.end());
			}

			tmpKssdSketchInfo.id = r.gid;
			tmpKssdSketchInfo.totalSeqLength = length;
			tmpKssdSketchInfo.use64 = use64;
			// Swap to transfer data while preserving capacity for reuse
			tmpKssdSketchInfo.hash32_arr.swap(reusable_hash32);
			tmpKssdSketchInfo.hash64_arr.swap(reusable_hash64);
			tmpKssdSketchInfo.sketchsize = use64 ? tmpKssdSketchInfo.hash64_arr.size() : tmpKssdSketchInfo.hash32_arr.size();
			sketches->push_back(std::move(tmpKssdSketchInfo));  // Move instead of copy!

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

bool sketchSequencesWithKssd(const string inputFile, const int minLen, const int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index){
	if(drlevel < 0 || drlevel > 3){
		cerr << "ERROR: RabbitSketch KSSD drlevel must be in [0, 3]" << endl;
		return false;
	}
	const int half_k = (kmerSize + 1) / 2;
	const int half_subk = 6;
	const bool use64 = half_k - drlevel > 8;
	auto parameters = std::make_shared<Sketch::kssd_parameter_t>(half_k, half_subk, drlevel);

	gzFile fp = gzopen(inputFile.c_str(), "r");
	if(fp == nullptr){
		cerr << "ERROR: cannot open genome file: " << inputFile << endl;
		return false;
	}
	kseq_t* ks = kseq_init(fp);
	int next_id = 0;
	#pragma omp parallel num_threads(threads)
	{
		while(true){
			int record_id = -1;
			string sequence;
			string name;
			string comment;
			#pragma omp critical(kssd_sequence_reader)
			{
				const int length = kseq_read(ks);
				if(length >= 0){
					record_id = next_id++;
					sequence.assign(ks->seq.s, static_cast<size_t>(length));
					name = ks->name.s ? ks->name.s : "noName";
					comment = ks->comment.s ? ks->comment.s : "noComment";
				}
			}
			if(record_id < 0) break;
			if(sequence.size() < static_cast<size_t>(minLen) ||
			   sequence.size() < static_cast<size_t>(2 * half_k)) continue;

			Sketch::Kssd sketch(parameters);
			sketch.update(sequence.c_str());
			sketch.finalize();
			KssdSketchInfo result;
			result.id = record_id;
			result.totalSeqLength = sequence.size();
			result.seqInfo = {name, comment, 0, static_cast<int>(sequence.size())};
			result.use64 = use64;
			if(use64) result.hash64_arr = sketch.hashList64;
			else result.hash32_arr = sketch.hashList;
			result.sketchsize = use64 ? result.hash64_arr.size() : result.hash32_arr.size();
			#pragma omp critical(kssd_sequence_results)
			sketches.push_back(std::move(result));
		}
	}
	kseq_destroy(ks);
	gzclose(fp);
	sort(sketches.begin(), sketches.end(), [](const KssdSketchInfo& lhs, const KssdSketchInfo& rhs){
		return lhs.id < rhs.id;
	});
	for(size_t i = 0; i < sketches.size(); ++i) sketches[i].id = static_cast<int>(i);
	initializeKssdMetadata(info, half_k, half_subk, drlevel, sketches.size());
	rebuildKssdInvertedIndex(sketches, use64, inverted_index);
	return true;
}


bool sketchSequences(string inputFile, int kmerSize, int sketchSize, int minLen, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads, MinHashInvertedIndex* inverted_index){
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
		if(length < kmerSize || length < minLen) continue;

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
		if(sketchFunc == "MinHash" || sketchFunc == "FastKMV" ||
		   sketchFunc == "ProbMinHash" || sketchFunc == "BinDash" ||
		   sketchFunc == "HLL" || sketchFunc == "SetSketch"){
			initializeGenericSketch(
				sketchInfo, sketchFunc, kmerSize, sketchSize,
				isContainment, containCompress, static_cast<uint64_t>(length));
		}
		else if(sketchFunc == "WMH"){
			Sketch::WMinHash * wmh = new Sketch::WMinHash(parameters);
			sketchInfo.WMinHash = wmh;
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
		threadArr[t] = new std::thread(std::bind(consumer_fasta_task, fastaPool, std::ref(queue1), kmerSize, sketchSize, minLen, sketchFunc, isContainment, containCompress, &parameters, &sketchesArr[t]));
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
		
		// Build inverted index if provided and sketchFunc is MinHash
		if(inverted_index != nullptr && sketchFunc == "MinHash" && tmpSketchInfo.minHash != nullptr) {
			vector<uint64_t> hashes = tmpSketchInfo.minHash->storeMinHashes();
			for(uint64_t hash : hashes) {
				inverted_index->insert_hash_safe(hash, index);
			}
		}
		
		index++;
	}//end while
	#endif
	#endif
	//cerr << "the number of sequence is: " << index << endl;
	gzclose(fp1);
	kseq_destroy(ks1);

	if(inverted_index != nullptr && sketchFunc == "MinHash") {
		// Thread-pool outputs may finish out of order. Build postings against the
		// current vector positions; Greedy remaps them after length ordering,
		// while MST consumes these positions directly. The legacy serial path may
		// already have populated the index, so do not insert duplicates.
		if(inverted_index->hash_map.empty()){
			for(size_t i = 0; i < sketches.size(); i++){
				if(sketches[i].minHash != nullptr) {
					vector<uint64_t> hashes = sketches[i].minHash->storeMinHashes();
					for(uint64_t hash : hashes) {
						inverted_index->hash_map[hash].push_back(i);
					}
				}
			}
		}
	} else {
		sort(sketches.begin(), sketches.end(), cmpSeqSize);
	}

	return true;
}

bool sketchFiles(string inputFile, uint64_t minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo>& sketches, int threads, MinHashInvertedIndex* inverted_index){
	fprintf(stderr, "-----input fileList, sketch by file\n");
	fstream fs(inputFile);
	if(!fs){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
		return false;
	}
	vector<string> fileList;
	string fileName;
	while(getline(fs, fileName)){
		if(!fileName.empty()) fileList.push_back(fileName);
	}

	vector<SketchInfo> built(fileList.size());
	vector<uint8_t> keep(fileList.size(), 0);
	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(size_t i = 0; i < fileList.size(); ++i){
		gzFile fp = gzopen(fileList[i].c_str(), "r");
		if(fp == nullptr){
			#pragma omp critical
			cerr << "ERROR: cannot open genome file: " << fileList[i] << endl;
			continue;
		}

		uint64_t length_hint = 0;
		if(isContainment){
			FILE* raw = fopen(fileList[i].c_str(), "rb");
			if(raw != nullptr){
				if(fileList[i].size() >= 2 && fileList[i].compare(fileList[i].size() - 2, 2, "gz") == 0){
					if(fseek(raw, -4, SEEK_END) == 0){
						uint32_t uncompressed = 0;
						if(fread(&uncompressed, sizeof(uncompressed), 1, raw) == 1)
							length_hint = uncompressed;
					}
				}else if(fseek(raw, 0, SEEK_END) == 0){
					long bytes = ftell(raw);
					if(bytes > 0) length_hint = static_cast<uint64_t>(bytes);
				}
				fclose(raw);
			}
		}

		SketchInfo info;
		info.id = static_cast<int>(i);
		info.fileName = fileList[i];
		initializeGenericSketch(
			info, sketchFunc, kmerSize, sketchSize,
			isContainment, containCompress, length_hint);
		vector<uint64_t> frequency_kmers;
		if(sketchFunc == "ProbMinHash" && length_hint > 0)
			frequency_kmers.reserve(static_cast<size_t>(length_hint / 8 + 1));
		uint64_t frequency_total = 0;

		kseq_t* ks = kseq_init(fp);
		while(true){
			const int length = kseq_read(ks);
			if(length < 0) break;
			info.totalSeqLength += static_cast<uint64_t>(length);
			if(info.fileSeqs.empty()){
				const string name = ks->name.s ? ks->name.s : "noName";
				const string comment = ks->comment.s ? ks->comment.s : "noComment";
				info.fileSeqs.push_back({name, comment, 0, length});
			}
			updateGenericSketch(
				info, sketchFunc, ks->seq.s, static_cast<uint64_t>(length),
				frequency_kmers, frequency_total);
		}
		kseq_destroy(ks);
		gzclose(fp);
		finalizeGenericSketch(info, sketchFunc, frequency_kmers, frequency_total);
		if(info.totalSeqLength >= minLen){
			built[i] = std::move(info);
			keep[i] = 1;
		}else{
			releaseSketchObjects(info);
		}
		if(i % 10000 == 0){
			#pragma omp critical
			cerr << "---finished sketching: " << i << " genomes" << endl;
		}
	}

	for(size_t i = 0; i < built.size(); ++i){
		if(!keep[i]) continue;
		const uint32_t actual_id = static_cast<uint32_t>(sketches.size());
		sketches.push_back(std::move(built[i]));
		if(inverted_index != nullptr && sketches.back().minHash != nullptr){
			for(uint64_t hash : sketches.back().minHash->storeMinHashes())
				inverted_index->hash_map[hash].push_back(actual_id);
		}
	}

	// With an index, keep the original file-list ID until the caller has
	// ordered the sketches and remapped the posting lists. Without an index,
	// retain the existing direct ordering behavior.
	if(inverted_index == nullptr) {
		sort(sketches.begin(), sketches.end(), cmpGenomeSize);
	}

	return true;
}

bool sketchFileWithKssd(const string inputFile, const uint64_t minLen, int kmerSize, const int drlevel, vector<KssdSketchInfo>& sketches, KssdParameters& info, int threads, KssdInvertedIndex* inverted_index){
	if(drlevel < 0 || drlevel > 3){
		cerr << "ERROR: RabbitSketch KSSD drlevel must be in [0, 3]" << endl;
		return false;
	}
	fprintf(stderr, "-----input fileList, sketch by file with RabbitSketch KSSD\n");
	ifstream input(inputFile);
	if(!input){
		cerr << "ERROR: cannot open input file list: " << inputFile << endl;
		return false;
	}
	vector<string> paths;
	string path;
	while(getline(input, path)) if(!path.empty()) paths.push_back(path);

	const int half_k = (kmerSize + 1) / 2;
	const int half_subk = 6;
	const bool use64 = half_k - drlevel > 8;
	auto parameters = std::make_shared<Sketch::kssd_parameter_t>(half_k, half_subk, drlevel);
	vector<KssdSketchInfo> built(paths.size());
	vector<uint8_t> keep(paths.size(), 0);

	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(size_t i = 0; i < paths.size(); ++i){
		gzFile fp = gzopen(paths[i].c_str(), "r");
		if(fp == nullptr){
			#pragma omp critical
			cerr << "ERROR: cannot open genome file: " << paths[i] << endl;
			continue;
		}
		Sketch::Kssd sketch(parameters);
		KssdSketchInfo result;
		result.id = static_cast<int>(i);
		result.fileName = paths[i];
		result.use64 = use64;
		kseq_t* ks = kseq_init(fp);
		while(true){
			const int length = kseq_read(ks);
			if(length < 0) break;
			result.totalSeqLength += static_cast<uint64_t>(length);
			if(result.fileSeqs.empty()){
				const string name = ks->name.s ? ks->name.s : "noName";
				const string comment = ks->comment.s ? ks->comment.s : "noComment";
				result.fileSeqs.push_back({name, comment, 0, length});
			}
			sketch.update(ks->seq.s);
		}
		kseq_destroy(ks);
		gzclose(fp);
		sketch.finalize();
		if(use64) result.hash64_arr = sketch.hashList64;
		else result.hash32_arr = sketch.hashList;
		result.sketchsize = use64 ? result.hash64_arr.size() : result.hash32_arr.size();
		if(result.totalSeqLength >= minLen){
			built[i] = std::move(result);
			keep[i] = 1;
		}
		if(i % 10000 == 0){
			#pragma omp critical
			cerr << "---finished sketching: " << i << " genomes" << endl;
		}
	}

	for(size_t i = 0; i < built.size(); ++i){
		if(!keep[i]) continue;
		built[i].id = static_cast<int>(sketches.size());
		sketches.push_back(std::move(built[i]));
	}
	initializeKssdMetadata(info, half_k, half_subk, drlevel, sketches.size());
	rebuildKssdInvertedIndex(sketches, use64, inverted_index);
	return true;

#if 0
	fprintf(stderr, "-----input fileList, sketch by file\n");
	ifstream ifs(inputFile);
	if(!ifs){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
		return false;
	}
	vector<string> fileList;
	string fileName;
	while(getline(ifs, fileName)){
		fileList.push_back(fileName);
	}

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
	
	// Initialize inverted index if provided
	if(inverted_index != nullptr) {
		size_t hashSize = use64 ? 0 : (1LLU << (4 * (half_k - drlevel)));
		inverted_index->init(use64, hashSize);
	}

	int comp_bittl = 64 - 4 * half_k;
	int half_outctx_len = half_k - half_subk;
	int rev_add_move = 4 * half_k - 2;
	//cout << "the comp_bittl is: " << comp_bittl << endl;

	uint64_t tupmask = _64MASK >> comp_bittl;
	uint64_t domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
	uint64_t undomask = (tupmask ^ domask) & tupmask;
	uint64_t undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
	uint64_t undomask0 = undomask ^ undomask1;

	// Use phmap for fast hash map
	phmap::flat_hash_map<uint32_t, int> shuffled_map;
	shuffled_map.reserve(dim_end - dim_start);  // Pre-allocate
	for(int t = 0; t < dim_size; t++){
		if(shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start){
			shuffled_map.emplace(t, shuffled_dim[t]);
		}
	}

	// Pre-compute constants to avoid repeated calculations in hot loop
	const int kmerSize_x2 = kmerSize * 2;
	const int drlevel_x4 = drlevel * 4;
	const int half_outctx_len_x2 = half_outctx_len * 2;
	const int half_outctx_len_x4 = half_outctx_len * 4;
	const int kmerSize_x2_minus_half_outctx_len_x4 = kmerSize_x2 - half_outctx_len_x4;

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
		phmap::flat_hash_set<uint32_t> hashValueSet;
		phmap::flat_hash_set<uint64_t> hashValueSet64;
		
		// Thread-local reusable vectors to avoid repeated allocation
		vector<uint32_t> hashArr32;
		vector<uint64_t> hashArr64;
		
		// Adaptive Pre-allocation for file-level sketching
		bool first_sequence = true;
		
		while(1){
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;
			
			// Adaptive allocation on first sequence
			if(first_sequence && length > 0) {
				first_sequence = false;
				int estimated_kmers = (length >= kmerSize) ? (length - kmerSize + 1) : 0;
				double dim_filter_ratio = (double)(dim_end - dim_start) / dim_size;
				double uniqueness_factor = 0.85;
				int estimated_capacity = (int)(estimated_kmers * dim_filter_ratio * uniqueness_factor * 1.2);
				
				if(use64) {
					hashValueSet64.reserve(estimated_capacity);
					hashArr64.reserve(estimated_capacity);
				} else {
					hashValueSet.reserve(estimated_capacity);
					hashArr32.reserve(estimated_capacity);
				}
			}
			
			// Optimize string operations: delay creation until needed
			const char* name_ptr = (ks1->name.s != NULL) ? ks1->name.s : "noName";
			const char* comment_ptr = (ks1->comment.s != NULL) ? ks1->comment.s : "noName";
			
			uint64_t tuple = 0LLU, rvs_tuple = 0LLU;
			int base = 1;
			const char* seq_ptr = ks1->seq.s;  // Cache pointer for better locality
			
			//slide window to generate k-mers and get hashes
			// Optimized: reduce branches, improve cache locality
			for(int j = 0; j < length; j++)
			{
				char ch = seq_ptr[j];
				int basenum = BaseMap[(unsigned char)ch];  // Use unsigned to avoid sign extension
				
				// Optimize: combine conditions to reduce branches
				if(__builtin_expect(basenum != -1, 1))  // Branch prediction hint
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^ 3LLU) << rev_add_move); 
					base++;
					
					// Check k-mer completion (moved inside valid base branch)
					if(__builtin_expect(base > kmerSize, 0))  // Unlikely branch
					{
						uint64_t uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
						uint32_t dim_id = (uni_tuple & domask) >> half_outctx_len_x2;

						// Optimized hash map lookup
						auto it = shuffled_map.find(dim_id);
						if(__builtin_expect(it == shuffled_map.end(), 0)) {  // Unlikely to miss
							continue;
						}
						int pfilter = it->second - dim_start;
						uint64_t dr_tuple = (((uni_tuple & undomask0) | 
						                      ((uni_tuple & undomask1) << kmerSize_x2_minus_half_outctx_len_x4)) 
						                      >> drlevel_x4) | (uint64_t)pfilter; 
						
						if(use64)
							hashValueSet64.emplace(dr_tuple);
						else
							hashValueSet.emplace(dr_tuple);
					}
				}
				else{
					base = 1;
					tuple = 0LLU;
					rvs_tuple = 0LLU;
				}
			}
			
			// Only save the info of the first sequence
			if(curFileSeqs.size() == 0) {
				string name(name_ptr);
				string comment(comment_ptr);
				SequenceInfo tmpSeq{name, comment, 0, length};
				curFileSeqs.push_back(tmpSeq);
			}
		}//end while, end sketch current file.

		// Reuse vector capacity: clear and reserve to avoid reallocation
		hashArr32.clear();
		hashArr64.clear();
		
		if(use64){
			if(hashArr64.capacity() < hashValueSet64.size()) {
				hashArr64.reserve(hashValueSet64.size());
			}
			hashArr64.assign(hashValueSet64.begin(), hashValueSet64.end());
			std::sort(hashArr64.begin(), hashArr64.end());
		}
		else{
			if(hashArr32.capacity() < hashValueSet.size()) {
				hashArr32.reserve(hashValueSet.size());
			}
			hashArr32.assign(hashValueSet.begin(), hashValueSet.end());
			std::sort(hashArr32.begin(), hashArr32.end());
		}

		uint32_t actual_id = 0;
		bool sketch_added = false;
		
		#pragma omp critical
		{
			KssdSketchInfo tmpKssdSketchInfo;

			tmpKssdSketchInfo.id = i;
			tmpKssdSketchInfo.fileName = fileList[i];
			tmpKssdSketchInfo.totalSeqLength = totalLength;
			tmpKssdSketchInfo.fileSeqs = curFileSeqs;
			tmpKssdSketchInfo.use64 = use64;
			tmpKssdSketchInfo.hash32_arr = std::move(hashArr32);  // Move instead of copy!
			tmpKssdSketchInfo.hash64_arr = std::move(hashArr64);  // Move instead of copy!
			tmpKssdSketchInfo.sketchsize = use64 ? tmpKssdSketchInfo.hash64_arr.size() : tmpKssdSketchInfo.hash32_arr.size();  // Fix: use correct size
			if(totalLength >= minLen) {  //filter the poor quality genome assemblies whose length less than minLen(fastANI paper)
				actual_id = sketches.size();  // Actual index in sketches vector
				tmpKssdSketchInfo.id = actual_id;
				sketches.push_back(std::move(tmpKssdSketchInfo));  // Move instead of copy!
				sketch_added = true;
				
				// Insert into inverted index while generating sketch (pipeline optimization)
				// Do this inside critical section to ensure thread safety and correct ID
				if(inverted_index != nullptr) {
					if(use64) {
						for(size_t j = 0; j < sketches[actual_id].hash64_arr.size(); j++){
							inverted_index->hash_map_64[sketches[actual_id].hash64_arr[j]].push_back(actual_id);
						}
					} else {
						for(size_t j = 0; j < sketches[actual_id].hash32_arr.size(); j++){
							uint32_t hash = sketches[actual_id].hash32_arr[j];
							inverted_index->hash_map_32[hash].push_back(actual_id);
						}
					}
				}
			}
			if(i % 10000 == 0)	cerr << "---finished sketching: " << i << " genomes" << endl;
		}

		gzclose(fp1);
		kseq_destroy(ks1);
	}//end for
	
	// Note: We don't sort here anymore. The clustering algorithms will sort if needed:
	// - Greedy algorithm sorts by sketch size in KssdGreedyClusterWithInvertedIndex
	// - MST algorithm doesn't require sorting
	// This avoids the expensive index rebuild after sorting.
	// The inverted index is already built correctly during sketch generation.
	
	// Just ensure IDs match positions (they should already be correct from generation)
	if(inverted_index != nullptr) {
		for(size_t i = 0; i < sketches.size(); i++){
			sketches[i].id = i;
		}
	}
	
	return true;
#endif
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
		phmap::flat_hash_map<uint64_t, vector<uint32_t>> hash_map_arr;
		for(size_t i = 0; i < sketches.size(); i++){
			//#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
			for(size_t j = 0; j < sketches[i].hash64_arr.size(); j++){
				uint64_t cur_hash = sketches[i].hash64_arr[j];
				// Optimized: use try_emplace to avoid double lookup
				auto [it, inserted] = hash_map_arr.try_emplace(cur_hash);
				it->second.push_back(i);
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
		// 32-bit: use phmap (sparse) and write same format as transSketchesFromIndex
		phmap::flat_hash_map<uint32_t, vector<uint32_t>> hash_map_arr;
		for(size_t i = 0; i < sketches.size(); i++){
			for(size_t j = 0; j < sketches[i].hash32_arr.size(); j++){
				uint32_t cur_hash = sketches[i].hash32_arr[j];
				hash_map_arr[cur_hash].push_back(i);
			}
		}
		double tt0 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of generate the idx (32-bit phmap) are: " << tt0 - t0 << endl;
		#endif
		size_t hash_number = hash_map_arr.size();
		cerr << "the hash_number (32-bit) is: " << hash_number << endl;
		size_t total_size = 0;
		uint32_t* hash_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
		FILE* fp_dict = fopen(cur_dict_file.c_str(), "w+");
		if(!fp_dict){
			cerr << "ERROR: transSketches, cannot open dictFile: " << cur_dict_file << endl;
			exit(1);
		}
		size_t cur_id = 0;
		for(const auto& x : hash_map_arr){
			hash_arr[cur_id] = x.first;
			fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
			hash_size_arr[cur_id] = x.second.size();
			total_size += x.second.size();
			cur_id++;
		}
		cerr << "the total size (32-bit) is: " << total_size << endl;
		fclose(fp_dict);
		string cur_index_file = folder_path + '/' + "kssd.sketch.index";
		FILE* fp_index = fopen(cur_index_file.c_str(), "w+");
		if(!fp_index){
			cerr << "ERROR: transSketches, cannot open indexFile: " << cur_index_file << endl;
			exit(1);
		}
		fwrite(&hash_number, sizeof(size_t), 1, fp_index);
		fwrite(hash_arr, sizeof(uint32_t), hash_number, fp_index);
		fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		free(hash_arr);
		free(hash_size_arr);
		double t2 = get_sec();
		#ifdef Timer_inner
		cerr << "the time of write output file (32-bit) is: " << t2 - tt0 << endl;
		#endif
	}

	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "the totalIndex is: " << totalIndex << endl;
}

// Write inverted index to files (pipeline version - index already built in memory)
void transSketchesFromIndex(const KssdInvertedIndex& inverted_index, const KssdParameters& info, const string folder_path){
	double t0 = get_sec();
	bool use64 = inverted_index.use64;
	
	if(use64){
		const auto& hash_map_arr = inverted_index.hash_map_64;
		size_t hash_number = hash_map_arr.size();
		cerr << "the hash_number is: " << hash_number << endl;
		
		size_t total_size = 0;
		uint64_t* hash_arr = (uint64_t*)malloc(hash_number * sizeof(uint64_t));
		uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		
		string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
		FILE* fp_dict = fopen((cur_dict_file).c_str(), "w+");
		if(!fp_dict){
			cerr << "ERROR: transSketchesFromIndex, cannot open dictFile: " << cur_dict_file << endl;
			exit(1);
		}
		size_t cur_id = 0;
		for(const auto& x : hash_map_arr){
			hash_arr[cur_id] = x.first;
			fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
			hash_size_arr[cur_id] = x.second.size();
			total_size += x.second.size();
			cur_id++;
		}
		cerr << "the total size is: " << total_size << endl;
		fclose(fp_dict);
		
		string cur_index_file = folder_path + '/' + "kssd.sketch.index";
		FILE* fp_index = fopen(cur_index_file.c_str(), "w+");
		if(!fp_index){
			cerr << "ERROR: transSketchesFromIndex, cannot open indexFile: " << cur_index_file << endl;
			exit(1);
		}
		fwrite(&hash_number, sizeof(size_t), 1, fp_index);
		fwrite(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		free(hash_arr);
		free(hash_size_arr);
	}
	else{
		// 32-bit: write sparse format (like 64-bit) for phmap
		const auto& hash_map_arr = inverted_index.hash_map_32;
		size_t hash_number = hash_map_arr.size();
		cerr << "the hash_number (32-bit) is: " << hash_number << endl;
		
		size_t total_size = 0;
		uint32_t* hash_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		
		string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
		FILE* fp_dict = fopen(cur_dict_file.c_str(), "w+");
		if(!fp_dict){
			cerr << "ERROR: transSketchesFromIndex, cannot open dictFile: " << cur_dict_file << endl;
			exit(1);
		}
		size_t cur_id = 0;
		for(const auto& x : hash_map_arr){
			hash_arr[cur_id] = x.first;
			fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
			hash_size_arr[cur_id] = x.second.size();
			total_size += x.second.size();
			cur_id++;
		}
		cerr << "the total size (32-bit) is: " << total_size << endl;
		fclose(fp_dict);
		
		string cur_index_file = folder_path + '/' + "kssd.sketch.index";
		FILE* fp_index = fopen(cur_index_file.c_str(), "w+");
		if(!fp_index){
			cerr << "ERROR: transSketchesFromIndex, cannot open indexFile: " << cur_index_file << endl;
			exit(1);
		}
		fwrite(&hash_number, sizeof(size_t), 1, fp_index);
		fwrite(hash_arr, sizeof(uint32_t), hash_number, fp_index);
		fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		free(hash_arr);
		free(hash_size_arr);
	}
	
	double t1 = get_sec();
	#ifdef Timer_inner
	cerr << "the time of writing index from pre-built inverted index is: " << t1 - t0 << endl;
	#endif
}

void saveMinHashIndex(const MinHashInvertedIndex& inverted_index, const string folder_path){
	double t0 = get_sec();
	const auto& hash_map = inverted_index.hash_map;
	size_t hash_number = hash_map.size();
	cerr << "-----saving MinHash inverted index: " << hash_number << " unique hashes" << endl;

	size_t total_size = 0;
	uint64_t* hash_arr = (uint64_t*)malloc(hash_number * sizeof(uint64_t));
	uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));

	string cur_dict_file = folder_path + '/' + "minhash.sketch.dict";
	FILE* fp_dict = fopen(cur_dict_file.c_str(), "w+");
	if(!fp_dict){
		cerr << "ERROR: saveMinHashIndex, cannot open dictFile: " << cur_dict_file << endl;
		exit(1);
	}
	size_t cur_id = 0;
	for(const auto& x : hash_map){
		hash_arr[cur_id] = x.first;
		fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
		hash_size_arr[cur_id] = x.second.size();
		total_size += x.second.size();
		cur_id++;
	}
	fclose(fp_dict);

	string cur_index_file = folder_path + '/' + "minhash.sketch.index";
	FILE* fp_index = fopen(cur_index_file.c_str(), "w+");
	if(!fp_index){
		cerr << "ERROR: saveMinHashIndex, cannot open indexFile: " << cur_index_file << endl;
		exit(1);
	}
	fwrite(&hash_number, sizeof(size_t), 1, fp_index);
	fwrite(hash_arr, sizeof(uint64_t), hash_number, fp_index);
	fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
	fclose(fp_index);
	free(hash_arr);
	free(hash_size_arr);

	double t1 = get_sec();
	cerr << "-----saved MinHash inverted index (" << hash_number << " hashes, " << total_size << " entries) in " << t1 - t0 << " seconds" << endl;
}

bool loadMinHashIndex(const string folder_path, MinHashInvertedIndex& inverted_index){
	string index_file = folder_path + '/' + "minhash.sketch.index";
	string dict_file = folder_path + '/' + "minhash.sketch.dict";
	FILE* fp_index = fopen(index_file.c_str(), "rb");
	if(!fp_index) return false;

	size_t hash_number;
	fread(&hash_number, sizeof(size_t), 1, fp_index);
	uint64_t* hash_arr = new uint64_t[hash_number];
	uint32_t* hash_size_arr = new uint32_t[hash_number];
	fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
	fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
	fclose(fp_index);

	FILE* fp_dict = fopen(dict_file.c_str(), "rb");
	if(!fp_dict){
		delete[] hash_arr;
		delete[] hash_size_arr;
		return false;
	}
	uint64_t total_size = 0;
	for(size_t i = 0; i < hash_number; i++) total_size += hash_size_arr[i];
	uint32_t* dict_buf = new uint32_t[total_size];
	fread(dict_buf, sizeof(uint32_t), total_size, fp_dict);
	fclose(fp_dict);

	inverted_index.hash_map.reserve(hash_number);
	size_t off = 0;
	for(size_t i = 0; i < hash_number; i++){
		size_t len = hash_size_arr[i];
		inverted_index.hash_map[hash_arr[i]] = vector<uint32_t>(dict_buf + off, dict_buf + off + len);
		off += len;
	}
	delete[] dict_buf;
	delete[] hash_arr;
	delete[] hash_size_arr;

	cerr << "-----loaded MinHash inverted index: " << hash_number << " unique hashes, " << total_size << " entries" << endl;
	return true;
}
