#ifdef GREEDY_CLUST
#include "greedy.h"
#include <map>
#include <unordered_map>  // For unordered_map (O(1) operations)
#include <unordered_set>  // For unordered_set
#include <vector>
#include <immintrin.h>
#include <algorithm> // For std::min, std::max, sorting
#include <limits>    // For numeric_limits
#include <omp.h>  // For OpenMP parallel processing
#include <iostream>  // For debugging/logging
#include <iomanip>   // For std::setprecision
#include <cmath>     // For log, exp
#include <cstring>   // For memset
#include "robin_hood.h"  // For faster hash maps (2-3x speedup over std::unordered_map)
using namespace std;




/* @brief									Generating clusters by greedy incremental algorthm.
 *
 * @param[in] sketches		sketch array including hash values and informations for each genome sketch.
 * @param[in] sketchFunc	sketch Function including MinHash and KSSD
 * @param[in] threshold		distance threshold for cluster, genomes with distance below this threshold are clustered together.
 * @param[in] threads			Thread number for multiThreading
 * @return								cluster result two-dimention array, each array in result is a cluster, and each element in a cluster
 * 												is a genome.
 */


bool aKssdcmpSketchSize(KssdSketchInfo s1, KssdSketchInfo s2){
	if(s1.sketchsize > s2.sketchsize)	return true;
	else if(s1.sketchsize == s2.sketchsize)	return s1.id < s2.id;
	else	return false;
}
uint64_t Kssdu32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3,
		uint64_t *i_a, uint64_t *i_b){
	uint64_t counter=0;
	const uint32_t *end1 = list1+size1, *end2 = list2+size2;
	*i_a = 0;
	*i_b = 0;
	//uint64_t stop = 0;
	// hard to get only the loop instructions, now only a tiny check at the top wrong
	while(list1 != end1 && list2 != end2 ){
		if(*list1 < *list2){
			list1++;
			(*i_a)++;
			size3--;
		}else if(*list1 > *list2){
			list2++; 
			(*i_b)++;
			size3--;
		}else{
			//result[counter++] = *list1;
			counter++;
			list1++; list2++; 
			(*i_a)++;
			(*i_b)++;
			size3--;
		}
		if(size3 == 0) break;
	}
	return counter;
}

uint32_t Kssdu32_intersect_vector_avx512(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t *i_a, uint64_t *i_b){
	//assert(size3 <= size1 + size2);
	uint64_t count=0;
	*i_a = 0;
	*i_b = 0;
	uint64_t st_a = (size1 / 16) * 16;
	uint64_t st_b = (size2 / 16) * 16;
	//	uint64_t stop = (size3 / 16) * 16;

	uint64_t i_a_s, i_b_s;

	if(size3 <= 16){
		count += Kssdu32_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
		return count;
	}

	uint64_t stop = size3 - 16;
	//cout << "stop: " << stop <<  endl;
	__m512i sv0   = _mm512_set_epi32(0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1); //u32_shuffle_vectors[0 ];
	__m512i sv1   = _mm512_set_epi32(1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2); //u32_shuffle_vectors[1 ];
	__m512i sv2   = _mm512_set_epi32(2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3); //u32_shuffle_vectors[2 ];
	__m512i sv3   = _mm512_set_epi32(3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4); //u32_shuffle_vectors[3 ];
	__m512i sv4   = _mm512_set_epi32(4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5); //u32_shuffle_vectors[4 ];
	__m512i sv5   = _mm512_set_epi32(5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6); //u32_shuffle_vectors[5 ];
	__m512i sv6   = _mm512_set_epi32(6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7); //u32_shuffle_vectors[6 ];
	__m512i sv7   = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8); //u32_shuffle_vectors[7 ];
	__m512i sv8   = _mm512_set_epi32(8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9); //u32_shuffle_vectors[8 ];
	__m512i sv9   = _mm512_set_epi32(9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10); //u32_shuffle_vectors[9 ];
	__m512i sv10  = _mm512_set_epi32(10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11); //u32_shuffle_vectors[10];
	__m512i sv11  = _mm512_set_epi32(11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12); //u32_shuffle_vectors[11];
	__m512i sv12  = _mm512_set_epi32(12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13); //u32_shuffle_vectors[12];
	__m512i sv13  = _mm512_set_epi32(13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14); //u32_shuffle_vectors[13];
	__m512i sv14  = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15); //u32_shuffle_vectors[14];
	//__m512i vzero = _mm512_setzero_epi32();
	while(*i_a < st_a && *i_b < st_b){



		uint32_t a_max = list1[*i_a+15];
		uint32_t b_max = list2[*i_b+15];

		__m512i v_a = _mm512_loadu_si512((__m512i*)&(list1[*i_a]));
		__m512i v_b = _mm512_loadu_si512((__m512i*)&(list2[*i_b]));

		*i_a += (a_max <= b_max) * 16;
		*i_b += (a_max >= b_max) * 16;

		__mmask16 cmp0 = _mm512_cmpeq_epu32_mask(v_a, v_b);
		__m512i rot0 = _mm512_permutexvar_epi32(sv0, v_b);
		__mmask16 cmp1 = _mm512_cmpeq_epu32_mask(v_a, rot0);
		__m512i rot1 = _mm512_permutexvar_epi32(sv1, v_b);
		__mmask16 cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot1);
		__m512i rot2 = _mm512_permutexvar_epi32(sv2, v_b);
		__mmask16 cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot2);
		cmp0 = _mm512_kor(_mm512_kor(cmp0, cmp1), _mm512_kor(cmp2, cmp3));

		__m512i rot3 = _mm512_permutexvar_epi32(sv3, v_b);
		__mmask16 cmp4 = _mm512_cmpeq_epu32_mask(v_a, rot3);
		__m512i rot4 = _mm512_permutexvar_epi32(sv4, v_b);
		__mmask16 cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot4);
		__m512i rot5 = _mm512_permutexvar_epi32(sv5, v_b);
		__mmask16 cmp6 = _mm512_cmpeq_epu32_mask(v_a, rot5);
		__m512i rot6 = _mm512_permutexvar_epi32(sv6, v_b);
		__mmask16 cmp7 = _mm512_cmpeq_epu32_mask(v_a, rot6);
		cmp4 = _mm512_kor(_mm512_kor(cmp4, cmp5), _mm512_kor(cmp6, cmp7));

		__m512i rot7 = _mm512_permutexvar_epi32(sv7, v_b);
		cmp1 = _mm512_cmpeq_epu32_mask(v_a, rot7);
		__m512i rot8 = _mm512_permutexvar_epi32(sv8, v_b);
		cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot8);
		__m512i rot9 = _mm512_permutexvar_epi32(sv9, v_b);
		cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot9);
		__m512i rot10 = _mm512_permutexvar_epi32(sv10, v_b);
		cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot10);
		cmp1 = _mm512_kor(_mm512_kor(cmp1, cmp2), _mm512_kor(cmp3, cmp5));

		__m512i rot11 = _mm512_permutexvar_epi32(sv11, v_b);
		cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot11);
		__m512i rot12 = _mm512_permutexvar_epi32(sv12, v_b);
		cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot12);
		__m512i rot13 = _mm512_permutexvar_epi32(sv13, v_b);
		cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot13);
		__m512i rot14 = _mm512_permutexvar_epi32(sv14, v_b);
		cmp6 = _mm512_cmpeq_epu32_mask(v_a, rot14);
		cmp2 = _mm512_kor(_mm512_kor(cmp2, cmp3), _mm512_kor(cmp5, cmp6));


		cmp0 = _mm512_kor(_mm512_kor(cmp0, cmp4), _mm512_kor(cmp1, cmp2));




		//_mm512_mask_compressstoreu_epi64(&result[count], cmp0, v_a);
		//__m512i vres = _mm512_mask_compress_epi64(_mm512_setzero_epi32(), cmp0, v_a);
		//if(cmp0 > 0)
		//	inspect(vres);
		count += _mm_popcnt_u32(cmp0);
		if(*i_a + *i_b - count >= stop){
			count -= _mm_popcnt_u32(cmp0);
			*i_a -= (a_max <= b_max) * 16;
			*i_b -= (a_max >= b_max) * 16;
			break;
		}

	}
	count += Kssdu32_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);

	*i_a += i_a_s;
	*i_b += i_b_s;
	return count;
}




size_t Kssdu32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t* i_a, uint64_t* i_b){
	//assert(size3 <= size1 + size2);
	uint64_t count=0;
	*i_a = 0;
	*i_b = 0;
	uint64_t st_a = (size1 / 8) * 8;
	uint64_t st_b = (size2 / 8) * 8;

	uint64_t i_a_s, i_b_s;

	if(size3 <= 16){
		count += Kssdu32_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
		return count;
	}

	uint64_t stop = size3 - 16;
	while(*i_a < st_a && *i_b < st_b){

		uint32_t a_max = list1[*i_a+7];
		uint32_t b_max = list2[*i_b+7];

		__m256i v_a = _mm256_loadu_si256((__m256i*)&(list1[*i_a]));
		__m256i v_b = _mm256_loadu_si256((__m256i*)&(list2[*i_b]));

		*i_a += (a_max <= b_max) * 8;
		*i_b += (a_max >= b_max) * 8;


		/*constexpr*/ const int32_t cyclic_shift = _MM_SHUFFLE(0,3,2,1); //rotating right
		/*constexpr*/ const int32_t cyclic_shift2= _MM_SHUFFLE(2,1,0,3); //rotating left
		/*constexpr*/ const int32_t cyclic_shift3= _MM_SHUFFLE(1,0,3,2); //between
		__m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a, v_b);
		__m256 rot1 = _mm256_permute_ps((__m256)v_b, cyclic_shift);
		__m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a, (__m256i)rot1);
		__m256 rot2 = _mm256_permute_ps((__m256)v_b, cyclic_shift3);
		__m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a, (__m256i)rot2);
		__m256 rot3 = _mm256_permute_ps((__m256)v_b, cyclic_shift2);
		__m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a, (__m256i)rot3);

		__m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);

		__m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
		__m256 rot5 = _mm256_permute_ps(rot4, cyclic_shift);
		__m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
		__m256 rot6 = _mm256_permute_ps(rot4, cyclic_shift3);
		__m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
		__m256 rot7 = _mm256_permute_ps(rot4, cyclic_shift2);
		__m256i cmp_mask8 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);

		__m256i cmp_mask = _mm256_or_si256(
				_mm256_or_si256(
					_mm256_or_si256(cmp_mask1, cmp_mask2),
					_mm256_or_si256(cmp_mask3, cmp_mask4)
					),
				_mm256_or_si256(
					_mm256_or_si256(cmp_mask5, cmp_mask6),
					_mm256_or_si256(cmp_mask7, cmp_mask8)
					)
				);
		int32_t mask = _mm256_movemask_ps((__m256)cmp_mask);

		count += _mm_popcnt_u32(mask);

		if(*i_a + *i_b - count >= stop){
			//count -= _mm_popcnt_u32(cmp0);
			//*i_a -= (a_max <= b_max) * 16;
			//*i_b -= (a_max >= b_max) * 16;
			break;
		}

	}
	count += Kssdu32_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);

	*i_a += i_a_s;
	*i_b += i_b_s;
	return count;
}




double jaccard(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry)
{
  uint64_t common_elements = 0;

  uint64_t i = 0, j = 0;
  uint32_t sizeRef = hashesRef.size();
  uint32_t sizeQry = hashesQry.size();
#if defined __AVX512F__ && defined __AVX512CD__
  common_elements = Kssdu32_intersect_vector_avx512((uint32_t*)hashesRef.data(), sizeRef, (uint32_t*)hashesQry.data(), sizeQry, sizeRef + sizeQry, &i, &j);

  uint64_t denom = i + j - common_elements;
  if (denom == 0) {
    return 0.0;
  }

  double  jaccard_value = (double)common_elements / denom;
  return jaccard_value;
#elif defined __AVX2__
  common_elements = Kssdu32_intersect_vector_avx2((uint32_t*)hashesRef.data(), sizeRef, (uint32_t*)hashesQry.data(), sizeQry, sizeRef + sizeQry, &i, &j);

  uint64_t denom = i + j - common_elements;
  if (denom == 0) {
    return 0.0;
  }

  double  jaccard_value = (double)common_elements / denom;
  return jaccard_value;
#else
  while (i < sizeRef && j < sizeQry) {
    if (hashesRef[i] < hashesQry[j]) {
      i++;
    } else if (hashesRef[i] > hashesQry[j]) {
      j++;
    } else {
      common_elements++;
      i++;
      j++;
    }
  }
  uint64_t union_elements = static_cast<uint64_t>(sizeRef) + sizeQry - common_elements;

  if (union_elements == 0) {
    return 0.0; 
  }

  return static_cast<double>(common_elements) / union_elements;
#endif
}



double distance(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry, double kmerSize)
{
  double jaccard_ = jaccard(hashesRef, hashesQry);

  if (jaccard_ == 1.0) {
    return 0.0;
  }

  double distance = -log(2 * jaccard_ / (1.0 + jaccard_)) / kmerSize;

  if (distance > 1.0) {
    distance = 1.0;
  }

  return distance;
}



double calculateMaxSizeRatio(double D, int k) {
    if (D < 0) {
        throw std::runtime_error("Mash distance cannot be negative.");
    }
    if (k <= 0) {
        throw std::runtime_error("k-mer size must be positive.");
    }
    
    // Formula derived from inverting Mash and Jaccard calculations:
    // R_max = 2 * e^(D * k) - 1
    return 2.0 * std::exp(D * k) - 1.0;
}


vector<std::vector<int>> KssdgreedyCluster(std::vector<KssdSketchInfo>& sketches, int sketch_func_id, double threshold, int threads)
{
    int numGenomes = sketches.size();
    if (numGenomes == 0) return std::vector<std::vector<int>>();

    // Use modern C++ memory management
    std::vector<int> clustLabels(numGenomes, 0);
    std::vector<std::vector<int>> cluster;
    std::vector<int> representiveArr;
    
    // Use unordered_map for O(1) average insert/lookup instead of O(log n) for map
    std::unordered_map<int, std::vector<int>> semiClust;
    
    double radio = calculateMaxSizeRatio(threshold, 19);
    
    // Optional: Debug output (controlled by environment variable or flag)
    #ifdef DEBUG_SKETCH_DISTRIBUTION
    std::map<int, int> sketchSizeDistribution;
    for (const auto& sketch : sketches) {
        int size = sketch.hash32_arr.size();
        int interval = (size / 1000) * 1000;
        sketchSizeDistribution[interval]++;
    }
    std::cout << "Sketch size distribution (in intervals of 1000):" << std::endl;
    for (const auto& pair : sketchSizeDistribution) {
        std::cout << "[" << pair.first << " - " << (pair.first + 999) << "]: " << pair.second << " sketches" << std::endl;
    }
    #endif

    representiveArr.push_back(0);
    semiClust[0] = std::vector<int>();
    
    // Pre-allocate thread-local storage for best matches to avoid critical section
    std::vector<std::pair<double, int>> thread_best_matches;

    for(int j = 1; j < numGenomes; j++){
        int sizeRef = sketches[j].hash32_arr.size();
        
        // Resize thread_best_matches for current number of threads
        thread_best_matches.assign(threads, {std::numeric_limits<double>::max(), -1});

#pragma omp parallel for num_threads(threads) schedule(dynamic, 64)
        for(int i = 0; i < representiveArr.size(); i++){
            int repId = representiveArr[i];
            int sizeQry = sketches[repId].hash32_arr.size();

            // Bi-directional size ratio filtering (skip if ratio too large in either direction)
            double ratio = (double)sizeQry / sizeRef;
            if (ratio > radio || ratio < 1.0 / radio) {
                continue;
            }

            // Calculate distance
            double dist = distance(sketches[repId].hash32_arr, sketches[j].hash32_arr, 19);
            
            if(dist <= threshold){
                // Use thread-local storage to avoid critical section
                int tid = omp_get_thread_num();
                if(dist < thread_best_matches[tid].first){
                    thread_best_matches[tid] = {dist, repId};
                }
            }
        }

        // Find best match across all threads (serial, but minimal work)
        double best_dist = std::numeric_limits<double>::max();
        int best_rep = -1;
        for(const auto& match : thread_best_matches){
            if(match.second != -1 && match.first < best_dist){
                best_dist = match.first;
                best_rep = match.second;
            }
        }

        if(best_rep != -1){
            // This genome belongs to an existing cluster
            clustLabels[j] = 1;
            semiClust[best_rep].push_back(j);
        }
        else{
            // This genome is a new representative
            representiveArr.push_back(j);
            semiClust[j] = std::vector<int>();
        }
        
        if(j % 10000 == 0) {
            std::cerr << "--- finished cluster: " << j << " | Active reps: " << representiveArr.size() << std::endl;
        }
    }

    // Build final clusters
    cluster.reserve(semiClust.size());
    for(auto const& [center, redundantArr] : semiClust){
        std::vector<int> curClust;
        curClust.reserve(1 + redundantArr.size());
        curClust.push_back(center);
        curClust.insert(curClust.end(), redundantArr.begin(), redundantArr.end());
        cluster.push_back(std::move(curClust));
    }
    
    return cluster;
}





vector<vector<int>> greedyCluster(vector<SketchInfo>& sketches, int sketch_func_id, double threshold, int threads)
{
  int numGenomes = sketches.size();
  int * clustLabels = new int[numGenomes];
  memset(clustLabels, 0, numGenomes*sizeof(int));
  vector<vector<int> > cluster;
  vector<int> representiveArr;
  map<int, vector<int> > semiClust;
  representiveArr.push_back(0);
  semiClust.insert({0, vector<int>()});

  for(int j = 1; j < numGenomes; j++){
    map<double, int> distMapCenter;
#pragma omp parallel for num_threads(threads)
    for(int i = 0; i < representiveArr.size(); i++){
      int repId = representiveArr[i];
      double dist;
      if(sketch_func_id == 0){
        if(sketches[repId].isContainment){
          dist = sketches[repId].minHash->containDistance(sketches[j].minHash);
          //std::cout << "Calculated distance: " << dist << std::endl;
        }
        //dist = 1.0 - sketches[repId].minHash->containJaccard(sketches[j].minHash);
        else
          dist = sketches[repId].minHash->distance(sketches[j].minHash);
      }
      else if(sketch_func_id == 1){
        dist = sketches[repId].KSSD->distance(sketches[j].KSSD);
      }
      else{
        cerr << "can only support MinHash and KSSD with greedy incremental clust" << endl;
        exit(1);
      }
      if(dist <= threshold){
        clustLabels[j] = 1;
#pragma omp critical
        {
          distMapCenter.insert({dist, repId});
        }
        //break;
      }
    }//end for i
    if(clustLabels[j] == 0){//this genome is a representative genome
      representiveArr.push_back(j);
      semiClust.insert({j, vector<int>()});
    }
    else{//this genome is a redundant genome, get the nearest representive genome as its center
      auto it = distMapCenter.begin();
      int repId = it->second;
      semiClust[repId].push_back(j);
    }
    map<double, int>().swap(distMapCenter);
    if(j % 10000 == 0) cerr << "---finished cluster: " << j << endl;

  }//end for j
  //cerr << "the representiveArr size is : " << representiveArr.size() << endl;

  for(auto x : semiClust){
    int center = x.first;
    vector<int> redundantArr = x.second;
    vector<int> curClust;
    curClust.push_back(center);
    curClust.insert(curClust.end(), redundantArr.begin(), redundantArr.end());
    cluster.push_back(curClust);
  }
  return cluster;
}


// ==================== KSSD Greedy Incremental Clustering with Inverted Index ====================

/**
 * Dynamic Inverted Index Manager
 * Maintains hash -> representative sequence IDs mapping
 * Used for fast intersection computation between query and all representatives
 */
class DynamicInvertedIndex {
public:
    // Inverted index: hash -> list of representative IDs containing this hash
    // Use robin_hood for 2-3x faster lookups compared to std::unordered_map
    // Made public for direct access in array-based counting optimization
    robin_hood::unordered_flat_map<uint32_t, std::vector<uint32_t>> index_map;
    
private:
    // Set of representative sequences for fast lookup
    robin_hood::unordered_set<int> representative_set;
    
public:
    DynamicInvertedIndex() {
        // Pre-allocate to avoid rehashing during construction
        index_map.reserve(100000);  // Expected unique hashes
        representative_set.reserve(10000);  // Expected representatives
    }
    
    // Add a new representative to the inverted index
    void add_representative(int rep_id, const std::vector<uint32_t>& hash_array) {
        representative_set.insert(rep_id);
        
        // Insert all hashes of this representative into the index
        for(uint32_t hash : hash_array) {
            index_map[hash].push_back(rep_id);
        }
    }
    
    // Check if a sequence is a representative
    bool is_representative(int id) const {
        return representative_set.count(id) > 0;
    }
    
    /**
     * Compute intersection sizes with all representatives (sparse version)
     * @param hash_array Query sequence's hash array
     * @param intersection_map Output: map[rep_id] = intersection size (only representatives)
     */
    void calculate_intersections_sparse(
        const std::vector<uint32_t>& hash_array,
        std::unordered_map<int, int>& intersection_map) const
    {
        // Clear previous results
        intersection_map.clear();
        
        // Iterate through each hash of current sequence
        for(uint32_t hash : hash_array) {
            auto it = index_map.find(hash);
            if(it != index_map.end()) {
                // Found representatives containing this hash, increment their counts
                // Only representatives are added to map, non-representatives use no memory
                for(uint32_t rep_id : it->second) {
                    intersection_map[rep_id]++;
                }
            }
        }
    }
    
    /**
     * Compute intersection sizes with all representatives (array version for compatibility)
     * @param hash_array Query sequence's hash array
     * @param intersection_counts Output: intersection_counts[rep_id] = intersection size
     */
    void calculate_intersections(
        const std::vector<uint32_t>& hash_array,
        std::vector<int>& intersection_counts) const
    {
        // Reset counters (only reset positions of representatives)
        for(int rep_id : representative_set) {
            intersection_counts[rep_id] = 0;
        }
        
        // Iterate through each hash of current sequence
        for(uint32_t hash : hash_array) {
            auto it = index_map.find(hash);
            if(it != index_map.end()) {
                // Found representatives containing this hash, increment their counts
                // Note: only representative counts will be >0, non-representatives stay 0
                for(uint32_t rep_id : it->second) {
                    intersection_counts[rep_id]++;
                }
            }
        }
    }
    
    size_t num_representatives() const {
        return representative_set.size();
    }
    
    const robin_hood::unordered_set<int>& get_representatives() const {
        return representative_set;
    }
    
    void clear() {
        index_map.clear();
        representative_set.clear();
    }
};

/**
 * Calculate Mash distance
 * Based on Jaccard similarity and k-mer size
 */
inline double calculate_mash_distance_fast(int common, int size0, int size1, int kmer_size) {
    int denom = size0 + size1 - common;
    
    if(size0 == 0 || size1 == 0 || denom == 0) {
        return 1.0;
    }
    
    double jaccard = (double)common / denom;
    
    if(jaccard == 1.0) {
        return 0.0;
    } else if(jaccard == 0.0) {
        return 1.0;
    } else {
        double mashD = (double)-1.0 / kmer_size * log((2 * jaccard) / (1.0 + jaccard));
        return mashD > 1.0 ? 1.0 : mashD;
    }
}

/**
 * @brief KSSD Greedy Incremental Clustering with Inverted Index
 * 
 * Core idea:
 * 1. Dynamically maintain an inverted index containing only representative hashes
 * 2. Outer loop serial processing (ensures correctness, avoids sequence dependency)
 * 3. Use inverted index to compute intersections with all representatives at once
 * 4. Inner loop parallel distance calculation and best match finding
 * 
 * Multi-threading strategy:
 * - Outer loop serial: process each sequence sequentially, ensure deterministic results
 * - Inverted index query serial: but complexity is only O(sketch_size)
 * - Distance calculation parallel: significant speedup when many representatives exist
 * 
 * @param sketches All sequence sketch information
 * @param sketch_func_id Sketch type ID (reserved parameter, currently only KSSD)
 * @param threshold Distance threshold
 * @param threads Number of threads
 * @param kmer_size K-mer size (default 19)
 * @return Clustering results
 */
vector<std::vector<int>> KssdGreedyClusterWithInvertedIndex(
    std::vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size = 19)
{
    int numGenomes = sketches.size();
    if(numGenomes == 0) return std::vector<std::vector<int>>();
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "KSSD Greedy Clustering with Inverted Index" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total genomes: " << numGenomes << std::endl;
    std::cerr << "Threshold: " << threshold << std::endl;
    std::cerr << "Threads: " << threads << std::endl;
    std::cerr << "K-mer size: " << kmer_size << std::endl;
    
    // Initialize data structures
    std::vector<int> clustLabels(numGenomes, 0);
    std::unordered_map<int, std::vector<int>> semiClust;
    std::vector<int> representativeArr;
    
    // Calculate size ratio filter (prefilter impossible matches)
    double size_ratio_filter = calculateMaxSizeRatio(threshold, kmer_size);
    std::cerr << "Size ratio filter: " << size_ratio_filter << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Initialize dynamic inverted index
    DynamicInvertedIndex dynamic_index;
    
    // First sequence as initial representative
    representativeArr.push_back(0);
    semiClust[0] = std::vector<int>();
    dynamic_index.add_representative(0, sketches[0].hash32_arr);
    
    // Thread-local storage to avoid critical sections
    // Store Jaccard instead of distance for best-match selection (avoid log computation)
    struct ThreadBestMatch {
        double jaccard;  // Compare Jaccard directly (monotonic with distance)
        int rep_id;
    };
    std::vector<ThreadBestMatch> thread_best_matches(threads);
    
    // Optimization: Array-based counting instead of unordered_map
    // This eliminates hashing overhead and provides cache-friendly sequential access
    // Memory: ~8-12 bytes per genome (cnt + mark), total ~22-32MB for 2.7M genomes
    std::vector<uint32_t> cnt(numGenomes, 0);      // Intersection count for each rep
    std::vector<uint32_t> mark(numGenomes, 0);     // Mark for current query (epoch)
    std::vector<uint32_t> touched;                 // Representatives with non-zero intersection
    touched.reserve(10000);                        // Pre-allocate for typical case
    uint32_t cur_mark = 0;                         // Current epoch marker
    
    // Candidates will reference touched (no need for separate vector)
    // common will be read from cnt[rep_id]
    
    // Statistics
    uint64_t total_comparisons = 0;
    uint64_t filtered_by_size = 0;
    uint64_t filtered_by_common = 0;  // Filtered by insufficient common hashes
    uint64_t skipped_no_intersection = 0;  // Skipped representatives with no intersection
    
    // Pre-calculate minimum jaccard threshold for common filtering
    // From Mash: D = -1/k * ln(2J/(1+J))
    // Solve for J: 2J/(1+J) = e^(-kD) = x  =>  J = x/(2-x)
    double x = std::exp(-threshold * kmer_size);
    double jaccard_min = x / (2.0 - x);
    
    // ===== Main loop: serial processing of each new sequence =====
    for(int j = 1; j < numGenomes; j++) {
        const std::vector<uint32_t>& query_sketch = sketches[j].hash32_arr;
        int sizeRef = query_sketch.size();
        
        // Step 1: Use inverted index with array-based counting (faster than unordered_map)
        // Time complexity: O(sizeRef), cache-friendly sequential access
        cur_mark++;  // New epoch for this query
        touched.clear();  // Clear touched list from previous query
        
        for(uint32_t hash : query_sketch) {
            auto it = dynamic_index.index_map.find(hash);
            if(it != dynamic_index.index_map.end()) {
                // Found representatives containing this hash
                for(uint32_t rep_id : it->second) {
                    if(mark[rep_id] != cur_mark) {
                        // First time seeing this rep in current query
                        mark[rep_id] = cur_mark;
                        cnt[rep_id] = 1;
                        touched.push_back(rep_id);
                    } else {
                        // Already seen this rep, increment count
                        cnt[rep_id]++;
                    }
                }
            }
        }
        
        // Step 2: Reset thread-local best matches (store Jaccard, not distance)
        for(int t = 0; t < threads; t++) {
            thread_best_matches[t] = {-1.0, -1};  // -1 means no match yet
        }
        
        // Step 3: Parallel best-match finding (compare Jaccard, not distance!)
        // Key insight: Jaccard is monotonic with distance, so we can avoid log computation
        // Larger Jaccard => Smaller distance
        
        // Count skipped representatives
        size_t num_candidates = touched.size();
        size_t num_total_reps = representativeArr.size();
        skipped_no_intersection += (num_total_reps - num_candidates);
        
        // No need to convert to candidates vector - touched is already a vector!
        // Just read common from cnt[rep_id] directly
        
        #pragma omp parallel for num_threads(threads) schedule(dynamic, 64)
        for(size_t i = 0; i < touched.size(); i++) {
            int repId = touched[i];
            int common = cnt[repId];  // Read from array, not from pair
            int sizeQry = sketches[repId].hash32_arr.size();
            
            // Optimization 1: Size filtering (skip if size ratio is too large)
            double ratio = (double)sizeQry / sizeRef;
            if(ratio > size_ratio_filter || ratio < 1.0 / size_ratio_filter) {
                #pragma omp atomic
                filtered_by_size++;
                continue;
            }
            
            // Optimization 2: Common threshold filtering
            // Calculate minimum required common for this pair
            // For standard Mash: jaccard = common / (sizeRef + sizeQry - common)
            // Required: common >= jaccard_min * (sizeRef + sizeQry) / (1 + jaccard_min)
            int common_min = (int)std::ceil(jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min));
            
            if(common < common_min) {
                #pragma omp atomic
                filtered_by_common++;
                continue;
            }
            
            #pragma omp atomic
            total_comparisons++;
            
            // Calculate Jaccard directly (no log/exp needed!)
            // Jaccard = common / (sizeRef + sizeQry - common)
            int denom = sizeRef + sizeQry - common;
            double jaccard = (denom == 0) ? 1.0 : (double)common / denom;
            
            // Compare Jaccard directly (larger Jaccard = smaller distance)
            // Only if Jaccard >= jaccard_min (already filtered by common_min above)
                int tid = omp_get_thread_num();
            if(jaccard > thread_best_matches[tid].jaccard) {
                thread_best_matches[tid].jaccard = jaccard;
                thread_best_matches[tid].rep_id = repId;
            }
        }
        
        // Step 4: Serial merge: find global best match (largest Jaccard)
        double best_jaccard = -1.0;
        int best_rep = -1;
        for(int t = 0; t < threads; t++) {
            if(thread_best_matches[t].rep_id != -1 && 
               thread_best_matches[t].jaccard > best_jaccard) {
                best_jaccard = thread_best_matches[t].jaccard;
                best_rep = thread_best_matches[t].rep_id;
            }
        }
        
        // Step 5: Update clustering
        if(best_rep != -1) {
            // Belongs to existing cluster
            clustLabels[j] = 1;
            semiClust[best_rep].push_back(j);
        } else {
            // Becomes new representative
            representativeArr.push_back(j);
            semiClust[j] = std::vector<int>();
            
            // Critical: add new representative to inverted index
            dynamic_index.add_representative(j, sketches[j].hash32_arr);
        }
        
        // Progress report
        if(j % 5000 == 0 || j == numGenomes - 1) {
            double clustering_rate = 100.0 * (j - representativeArr.size() + 1) / j;
            uint64_t total_filtered = filtered_by_size + filtered_by_common;
            double avg_candidates = (j > 1) ? (double)(total_comparisons + total_filtered) / (j - 1) : 0;
            std::cerr << "Progress: " << j << "/" << numGenomes 
                     << " | Reps: " << representativeArr.size()
                     << " | Clustering: " << std::fixed << std::setprecision(2) << clustering_rate << "%"
                     << " | Comparisons: " << total_comparisons
                     << " | AvgCandidates: " << std::fixed << std::setprecision(0) << avg_candidates
                     << " | SkippedNoInt: " << skipped_no_intersection
                     << " | FilteredByCommon: " << filtered_by_common
                     << std::endl;
        }
    }
    
    // ===== Build final clustering results =====
    std::vector<std::vector<int>> cluster;
    cluster.reserve(semiClust.size());
    
    for(const auto& [center, members] : semiClust) {
        std::vector<int> curClust;
        curClust.reserve(1 + members.size());
        curClust.push_back(center);
        curClust.insert(curClust.end(), members.begin(), members.end());
        cluster.push_back(std::move(curClust));
    }
    
    // ===== Output statistics =====
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "Clustering Completed!" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total clusters: " << cluster.size() << std::endl;
    std::cerr << "Final clustering rate: " 
             << std::fixed << std::setprecision(2)
             << (100.0 * (numGenomes - cluster.size()) / numGenomes) << "%" << std::endl;
    std::cerr << "Total distance comparisons: " << total_comparisons << std::endl;
    std::cerr << "Filtered by size ratio: " << filtered_by_size << std::endl;
    std::cerr << "Filtered by common threshold: " << filtered_by_common << std::endl;
    std::cerr << "Skipped (no intersection): " << skipped_no_intersection << std::endl;
    
    uint64_t total_candidates = skipped_no_intersection + filtered_by_size + filtered_by_common + total_comparisons;
    uint64_t total_skipped = skipped_no_intersection + filtered_by_size + filtered_by_common;
    double skip_rate = total_candidates > 0 ? (100.0 * total_skipped / total_candidates) : 0;
    std::cerr << "Overall skip rate (inverted index + filters): " 
             << std::fixed << std::setprecision(2) << skip_rate << "%" << std::endl;
    
    std::cerr << "Average candidates per query: " 
             << std::fixed << std::setprecision(1)
             << (double)total_candidates / (numGenomes - 1) << std::endl;
    std::cerr << "Average comparisons per query: " 
             << std::fixed << std::setprecision(1)
             << (double)total_comparisons / (numGenomes - 1) << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Cleanup
    dynamic_index.clear();
    
    return cluster;
}


/**
 * @brief MinHash Greedy Incremental Clustering with Inverted Index
 *
 * Core idea:
 * 1. Dynamically maintain an inverted index containing only representative MinHash values
 * 2. Outer loop serial processing (ensures correctness, avoids sequence dependency)
 * 3. Use inverted index to compute hash intersections with all representatives at once
 * 4. Inner loop parallel distance calculation and best match finding
 *
 * Multi-threading strategy:
 * - Outer loop serial: process each sequence sequentially, ensure deterministic results
 * - Inverted index query serial: but complexity is only O(sketch_size)
 * - Distance calculation parallel: significant speedup when many representatives exist
 *
 * @param sketches All sequence sketch information (MinHash)
 * @param sketch_func_id Sketch type ID (should be 0 for MinHash)
 * @param threshold Distance threshold
 * @param threads Number of threads
 * @param kmer_size K-mer size (for size ratio filtering)
 * @return Clustering results
 */
vector<std::vector<int>> MinHashGreedyClusterWithInvertedIndex(
    std::vector<SketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size = 21)
{
    int numGenomes = sketches.size();
    if(numGenomes == 0) return std::vector<std::vector<int>>();

    std::cerr << "\n========================================" << std::endl;
    std::cerr << "MinHash Greedy Clustering with Inverted Index" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total genomes: " << numGenomes << std::endl;
    std::cerr << "Threshold: " << threshold << std::endl;
    std::cerr << "Threads: " << threads << std::endl;
    std::cerr << "K-mer size: " << kmer_size << std::endl;

    // Initialize data structures
    std::vector<int> clustLabels(numGenomes, 0);
    std::unordered_map<int, std::vector<int>> semiClust;
    std::vector<int> representativeArr;

    // Calculate size ratio filter (prefilter impossible matches)
    double size_ratio_filter = calculateMaxSizeRatio(threshold, kmer_size);
    std::cerr << "Size ratio filter: " << size_ratio_filter << std::endl;
    std::cerr << "========================================\n" << std::endl;

    // Dynamic inverted index for MinHash (using uint64_t hashes)
    class MinHashInvertedIndex {
    public:
        // Inverted index: hash -> list of representative IDs containing this hash
        // Use robin_hood for 2-3x faster lookups compared to std::unordered_map
        // Made public for direct access in array-based counting optimization
        robin_hood::unordered_flat_map<uint64_t, std::vector<uint32_t>> index_map;

    private:
        // Set of representative sequences for fast lookup
        robin_hood::unordered_set<int> representative_set;

    public:
        MinHashInvertedIndex() {
            // Pre-allocate to avoid rehashing during construction
            index_map.reserve(100000);  // Expected unique hashes
            representative_set.reserve(10000);  // Expected representatives
        }

        // Add a new representative to the inverted index
        void add_representative(int rep_id, const std::vector<uint64_t>& hash_array) {
            representative_set.insert(rep_id);

            // Insert all hashes of this representative into the index
            for(uint64_t hash : hash_array) {
                index_map[hash].push_back(rep_id);
            }
        }

        // Check if a sequence is a representative
        bool is_representative(int id) const {
            return representative_set.count(id) > 0;
        }

        /**
         * Compute intersection sizes with all representatives (sparse version)
         * @param hash_array Query sequence's hash array
         * @param intersection_map Output: map[rep_id] = intersection size (only representatives)
         */
        void calculate_intersections_sparse(
            const std::vector<uint64_t>& hash_array,
            std::unordered_map<int, int>& intersection_map) const
        {
            // Clear previous results
            intersection_map.clear();

            // Iterate through each hash of current sequence
            for(uint64_t hash : hash_array) {
                auto it = index_map.find(hash);
                if(it != index_map.end()) {
                    // Found representatives containing this hash, increment their counts
                    for(uint32_t rep_id : it->second) {
                        intersection_map[rep_id]++;
                    }
                }
            }
        }

        size_t num_representatives() const {
            return representative_set.size();
        }

        const robin_hood::unordered_set<int>& get_representatives() const {
            return representative_set;
        }

        void clear() {
            index_map.clear();
            representative_set.clear();
        }
    };

    // Initialize dynamic inverted index
    MinHashInvertedIndex dynamic_index;

    // First sequence as initial representative
    representativeArr.push_back(0);
    semiClust[0] = std::vector<int>();
    std::vector<uint64_t> first_hashes = sketches[0].minHash->storeMinHashes();
    dynamic_index.add_representative(0, first_hashes);

    // Thread-local storage to avoid critical sections
    // For standard MinHash, we compare common (or Jaccard) directly instead of distance
    // because Jaccard is monotonic with distance: larger J => smaller dist
    struct ThreadBestMatch {
        int common;      // For standard MinHash: use common directly
        double distance; // For containment: still need distance
        int rep_id;
    };
    std::vector<ThreadBestMatch> thread_best_matches(threads);

    // Pre-compute fixed parameters for standard MinHash (sketch size is constant)
    int fixed_sketch_size = sketches[0].minHash->getSketchSize();
    bool all_fixed_size = true;
    bool all_standard_mode = !sketches[0].isContainment;
    
    // Check if all sketches use standard mode with fixed size
    for(int i = 1; i < std::min(100, numGenomes); i++) { // Sample check
        if(sketches[i].isContainment || sketches[i].minHash->getSketchSize() != fixed_sketch_size) {
            all_fixed_size = false;
            all_standard_mode = false;
            break;
        }
    }
    
    // Pre-compute common_min for fixed sketch size (standard MinHash only)
    int fixed_common_min = 0;
    if(all_fixed_size && all_standard_mode) {
        int actual_kmer_size = sketches[0].minHash->getKmerSize();
        double x = std::exp(-threshold * actual_kmer_size);
        double jaccard_min = x / (2.0 - x);
        // For standard Mash with fixed size: common >= jaccard_min * 2*size / (1 + jaccard_min)
        fixed_common_min = (int)std::ceil(jaccard_min * (2 * fixed_sketch_size) / (1.0 + jaccard_min));
        std::cerr << "Optimization: Fixed sketch size detected (size=" << fixed_sketch_size 
                  << "), pre-computed common_min=" << fixed_common_min << std::endl;
    }
    
    // Optimization: Array-based counting instead of unordered_map
    // This eliminates hashing overhead and provides cache-friendly sequential access
    // Memory: ~8-12 bytes per genome (cnt + mark), total ~22-32MB for 2.7M genomes
    std::vector<uint32_t> cnt(numGenomes, 0);      // Intersection count for each rep
    std::vector<uint32_t> mark(numGenomes, 0);     // Mark for current query (epoch)
    std::vector<uint32_t> touched;                 // Representatives with non-zero intersection
    touched.reserve(10000);                        // Pre-allocate for typical case
    uint32_t cur_mark = 0;                         // Current epoch marker
    
    // Optimization: Reusable query_hashes vector
    std::vector<uint64_t> query_hashes;
    query_hashes.reserve(fixed_sketch_size > 0 ? fixed_sketch_size : 1000);

    // Statistics
    uint64_t total_comparisons = 0;
    uint64_t filtered_by_size = 0;
    uint64_t filtered_by_common = 0;  // Filtered by insufficient common hashes
    uint64_t skipped_no_intersection = 0;  // Skipped representatives with no intersection

    // ===== Main loop: serial processing of each new sequence =====
    for(int j = 1; j < numGenomes; j++) {
        // Reuse query_hashes vector (still need to call storeMinHashes unfortunately)
        query_hashes = sketches[j].minHash->storeMinHashes();
        int sizeRef = query_hashes.size();
        bool isContainment = sketches[j].isContainment;

        // Step 1: Use inverted index with array-based counting (faster than unordered_map)
        // Time complexity: O(sizeRef), cache-friendly sequential access
        cur_mark++;  // New epoch for this query
        touched.clear();  // Clear touched list from previous query
        
        for(uint64_t hash : query_hashes) {
            auto it = dynamic_index.index_map.find(hash);
            if(it != dynamic_index.index_map.end()) {
                // Found representatives containing this hash
                for(uint32_t rep_id : it->second) {
                    if(mark[rep_id] != cur_mark) {
                        // First time seeing this rep in current query
                        mark[rep_id] = cur_mark;
                        cnt[rep_id] = 1;
                        touched.push_back(rep_id);
                    } else {
                        // Already seen this rep, increment count
                        cnt[rep_id]++;
                    }
                }
            }
        }

        // Step 2: Reset thread-local best matches
        for(int t = 0; t < threads; t++) {
            thread_best_matches[t] = {-1, std::numeric_limits<double>::max(), -1};
        }

        // Step 3: Parallel best-match finding (optimized for fixed-size MinHash)
        // Key insight: For standard MinHash with fixed size, Jaccard is monotonic with common
        // So we can compare common directly instead of computing distance!

        // Count skipped representatives
        size_t num_candidates = touched.size();
        size_t num_total_reps = representativeArr.size();
        skipped_no_intersection += (num_total_reps - num_candidates);
        
        // No need to convert to candidates vector - touched is already a vector!
        // Just read common from cnt[rep_id] directly

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 64)
        for(size_t i = 0; i < touched.size(); i++) {
            int repId = touched[i];
            int common = cnt[repId];  // Read from array, not from pair
            int sizeQry = sketches[repId].minHash->getSketchSize();
            bool repIsContainment = sketches[repId].isContainment;

            // Optimization 1: Size filtering (only for containment mode)
            if(repIsContainment || isContainment) {
                double ratio = (double)sizeQry / sizeRef;
                if(ratio > size_ratio_filter || ratio < 1.0 / size_ratio_filter) {
                    #pragma omp atomic
                    filtered_by_size++;
                    continue;
                }
            }

            // Optimization 2: Common threshold filtering
            int common_min;
            if(all_fixed_size && all_standard_mode && !repIsContainment && !isContainment) {
                // Use pre-computed common_min for standard MinHash with fixed size
                common_min = fixed_common_min;
            } else {
                // Dynamic calculation for containment or variable size
                int actual_kmer_size = sketches[repId].minHash->getKmerSize();
                double x = std::exp(-threshold * actual_kmer_size);
                double jaccard_min = x / (2.0 - x);
                
                if(repIsContainment) {
                    common_min = (int)std::ceil(jaccard_min * std::min(sizeRef, sizeQry));
                } else {
                    common_min = (int)std::ceil(jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min));
                }
            }
            
            if(common < common_min) {
                #pragma omp atomic
                filtered_by_common++;
                continue;
            }

            #pragma omp atomic
            total_comparisons++;

            // Best-match selection strategy:
            // For standard MinHash with fixed size: compare common directly (faster!)
            // For containment: still need to compute distance
            int tid = omp_get_thread_num();
            
            if(all_fixed_size && all_standard_mode && !repIsContainment && !isContainment) {
                // Fast path: For fixed-size standard MinHash, larger common => larger Jaccard => smaller distance
                // So just compare common directly, no need to compute distance!
                if(common > thread_best_matches[tid].common) {
                    thread_best_matches[tid].common = common;
                    thread_best_matches[tid].rep_id = repId;
                }
            } else {
                // Slow path: For containment or variable size, need to compute actual distance
                int actual_kmer_size = sketches[repId].minHash->getKmerSize();
                double dist;
                
                if(repIsContainment) {
                    int minSize = std::min(sizeRef, sizeQry);
                    if(minSize == 0) {
                        dist = 1.0;
                    } else {
                        double jaccard = (double)common / minSize;
                        if(jaccard >= 1.0) {
                            dist = 0.0;
                        } else if(jaccard <= 0.0) {
                            dist = 1.0;
                        } else {
                            dist = -log(2.0 * jaccard / (1.0 + jaccard)) / actual_kmer_size;
                            if(dist > 1.0) dist = 1.0;
                        }
                    }
                } else {
                    int denom = sizeRef + sizeQry - common;
                    if(denom == 0) {
                        dist = 0.0;
                    } else {
                        double jaccard = (double)common / denom;
                        if(jaccard >= 1.0) {
                            dist = 0.0;
                        } else if(jaccard <= 0.0) {
                            dist = 1.0;
                        } else {
                            dist = -log(2.0 * jaccard / (1.0 + jaccard)) / actual_kmer_size;
                            if(dist > 1.0) dist = 1.0;
                        }
                    }
                }
                
                if(dist <= threshold && dist < thread_best_matches[tid].distance) {
                    thread_best_matches[tid].distance = dist;
                    thread_best_matches[tid].common = common;
                    thread_best_matches[tid].rep_id = repId;
                }
            }
        }
        
        // Step 4: Serial merge: find global best match
        int best_common = -1;
        double best_dist = std::numeric_limits<double>::max();
        int best_rep = -1;
        
        for(int t = 0; t < threads; t++) {
            if(thread_best_matches[t].rep_id != -1) {
                if(all_fixed_size && all_standard_mode && !isContainment) {
                    // Compare common directly for fixed-size standard MinHash
                    if(thread_best_matches[t].common > best_common) {
                        best_common = thread_best_matches[t].common;
                        best_rep = thread_best_matches[t].rep_id;
                    }
                } else {
                    // Compare distance for containment mode
                    if(thread_best_matches[t].distance < best_dist) {
                best_dist = thread_best_matches[t].distance;
                best_rep = thread_best_matches[t].rep_id;
                    }
                }
            }
        }
        
        // Step 5: Update clustering
        if(best_rep != -1) {
            // Belongs to existing cluster
            clustLabels[j] = 1;
            semiClust[best_rep].push_back(j);
        } else {
            // Becomes new representative
            representativeArr.push_back(j);
            semiClust[j] = std::vector<int>();
            
            // Critical: add new representative to inverted index
            std::vector<uint64_t> new_rep_hashes = sketches[j].minHash->storeMinHashes();
            dynamic_index.add_representative(j, new_rep_hashes);
        }
        
        // Progress report
        if(j % 5000 == 0 || j == numGenomes - 1) {
            double clustering_rate = 100.0 * (j - representativeArr.size() + 1) / j;
            uint64_t total_filtered = filtered_by_size + filtered_by_common;
            double avg_candidates = (j > 1) ? (double)(total_comparisons + total_filtered) / (j - 1) : 0;
            std::cerr << "Progress: " << j << "/" << numGenomes 
                     << " | Reps: " << representativeArr.size()
                     << " | Clustering: " << std::fixed << std::setprecision(2) << clustering_rate << "%"
                     << " | Comparisons: " << total_comparisons
                     << " | AvgCandidates: " << std::fixed << std::setprecision(0) << avg_candidates
                     << " | SkippedNoInt: " << skipped_no_intersection
                     << " | FilteredByCommon: " << filtered_by_common
                     << std::endl;
        }
    }
    
    // ===== Build final clustering results =====
    std::vector<std::vector<int>> cluster;
    cluster.reserve(semiClust.size());
    
    for(const auto& [center, members] : semiClust) {
        std::vector<int> curClust;
        curClust.reserve(1 + members.size());
        curClust.push_back(center);
        curClust.insert(curClust.end(), members.begin(), members.end());
        cluster.push_back(std::move(curClust));
    }
    
    // ===== Output statistics =====
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "Clustering Completed!" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total clusters: " << cluster.size() << std::endl;
    std::cerr << "Final clustering rate: " 
             << std::fixed << std::setprecision(2)
             << (100.0 * (numGenomes - cluster.size()) / numGenomes) << "%" << std::endl;
    std::cerr << "Total distance comparisons: " << total_comparisons << std::endl;
    std::cerr << "Filtered by size ratio: " << filtered_by_size << std::endl;
    std::cerr << "Filtered by common threshold: " << filtered_by_common << std::endl;
    std::cerr << "Skipped (no intersection): " << skipped_no_intersection << std::endl;
    
    uint64_t total_candidates = skipped_no_intersection + filtered_by_size + filtered_by_common + total_comparisons;
    uint64_t total_skipped = skipped_no_intersection + filtered_by_size + filtered_by_common;
    double skip_rate = total_candidates > 0 ? (100.0 * total_skipped / total_candidates) : 0;
    std::cerr << "Overall skip rate (inverted index + filters): "
             << std::fixed << std::setprecision(2) << skip_rate << "%" << std::endl;
    
    std::cerr << "Average candidates per query: " 
             << std::fixed << std::setprecision(1)
             << (double)total_candidates / (numGenomes - 1) << std::endl;
    std::cerr << "Average comparisons per query: "
             << std::fixed << std::setprecision(1)
             << (double)total_comparisons / (numGenomes - 1) << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Cleanup
    dynamic_index.clear();
    
    return cluster;
}

/**
 * @brief Batched version: process multiple sequences at once (experimental)
 * 
 * Note: This version slightly changes algorithm semantics, results may differ from fully serial version
 * Advantage: Sequences within a batch can be processed in parallel, improving parallelism
 * Disadvantage: Sequences within a batch may be suitable as each other's representatives, requires conflict resolution
 * 
 * Use case: When the number of representatives is small initially, better parallel performance
 * 
 * @param batch_size Batch size, recommended 32-128
 */
vector<std::vector<int>> KssdGreedyClusterWithInvertedIndexBatched(
    std::vector<KssdSketchInfo>& sketches,
    int sketch_func_id,
    double threshold,
    int threads,
    int kmer_size = 19,
    int batch_size = 64)
{
    int numGenomes = sketches.size();
    if(numGenomes == 0) return std::vector<std::vector<int>>();
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "BATCHED KSSD Greedy Clustering" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total genomes: " << numGenomes << std::endl;
    std::cerr << "Batch size: " << batch_size << std::endl;
    std::cerr << "Threshold: " << threshold << std::endl;
    std::cerr << "Threads: " << threads << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Initialize
    std::vector<int> clustLabels(numGenomes, 0);
    std::unordered_map<int, std::vector<int>> semiClust;
    DynamicInvertedIndex dynamic_index;
    
    double size_ratio_filter = calculateMaxSizeRatio(threshold, kmer_size);
    
    // First sequence as initial representative
    semiClust[0] = std::vector<int>();
    dynamic_index.add_representative(0, sketches[0].hash32_arr);
    
    // Batch processing result structure
    struct BatchResult {
        int genome_id;
        double best_dist;
        int best_rep;
    };
    
    // Process by batches
    for(int batch_start = 1; batch_start < numGenomes; batch_start += batch_size) {
        int batch_end = std::min(batch_start + batch_size, numGenomes);
        int current_batch_size = batch_end - batch_start;
        
        std::vector<BatchResult> batch_results(current_batch_size);
        
        // ===== Parallel processing within batch =====
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for(int idx = 0; idx < current_batch_size; idx++) {
            int j = batch_start + idx;
            batch_results[idx].genome_id = j;
            batch_results[idx].best_dist = std::numeric_limits<double>::max();
            batch_results[idx].best_rep = -1;
            
            int sizeRef = sketches[j].hash32_arr.size();
            // Optimization: use sparse map, only store representatives, completely skip non-representatives
            std::unordered_map<int, int> intersection_map;
            
            // Compute intersections (note: each thread computes independently, with duplication)
            dynamic_index.calculate_intersections_sparse(sketches[j].hash32_arr, intersection_map);
            
            // Find best representative (directly traverse map, more efficient)
            for(const auto& [repId, common] : intersection_map) {
                int sizeQry = sketches[repId].hash32_arr.size();
                double ratio = (double)sizeQry / sizeRef;
                if(ratio > size_ratio_filter || ratio < 1.0 / size_ratio_filter) {
                    continue;
                }
                
                double dist = calculate_mash_distance_fast(common, sizeRef, sizeQry, kmer_size);
                
                if(dist <= threshold && dist < batch_results[idx].best_dist) {
                    batch_results[idx].best_dist = dist;
                    batch_results[idx].best_rep = repId;
                }
            }
        }
        
        // ===== Serial update of results =====
        // Strategy: sort by distance, larger distance gets priority to become representative
        // This reduces the impact of conflicts within the batch
        std::sort(batch_results.begin(), batch_results.end(),
                 [](const BatchResult& a, const BatchResult& b) {
                     return a.best_dist > b.best_dist;
                 });
        
        for(const auto& result : batch_results) {
            int j = result.genome_id;
            if(result.best_rep != -1) {
                clustLabels[j] = 1;
                semiClust[result.best_rep].push_back(j);
            } else {
                semiClust[j] = std::vector<int>();
                dynamic_index.add_representative(j, sketches[j].hash32_arr);
            }
        }
        
        // Progress report
        if(batch_start % 10000 < batch_size) {
            std::cerr << "Progress: " << batch_end << "/" << numGenomes 
                     << " | Representatives: " << dynamic_index.num_representatives() << std::endl;
        }
    }
    
    // Build final results
    std::vector<std::vector<int>> cluster;
    cluster.reserve(semiClust.size());
    for(const auto& [center, members] : semiClust) {
        std::vector<int> curClust;
        curClust.push_back(center);
        curClust.insert(curClust.end(), members.begin(), members.end());
        cluster.push_back(std::move(curClust));
    }
    
    std::cerr << "\nBatched clustering completed!" << std::endl;
    std::cerr << "Total clusters: " << cluster.size() << std::endl;
    std::cerr << "Clustering rate: " 
             << std::fixed << std::setprecision(2)
             << (100.0 * (numGenomes - cluster.size()) / numGenomes) << "%" << std::endl;
    
    dynamic_index.clear();
    return cluster;
}


#endif

