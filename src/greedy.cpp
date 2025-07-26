#ifdef GREEDY_CLUST
#include "greedy.h"
#include <map>
#include <vector>
#include <immintrin.h>
#include <algorithm> // For std::min, std::max, sorting
#include <omp.h>  // For OpenMP parallel processing
#include <iostream>  // For debugging/logging
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


//bool KssdcmpSketchSize(KssdSketchInfo s1, KssdSketchInfo s2){
//	if(s1.sketchsize > s2.sketchsize)	return true;
//	else if(s1.sketchsize == s2.sketchsize)	return s1.id < s2.id;
//	else	return false;
//}
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







//vector<vector<int>> KssdgreedyCluster(vector<KssdSketchInfo>& sketches, int sketch_func_id, double threshold, int threads)
//{
//  int numGenomes = sketches.size();
//  int * clustLabels = new int[numGenomes];
//  memset(clustLabels, 0, numGenomes*sizeof(int));
//  vector<vector<int> > cluster;
//  vector<int> representiveArr;
//  map<int, vector<int> > semiClust;
//  representiveArr.push_back(0);
//  semiClust.insert({0, vector<int>()});
//
//  for(int j = 1; j < numGenomes; j++){
//    map<double, int> distMapCenter;
//    int sizeRef = sketches[j].hash32_arr.size();
//    //std::sort(representiveArr.begin(), representiveArr.end(),KssdcmpSketchSize);
//#pragma omp parallel for num_threads(threads)
//    for(int i = 0; i < representiveArr.size(); i++){
//      int repId = representiveArr[i];
//      double dist;
//      int sizeQry = sketches[repId].hash32_arr.size();
//      if(sizeRef/sizeQry > 1.5 || sizeQry/sizeRef > 1.5)
//        continue;
//
//      dist = distance(sketches[repId].hash32_arr, sketches[j].hash32_arr, 19);
//      if(dist <= threshold){
//        clustLabels[j] = 1;
//#pragma omp critical
//        {
//          distMapCenter.insert({dist, repId});
//        }
//        //break;
//      }
//    }//end for i
//    if(clustLabels[j] == 0){//this genome is a representative genome
//      representiveArr.push_back(j);
//      semiClust.insert({j, vector<int>()});
//    }
//    else{//this genome is a redundant genome, get the nearest representive genome as its center
//      auto it = distMapCenter.begin();
//      int repId = it->second;
//      semiClust[repId].push_back(j);
//    }
//    map<double, int>().swap(distMapCenter);
//    if(j % 10000 == 0) cerr << "---finished cluster: " << j << endl;
//
//  }//end for j
//  //cerr << "the representiveArr size is : " << representiveArr.size() << endl;
//
//  for(auto x : semiClust){
//    int center = x.first;
//    vector<int> redundantArr = x.second;
//    vector<int> curClust;
//    curClust.push_back(center);
//    curClust.insert(curClust.end(), redundantArr.begin(), redundantArr.end());
//    cluster.push_back(curClust);
//  }
//  return cluster;
//}

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
    int * clustLabels = new int[numGenomes];
    memset(clustLabels, 0, numGenomes * sizeof(int));
    std::vector<std::vector<int>> cluster;
    std::vector<int> representiveArr;
    std::map<int, std::vector<int>> semiClust;
   	int radio = calculateMaxSizeRatio(threshold, 19); 
    if (numGenomes == 0) return cluster;
    cerr<< "aaaaaaaaaaaaaaaaaaa"<<radio<<endl;
    representiveArr.push_back(0);
    semiClust.insert({0, std::vector<int>()});
int count_to_print = std::min(100, (int)sketches.size());

for (int i = 0; i < count_to_print; ++i) {
    std::cout << "Sketch[" << i << "].hash32_arr.size() = " 
                  << sketches[i].hash32_arr.size() << std::endl;
                  }


    for(int j = 1; j < numGenomes; j++){
        std::map<double, int> distMapCenter;
        int sizeRef = sketches[j].hash32_arr.size();
        
        std::vector<int> reps_to_blacklist; 

#pragma omp parallel for num_threads(threads)
        for(int i = 0; i < representiveArr.size(); i++){
            int repId = representiveArr[i];
            double dist;
            int sizeQry = sketches[repId].hash32_arr.size();

            if ( std::max(sizeRef, sizeQry) > radio * std::min(sizeRef, sizeQry))
            {
                #pragma omp critical
                {
                    reps_to_blacklist.push_back(repId);
                }
                continue;
            }

            dist = distance(sketches[repId].hash32_arr, sketches[j].hash32_arr, 19);
            if(dist <= threshold){
                clustLabels[j] = 1;
                #pragma omp critical
                {
                    distMapCenter.insert({dist, repId});
                }
            }
        }

        if(clustLabels[j] == 0){
            representiveArr.push_back(j);
            semiClust.insert({j, std::vector<int>()});
        }
        else{
            auto it = distMapCenter.begin();
            int repId = it->second;
            semiClust[repId].push_back(j);
        }
        
        std::map<double, int>().swap(distMapCenter);
        if(j % 10000 == 0) std::cerr << "--- finished cluster: " << j << " | Active reps: " << representiveArr.size() << std::endl;

        if (!reps_to_blacklist.empty()) {
            std::unordered_set<int> blacklist_set(reps_to_blacklist.begin(), reps_to_blacklist.end());

            representiveArr.erase(
                std::remove_if(representiveArr.begin(), representiveArr.end(),
                               [&blacklist_set](int repId) {
                                   return blacklist_set.count(repId);
                               }),
                representiveArr.end()
            );
        }
    }

    delete[] clustLabels;

    for(auto const& [center, redundantArr] : semiClust){
        std::vector<int> curClust;
        curClust.push_back(center);
        curClust.insert(curClust.end(), redundantArr.begin(), redundantArr.end());
        cluster.push_back(curClust);
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
    if(j % 10000 == 0) cerr << "---finished cluster: " << j << " | Active reps: " << representiveArr.size() <<  endl;

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





#endif
