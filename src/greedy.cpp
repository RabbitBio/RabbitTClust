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
#include <fstream>   // For ifstream, ofstream (binary I/O)
#include <mutex>
#include "phmap.h"  // Google's Swiss Tables (fastest hash map, 3x speedup over std::unordered_map)
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





double jaccard(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry)
{
  uint64_t common_elements = 0;
  uint64_t i = 0, j = 0;
  uint32_t sizeRef = hashesRef.size();
  uint32_t sizeQry = hashesQry.size();
  
  // Use scalar implementation (removed vectorized versions for compatibility)
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
}



double calculate_mash_distance(const std::vector<uint32_t>& hashesRef, const std::vector<uint32_t>& hashesQry, double kmerSize)
{
  double jaccard_ = jaccard(hashesRef, hashesQry);

  if (jaccard_ == 1.0) {
    return 0.0;
  }

  double dist = -log(2 * jaccard_ / (1.0 + jaccard_)) / kmerSize;

  if (dist > 1.0) {
    dist = 1.0;
  }

  return dist;
}

/** KSSD sketch length (hash32_arr or hash64_arr depending on use64). */
static inline size_t kssd_hash_count(const KssdSketchInfo& s) {
  return s.use64 ? s.hash64_arr.size() : s.hash32_arr.size();
}

double jaccard_u64(const std::vector<uint64_t>& hashesRef, const std::vector<uint64_t>& hashesQry)
{
  uint64_t common_elements = 0;
  size_t i = 0, j = 0;
  size_t sizeRef = hashesRef.size();
  size_t sizeQry = hashesQry.size();
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
}

double calculate_mash_distance(const std::vector<uint64_t>& hashesRef, const std::vector<uint64_t>& hashesQry, double kmerSize)
{
  double jaccard_ = jaccard_u64(hashesRef, hashesQry);
  if (jaccard_ == 1.0) {
    return 0.0;
  }
  double dist = -log(2 * jaccard_ / (1.0 + jaccard_)) / kmerSize;
  if (dist > 1.0) {
    dist = 1.0;
  }
  return dist;
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
            double dist = calculate_mash_distance(sketches[repId].hash32_arr, sketches[j].hash32_arr, 19);
            
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
    for(auto const& pair : semiClust){
        int center = pair.first;
        const std::vector<int>& redundantArr = pair.second;
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
    // uint64_t keys unify 32-bit (zero-extended) and 64-bit KSSD sketches
    phmap::flat_hash_map<uint64_t, std::vector<uint32_t>> index_map;
    
private:
    // Set of representative sequences for fast lookup
    phmap::flat_hash_set<int> representative_set;
    // Track minimum query size seen for correct monotonic pruning (handles data fluctuations)
    int min_query_size_seen = INT_MAX;
    
public:
    DynamicInvertedIndex() {
        // Pre-allocate to avoid rehashing during construction
        // For large datasets (e.g., 2.7M genomes), unique hashes can be millions
        // Estimate: ~50K reps × ~1000 hashes/sketch × ~20% uniqueness = ~10M unique hashes
        // Conservative reserve to avoid expensive rehashing
        index_map.reserve(10000000);  // 10M unique hashes (adjust based on your dataset)
        representative_set.reserve(50000);  // Expected representatives (~2% of genomes)
    }
    
    // Add a new representative (32- or 64-bit sketch hashes)
    template<typename T>
    void add_representative(int rep_id, const std::vector<T>& hash_array) {
        representative_set.insert(rep_id);
        for(T hash : hash_array) {
            uint64_t key = static_cast<uint64_t>(hash);
            auto result = index_map.try_emplace(key);
            result.first->second.push_back(static_cast<uint32_t>(rep_id));
        }
    }
    
    // Check if a sequence is a representative
    bool is_representative(int id) const {
        return representative_set.count(id) > 0;
    }
    
    /**
     * Prune representatives whose size makes them theoretically impossible to match future sequences
     * Optimized for decreasing size order: sequences are sorted large->small (but may fluctuate)
     * @param current_size Current query sequence size
     * @param jaccard_min Minimum jaccard threshold
     * @param sketches Reference to sketch array to get sizes
     * @return Number of representatives pruned
     */
    template<typename SketchType>
    size_t prune_too_large_monotonic(int current_size, double jaccard_min, const std::vector<SketchType>& sketches) {
        // Update minimum query size seen (handles fluctuations in sorted data)
        int old_min = min_query_size_seen;
        min_query_size_seen = std::min(min_query_size_seen, current_size);
        
        // For decreasing order: future sequences will be <= min_query_size_seen
        // Max acceptable ratio: sizeQry / sizeRef <= 1 / jaccard_min
        // So: sizeRef <= min_query_size_seen / jaccard_min
        // Add safety margin (1.25x = 0.8 factor)
        int max_acceptable_size = (int)(min_query_size_seen / (jaccard_min * 0.8));
        
        // DEBUG: Scan all representatives to see their size distribution
        int max_rep_size = 0;
        int min_rep_size = INT_MAX;
        int count_above_threshold = 0;
        
        std::vector<int> to_remove;
        
        // Find representatives that are too large (will never match current or future sequences)
        for(int rep_id : representative_set) {
            int rep_size = (int)kssd_hash_count(sketches[rep_id]);
            max_rep_size = std::max(max_rep_size, rep_size);
            min_rep_size = std::min(min_rep_size, rep_size);
            
            if(rep_size > max_acceptable_size) {
                to_remove.push_back(rep_id);
                count_above_threshold++;
            }
        }
        
        // DEBUG: Always print to diagnose (with rep size distribution)
        std::cerr << "  [DEBUG] Rep sizes: min=" << min_rep_size << ", max=" << max_rep_size 
                 << " | Threshold=" << max_acceptable_size 
                 << " | Above threshold: " << count_above_threshold << "/" << representative_set.size() << std::endl;
        
        if(to_remove.empty()) {
            return 0;
        }
        
        // Remove from representative set
        for(int rep_id : to_remove) {
            representative_set.erase(rep_id);
        }
        
        // Remove from inverted index (expensive but necessary)
        for(int rep_id : to_remove) {
            const KssdSketchInfo& sk = sketches[rep_id];
            auto erase_key = [&](uint64_t key) {
                auto it = index_map.find(key);
                if(it != index_map.end()) {
                    auto& posting_list = it->second;
                    uint32_t ru = static_cast<uint32_t>(rep_id);
                    posting_list.erase(
                        std::remove(posting_list.begin(), posting_list.end(), ru),
                        posting_list.end()
                    );
                    if(posting_list.empty()) {
                        index_map.erase(it);
                    }
                }
            };
            if(sk.use64) {
                for(uint64_t h : sk.hash64_arr) erase_key(h);
            } else {
                for(uint32_t h : sk.hash32_arr) erase_key(static_cast<uint64_t>(h));
            }
        }
        
        return to_remove.size();
    }
    
    /**
     * DEPRECATED: Only kept for batched version compatibility
     * Main inverted index functions use array-based counting instead
     * 
     * Compute intersection sizes with all representatives (sparse version)
     * @param hash_array Query sequence's hash array
     * @param intersection_map Output: map[rep_id] = intersection size (only representatives)
     */
    template<typename Vec>
    void calculate_intersections_sparse(
        const Vec& hash_array,
        std::unordered_map<int, int>& intersection_map) const
    {
        intersection_map.clear();
        for(auto hash : hash_array) {
            uint64_t key = static_cast<uint64_t>(hash);
            auto it = index_map.find(key);
            if(it != index_map.end()) {
                for(uint32_t rep_id : it->second) {
                    intersection_map[static_cast<int>(rep_id)]++;
                }
            }
        }
    }
    
    size_t num_representatives() const {
        return representative_set.size();
    }
    
    const phmap::flat_hash_set<int>& get_representatives() const {
        return representative_set;
    }
    
    int get_min_query_size_seen() const {
        return min_query_size_seen;
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
    const bool use64 = sketches[0].use64;
    for(int vi = 1; vi < numGenomes; vi++) {
        if(sketches[vi].use64 != use64) {
            std::cerr << "ERROR: Mixed KSSD hash width (use64) in dataset at genome " << vi << std::endl;
            exit(1);
        }
    }
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "KSSD Greedy Clustering with Inverted Index" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total genomes: " << numGenomes << std::endl;
    std::cerr << "Threshold: " << threshold << std::endl;
    std::cerr << "Threads: " << threads << std::endl;
    std::cerr << "K-mer size: " << kmer_size << std::endl;
    std::cerr << "KSSD hash storage: " << (use64 ? "64-bit" : "32-bit") << std::endl;
    
    // ===== Sort sketches by size (descending) for monotonic pruning =====
    std::cerr << "\nSorting sketches by size (descending order)..." << std::endl;
    std::sort(sketches.begin(), sketches.end(), 
              [](const KssdSketchInfo& a, const KssdSketchInfo& b) {
                  return kssd_hash_count(a) > kssd_hash_count(b);
              });
    std::cerr << "Sorting completed: [" << kssd_hash_count(sketches[0])
             << " ... " << kssd_hash_count(sketches[numGenomes-1]) << "]" << std::endl;
    
    // Initialize data structures
    std::vector<int> clustLabels(numGenomes, 0);
    
    // Optimized cluster storage: use array + vector<vector<int>> instead of unordered_map
    // rep2cid[genome_id] = cluster_id (-1 if not a representative)
    std::vector<int> rep2cid(numGenomes, -1);
    std::vector<int> representativeArr;  // List of representative genome IDs
    std::vector<std::vector<int>> cluster_members;  // cluster_members[cid] = member list (without center)
    
    std::cerr << "========================================\n" << std::endl;
    
    // Initialize dynamic inverted index
    DynamicInvertedIndex dynamic_index;
    
    // First sequence as initial representative
    representativeArr.push_back(0);
    rep2cid[0] = 0;  // First representative gets cluster ID 0
    cluster_members.emplace_back();  // Empty member list for first cluster
    if(use64)
        dynamic_index.add_representative(0, sketches[0].hash64_arr);
    else
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
    uint64_t filtered_by_common = 0;  // Filtered by insufficient common hashes
    uint64_t skipped_no_intersection = 0;  // Skipped representatives with no intersection
    
    // Pre-calculate minimum jaccard threshold for common filtering
    // From Mash: D = -1/k * ln(2J/(1+J))
    // Solve for J: 2J/(1+J) = e^(-kD) = x  =>  J = x/(2-x)
    double x = std::exp(-threshold * kmer_size);
    double jaccard_min = x / (2.0 - x);
    
    // Adaptive pruning interval based on dataset size
    // Small datasets (~200K): prune every 100K (more frequent)
    // Large datasets (~2M+): prune every 1M (less overhead)
    int prune_interval = (numGenomes < 500000) ? 100000 : 1000000;
    int prune_start = prune_interval;  // Start pruning after first interval
    
    // ===== DEBUG: Print sketch sizes to understand data distribution =====
    std::cerr << "\n=== Data Distribution Analysis ===" << std::endl;
    std::cerr << "First sequence sketch size: " << kssd_hash_count(sketches[0]) << std::endl;
    if(numGenomes >= 100) {
        std::cerr << "Sequence 100 sketch size: " << kssd_hash_count(sketches[99]) << std::endl;
    }
    if(numGenomes >= 1000) {
        std::cerr << "Sequence 1000 sketch size: " << kssd_hash_count(sketches[999]) << std::endl;
    }
    if(numGenomes >= 10000) {
        std::cerr << "Sequence 10000 sketch size: " << kssd_hash_count(sketches[9999]) << std::endl;
    }
    if(numGenomes >= 100000) {
        std::cerr << "Sequence 100000 sketch size: " << kssd_hash_count(sketches[99999]) << std::endl;
    }
    if(numGenomes >= 1000000) {
        std::cerr << "Sequence 1000000 sketch size: " << kssd_hash_count(sketches[999999]) << std::endl;
    }
    std::cerr << "Last sequence sketch size: " << kssd_hash_count(sketches[numGenomes-1]) << std::endl;
    std::cerr << "===================================\n" << std::endl;
    
    // ===== Main loop: serial processing of each new sequence =====
    for(int j = 1; j < numGenomes; j++) {
        int sizeRef = (int)kssd_hash_count(sketches[j]);
        
        // Periodic pruning: remove representatives that are too large (exploiting decreasing size order)
        // Since sequences are sorted large->small, representatives that are too large for current
        // sequence will NEVER match any future sequence (which are even smaller)
        // Frequency adapts to dataset size for optimal performance
        if(j >= prune_start && j % prune_interval == 0) {
            int before_prune = dynamic_index.num_representatives();
            size_t pruned = dynamic_index.prune_too_large_monotonic(sizeRef, jaccard_min, sketches);
            
            // Always report pruning (even if pruned == 0, helps debugging)
            int min_seen = dynamic_index.get_min_query_size_seen();
            int max_acceptable = (int)(min_seen / (jaccard_min * 0.8));
            std::cerr << "\n=== Monotonic Pruning at iter " << j << " ===" << std::endl;
            std::cerr << "  Before: " << before_prune << " reps | After: " << dynamic_index.num_representatives() << " reps" << std::endl;
            std::cerr << "  Removed: " << pruned << " too-large reps" << std::endl;
            std::cerr << "  Threshold size: " << max_acceptable << " (min_seen=" << min_seen << ", current=" << sizeRef << ")" << std::endl;
            std::cerr << "  Jaccard_min: " << std::fixed << std::setprecision(6) << jaccard_min << std::endl;
        }
        
        // Step 1: Use inverted index with array-based counting (faster than unordered_map)
        // Time complexity: O(sizeRef), cache-friendly sequential access
        cur_mark++;  // New epoch for this query
        touched.clear();  // Clear touched list from previous query
        
        auto touch_from_posting = [&](uint32_t rep_id) {
            if(mark[rep_id] != cur_mark) {
                mark[rep_id] = cur_mark;
                cnt[rep_id] = 1;
                touched.push_back(static_cast<int>(rep_id));
            } else {
                cnt[rep_id]++;
            }
        };
        if(use64) {
            for(uint64_t hash : sketches[j].hash64_arr) {
                auto it = dynamic_index.index_map.find(hash);
                if(it != dynamic_index.index_map.end()) {
                    for(uint32_t rep_id : it->second) touch_from_posting(rep_id);
                }
            }
        } else {
            for(uint32_t h : sketches[j].hash32_arr) {
                auto it = dynamic_index.index_map.find(static_cast<uint64_t>(h));
                if(it != dynamic_index.index_map.end()) {
                    for(uint32_t rep_id : it->second) touch_from_posting(rep_id);
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
        size_t num_total_reps = dynamic_index.num_representatives();
        skipped_no_intersection += (num_total_reps - num_candidates);
        
        // No need to convert to candidates vector - touched is already a vector!
        // Just read common from cnt[rep_id] directly
        
        // Thread-local counters to avoid atomic contention in hot loop
        #pragma omp parallel num_threads(threads)
        {
            // Get thread ID once per parallel region (not per iteration!)
            int tid = omp_get_thread_num();
            
            uint64_t local_comparisons = 0;
            uint64_t local_filtered_common = 0;
            
            // Thread-local best match
            double local_best_jaccard = -1.0;
            int local_best_rep = -1;
            
            #pragma omp for schedule(dynamic, 64)
            for(size_t i = 0; i < touched.size(); i++) {
                int repId = touched[i];
                int common = cnt[repId];  // Read from array, not from pair
                int sizeQry = (int)kssd_hash_count(sketches[repId]);
                
                // Optimization: Common threshold filtering
                // Calculate minimum required common for this pair
                // For standard Mash: jaccard = common / (sizeRef + sizeQry - common)
                // Required: common >= jaccard_min * (sizeRef + sizeQry) / (1 + jaccard_min)
                int common_min = (int)std::ceil(jaccard_min * (sizeRef + sizeQry) / (1.0 + jaccard_min));
                
                if(common < common_min) {
                    local_filtered_common++;
                    continue;
                }
                
                local_comparisons++;
                
                // Calculate Jaccard directly (no log/exp needed!)
                // Jaccard = common / (sizeRef + sizeQry - common)
                int denom = sizeRef + sizeQry - common;
                double jaccard = (denom == 0) ? 1.0 : (double)common / denom;
                
                // Compare Jaccard directly (larger Jaccard = smaller distance)
                // Update thread-local best (no array access!)
                if(jaccard > local_best_jaccard) {
                    local_best_jaccard = jaccard;
                    local_best_rep = repId;
                }
            }
            
            // Store thread-local best to shared array (only once per thread!)
            thread_best_matches[tid].jaccard = local_best_jaccard;
            thread_best_matches[tid].rep_id = local_best_rep;
            
            // Aggregate thread-local counters (only 2 atomic ops per thread, not per candidate!)
            #pragma omp atomic
            total_comparisons += local_comparisons;
            #pragma omp atomic
            filtered_by_common += local_filtered_common;
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
            int cid = rep2cid[best_rep];  // O(1) lookup
            cluster_members[cid].push_back(j);  // O(1) append
        } else {
            // Becomes new representative
            int new_cid = representativeArr.size();  // Next cluster ID
            representativeArr.push_back(j);
            rep2cid[j] = new_cid;  // O(1) assignment
            cluster_members.emplace_back();  // O(1) create empty member list
            
            // Critical: add new representative to inverted index
            if(use64)
                dynamic_index.add_representative(j, sketches[j].hash64_arr);
            else
                dynamic_index.add_representative(j, sketches[j].hash32_arr);
        }
        
        // Progress report (reduced frequency to minimize stderr overhead)
        if(j % 10000 == 0 || j == numGenomes - 1) {
            double clustering_rate = 100.0 * (j - representativeArr.size() + 1) / j;
            uint64_t total_filtered = filtered_by_common;
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
    cluster.reserve(representativeArr.size());
    
    for(size_t cid = 0; cid < representativeArr.size(); cid++) {
        std::vector<int> curClust;
        int center = representativeArr[cid];
        const auto& members = cluster_members[cid];

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
    std::cerr << "Filtered by common threshold: " << filtered_by_common << std::endl;
    std::cerr << "Skipped (no intersection): " << skipped_no_intersection << std::endl;
    
    uint64_t total_candidates = skipped_no_intersection + filtered_by_common + total_comparisons;
    uint64_t total_skipped = skipped_no_intersection + filtered_by_common;
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


KssdClusterState KssdInitialClusterWithState(
    vector<KssdSketchInfo>& sketches,
    const KssdParameters& params,
    double threshold,
    int threads,
    int kmer_size)
{
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "Initial Clustering with State Saving" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    

    vector<vector<int>> clusters = KssdGreedyClusterWithInvertedIndex(
        sketches, 0, threshold, threads, kmer_size);
    

    KssdClusterState state;
    state.all_sketches = sketches;
    state.clusters = clusters;
    state.params = params;
    state.threshold = threshold;
    state.kmer_size = kmer_size;
    

    state.representative_ids.reserve(clusters.size());
    state.representatives.reserve(clusters.size());
    
    for (const auto& cluster : clusters) {
        if (!cluster.empty()) {
            int rep_id = cluster[0];  // 簇的第一个成员是代表
            state.representative_ids.push_back(rep_id);
            state.representatives.push_back(sketches[rep_id]);
        }
    }
    

    std::cerr << "\nBuilding inverted index from " << state.representatives.size() 
             << " representatives..." << std::endl;
    
    state.inverted_index.clear();
    state.inverted_index.reserve(10000000); 
    
    for (size_t rep_idx = 0; rep_idx < state.representatives.size(); rep_idx++) {
        const auto& rep = state.representatives[rep_idx];
        if(rep.use64) {
            for(uint64_t h : rep.hash64_arr) {
                state.inverted_index[h].push_back(static_cast<int>(rep_idx));
            }
        } else {
            for(uint32_t h : rep.hash32_arr) {
                state.inverted_index[static_cast<uint64_t>(h)].push_back(static_cast<int>(rep_idx));
            }
        }
    }
    
    std::cerr << "Inverted index built: " << state.inverted_index.size() << " unique hashes" << std::endl;
    std::cerr << "State ready for incremental updates!" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    return state;
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
    
    // Optimized cluster storage: use array + vector<vector<int>> instead of unordered_map
    std::vector<int> rep2cid(numGenomes, -1);
    std::vector<int> representativeArr;
    std::vector<std::vector<int>> cluster_members;

    std::cerr << "========================================\n" << std::endl;

    // Dynamic inverted index for MinHash (using uint64_t hashes)
    class MinHashInvertedIndex {
    public:
        // Inverted index: hash -> list of representative IDs containing this hash
        // Use phmap (Google's Swiss Tables) for 3x faster lookups compared to std::unordered_map
        // Made public for direct access in array-based counting optimization
        phmap::flat_hash_map<uint64_t, std::vector<uint32_t>> index_map;

    private:
        // Set of representative sequences for fast lookup
        phmap::flat_hash_set<int> representative_set;

    public:
        MinHashInvertedIndex() {
            // Pre-allocate to avoid rehashing during construction
            // For large datasets (e.g., 2.7M genomes), unique hashes can be millions
            // Estimate: ~50K reps × ~1000 hashes/sketch × ~20% uniqueness = ~10M unique hashes
            // Conservative reserve to avoid expensive rehashing
            index_map.reserve(10000000);  // 10M unique hashes (adjust based on your dataset)
            representative_set.reserve(50000);  // Expected representatives (~2% of genomes)
        }

        // Add a new representative to the inverted index
        void add_representative(int rep_id, const std::vector<uint64_t>& hash_array) {
            representative_set.insert(rep_id);

            // Insert all hashes of this representative into the index
            // Optimization: use try_emplace to avoid double lookup (operator[] does find + insert)
            for(uint64_t hash : hash_array) {
                auto result = index_map.try_emplace(hash);
            auto it = result.first;
            // bool inserted = result.second;  // unused
                it->second.push_back(rep_id);
            }
        }

        // Check if a sequence is a representative
        bool is_representative(int id) const {
            return representative_set.count(id) > 0;
        }
        
        // Note: MinHash standard mode uses fixed sketch size, so size-based pruning is not applicable

        size_t num_representatives() const {
            return representative_set.size();
        }

        const phmap::flat_hash_set<int>& get_representatives() const {
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
    rep2cid[0] = 0;
    cluster_members.emplace_back();
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
    uint64_t filtered_by_common = 0;  // Filtered by insufficient common hashes
    uint64_t skipped_no_intersection = 0;  // Skipped representatives with no intersection

    // ===== Main loop: serial processing of each new sequence =====
    for(int j = 1; j < numGenomes; j++) {
        // Reuse query_hashes vector (still need to call storeMinHashes unfortunately)
        query_hashes = sketches[j].minHash->storeMinHashes();
        int sizeRef = query_hashes.size();
        bool isContainment = sketches[j].isContainment;

        // Note: MinHash standard mode uses fixed sketch size, no pruning needed
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
        size_t num_total_reps = dynamic_index.num_representatives();
        skipped_no_intersection += (num_total_reps - num_candidates);
        
        // No need to convert to candidates vector - touched is already a vector!
        // Just read common from cnt[rep_id] directly

        // Thread-local counters to avoid atomic contention in hot loop
        #pragma omp parallel num_threads(threads)
        {
            // Get thread ID once per parallel region (not per iteration!)
            int tid = omp_get_thread_num();
            
            uint64_t local_comparisons = 0;
            uint64_t local_filtered_common = 0;
            
            // Thread-local best match (avoids array access in hot loop!)
            int local_best_common = -1;
            double local_best_dist = std::numeric_limits<double>::max();
            int local_best_rep = -1;
            
            #pragma omp for schedule(dynamic, 64)
            for(size_t i = 0; i < touched.size(); i++) {
                int repId = touched[i];
                int common = cnt[repId];  // Read from array, not from pair
                int sizeQry = sketches[repId].minHash->getSketchSize();
                bool repIsContainment = sketches[repId].isContainment;

                // Optimization: Common threshold filtering
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
                    local_filtered_common++;
                    continue;
                }

                local_comparisons++;

                // Best-match selection strategy:
                // For standard MinHash with fixed size: compare common directly (faster!)
                // For containment: still need to compute distance
                
                if(all_fixed_size && all_standard_mode && !repIsContainment && !isContainment) {
                    // Fast path: For fixed-size standard MinHash, larger common => larger Jaccard => smaller distance
                    // So just compare common directly, no need to compute distance!
                    if(common > local_best_common) {
                        local_best_common = common;
                        local_best_rep = repId;
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
                    
                    if(dist <= threshold && dist < local_best_dist) {
                        local_best_dist = dist;
                        local_best_common = common;
                        local_best_rep = repId;
                    }
                }
            }  // End of for loop
            
            // Store thread-local best to shared array (only once per thread!)
            thread_best_matches[tid].common = local_best_common;
            thread_best_matches[tid].distance = local_best_dist;
            thread_best_matches[tid].rep_id = local_best_rep;
            
            // Aggregate thread-local counters (only 2 atomic ops per thread, not per candidate!)
            #pragma omp atomic
            total_comparisons += local_comparisons;
            #pragma omp atomic
            filtered_by_common += local_filtered_common;
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
            int cid = rep2cid[best_rep];
            cluster_members[cid].push_back(j);
        } else {
            // Becomes new representative
            int new_cid = representativeArr.size();
            representativeArr.push_back(j);
            rep2cid[j] = new_cid;
            cluster_members.emplace_back();
            
            // Critical: add new representative to inverted index
            std::vector<uint64_t> new_rep_hashes = sketches[j].minHash->storeMinHashes();
            dynamic_index.add_representative(j, new_rep_hashes);
        }
        
        // Progress report (reduced frequency to minimize stderr overhead)
        if(j % 10000 == 0 || j == numGenomes - 1) {
            double clustering_rate = 100.0 * (j - representativeArr.size() + 1) / j;
            uint64_t total_filtered = filtered_by_common;
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
    cluster.reserve(representativeArr.size());
    
    for(size_t cid = 0; cid < representativeArr.size(); cid++) {
        int center = representativeArr[cid];
        const auto& members = cluster_members[cid];
        
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
    std::cerr << "Filtered by common threshold: " << filtered_by_common << std::endl;
    std::cerr << "Skipped (no intersection): " << skipped_no_intersection << std::endl;
    
    uint64_t total_candidates = skipped_no_intersection + filtered_by_common + total_comparisons;
    uint64_t total_skipped = skipped_no_intersection + filtered_by_common;
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
    const bool use64 = sketches[0].use64;
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "BATCHED KSSD Greedy Clustering" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Total genomes: " << numGenomes << std::endl;
    std::cerr << "Batch size: " << batch_size << std::endl;
    std::cerr << "Threshold: " << threshold << std::endl;
    std::cerr << "Threads: " << threads << std::endl;
    std::cerr << "KSSD hash storage: " << (use64 ? "64-bit" : "32-bit") << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Initialize
    std::vector<int> clustLabels(numGenomes, 0);
    std::unordered_map<int, std::vector<int>> semiClust;
    DynamicInvertedIndex dynamic_index;
    
    // First sequence as initial representative
    semiClust[0] = std::vector<int>();
    if(use64)
        dynamic_index.add_representative(0, sketches[0].hash64_arr);
    else
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
            
            int sizeRef = (int)kssd_hash_count(sketches[j]);
            // Optimization: use sparse map, only store representatives, completely skip non-representatives
            std::unordered_map<int, int> intersection_map;
            
            // Compute intersections (note: each thread computes independently, with duplication)
            if(use64)
                dynamic_index.calculate_intersections_sparse(sketches[j].hash64_arr, intersection_map);
            else
                dynamic_index.calculate_intersections_sparse(sketches[j].hash32_arr, intersection_map);
            
            // Find best representative (directly traverse map, more efficient)
            for(const auto& pair : intersection_map) {
                int repId = pair.first;
                int common = pair.second;
                int sizeQry = (int)kssd_hash_count(sketches[repId]);
                
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
                if(use64)
                    dynamic_index.add_representative(j, sketches[j].hash64_arr);
                else
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
    for(const auto& pair : semiClust) {
        int center = pair.first;
        const std::vector<int>& members = pair.second;
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


bool KssdClusterState::save(const string& filepath) const {
    std::ofstream ofs(filepath, std::ios::binary);
    if (!ofs) {
        std::cerr << "ERROR: Cannot open file for writing: " << filepath << std::endl;
        return false;
    }
    

    ofs.write((char*)&threshold, sizeof(double));
    ofs.write((char*)&kmer_size, sizeof(int));
    ofs.write((char*)&params.half_k, sizeof(int));
    ofs.write((char*)&params.half_subk, sizeof(int));
    ofs.write((char*)&params.drlevel, sizeof(int));
    ofs.write((char*)&params.genomeNumber, sizeof(int));
    

    size_t rep_count = representative_ids.size();
    ofs.write((char*)&rep_count, sizeof(size_t));
    ofs.write((char*)representative_ids.data(), sizeof(int) * rep_count);
    

    size_t sketch_count = all_sketches.size();
    ofs.write((char*)&sketch_count, sizeof(size_t));
    for (const auto& sketch : all_sketches) {

        ofs.write((char*)&sketch.id, sizeof(int));
        ofs.write((char*)&sketch.totalSeqLength, sizeof(uint64_t));
        ofs.write((char*)&sketch.use64, sizeof(bool));
        ofs.write((char*)&sketch.sketchsize, sizeof(uint32_t));
        

        size_t hash32_size = sketch.hash32_arr.size();
        size_t hash64_size = sketch.hash64_arr.size();
        ofs.write((char*)&hash32_size, sizeof(size_t));
        ofs.write((char*)&hash64_size, sizeof(size_t));
        if (hash32_size > 0) {
            ofs.write((char*)sketch.hash32_arr.data(), sizeof(uint32_t) * hash32_size);
        }
        if (hash64_size > 0) {
            ofs.write((char*)sketch.hash64_arr.data(), sizeof(uint64_t) * hash64_size);
        }
        

        size_t name_len = sketch.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sketch.fileName.c_str(), name_len);
    }

    size_t cluster_count = clusters.size();
    ofs.write((char*)&cluster_count, sizeof(size_t));
    for (const auto& cluster : clusters) {
        size_t member_count = cluster.size();
        ofs.write((char*)&member_count, sizeof(size_t));
        ofs.write((char*)cluster.data(), sizeof(int) * member_count);
    }
    

    // Inverted index: marker + 64-bit keys (supports 32- and 64-bit KSSD sketches)
    const char kssd_inv_magic[8] = {'K','S','S','I','0','2','\0','\0'};
    ofs.write(kssd_inv_magic, 8);
    size_t index_size = inverted_index.size();
    ofs.write((char*)&index_size, sizeof(size_t));
    std::cerr << "Saving inverted index: " << index_size << " unique hashes..." << std::endl;
    
    for (const auto& pair : inverted_index) {
        uint64_t hash = pair.first;
        const std::vector<int>& rep_list = pair.second;
        ofs.write((char*)&hash, sizeof(uint64_t));
        size_t list_size = rep_list.size();
        ofs.write((char*)&list_size, sizeof(size_t));
        ofs.write((char*)rep_list.data(), sizeof(int) * list_size);
    }
    
    ofs.close();
    std::cerr << "Saved clustering state to: " << filepath << std::endl;
    std::cerr << "  - " << sketch_count << " genomes" << std::endl;
    std::cerr << "  - " << rep_count << " clusters (representatives)" << std::endl;
    std::cerr << "  - " << index_size << " unique hashes in inverted index" << std::endl;
    return true;
}


bool KssdClusterState::load(const string& filepath) {
    std::ifstream ifs(filepath, std::ios::binary);
    if (!ifs) {
        std::cerr << "ERROR: Cannot open file for reading: " << filepath << std::endl;
        return false;
    }
    

    ifs.read((char*)&threshold, sizeof(double));
    ifs.read((char*)&kmer_size, sizeof(int));
    ifs.read((char*)&params.half_k, sizeof(int));
    ifs.read((char*)&params.half_subk, sizeof(int));
    ifs.read((char*)&params.drlevel, sizeof(int));
    ifs.read((char*)&params.genomeNumber, sizeof(int));
    

    size_t rep_count;
    ifs.read((char*)&rep_count, sizeof(size_t));
    representative_ids.resize(rep_count);
    ifs.read((char*)representative_ids.data(), sizeof(int) * rep_count);
    

    size_t sketch_count;
    ifs.read((char*)&sketch_count, sizeof(size_t));
    all_sketches.resize(sketch_count);
    for (auto& sketch : all_sketches) {

        ifs.read((char*)&sketch.id, sizeof(int));
        ifs.read((char*)&sketch.totalSeqLength, sizeof(uint64_t));
        ifs.read((char*)&sketch.use64, sizeof(bool));
        ifs.read((char*)&sketch.sketchsize, sizeof(uint32_t));
        

        size_t hash32_size, hash64_size;
        ifs.read((char*)&hash32_size, sizeof(size_t));
        ifs.read((char*)&hash64_size, sizeof(size_t));
        if (hash32_size > 0) {
            sketch.hash32_arr.resize(hash32_size);
            ifs.read((char*)sketch.hash32_arr.data(), sizeof(uint32_t) * hash32_size);
        }
        if (hash64_size > 0) {
            sketch.hash64_arr.resize(hash64_size);
            ifs.read((char*)sketch.hash64_arr.data(), sizeof(uint64_t) * hash64_size);
        }
        

        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        sketch.fileName.resize(name_len);
        ifs.read(&sketch.fileName[0], name_len);
    }
    

    representatives.clear();
    for (int rep_id : representative_ids) {
        representatives.push_back(all_sketches[rep_id]);
    }
    
    size_t cluster_count;
    ifs.read((char*)&cluster_count, sizeof(size_t));
    clusters.resize(cluster_count);
    for (auto& cluster : clusters) {
        size_t member_count;
        ifs.read((char*)&member_count, sizeof(size_t));
        cluster.resize(member_count);
        ifs.read((char*)cluster.data(), sizeof(int) * member_count);
    }
    

    std::streampos inv_start = ifs.tellg();
    char inv_probe[8];
    ifs.read(inv_probe, 8);
    bool inv64 = (std::strncmp(inv_probe, "KSSI02", 6) == 0 && inv_probe[6] == '\0');
    if(!inv64) {
        ifs.seekg(inv_start);
    }
    size_t index_size;
    ifs.read((char*)&index_size, sizeof(size_t));
    std::cerr << "Loading inverted index: " << index_size << " unique hashes (" 
              << (inv64 ? "64-bit keys" : "legacy 32-bit keys") << ")..." << std::endl;
    
    inverted_index.clear();
    inverted_index.reserve(index_size);
    
    for (size_t i = 0; i < index_size; i++) {
        uint64_t hash64;
        if(inv64) {
            ifs.read((char*)&hash64, sizeof(uint64_t));
        } else {
            uint32_t hash32;
            ifs.read((char*)&hash32, sizeof(uint32_t));
            hash64 = static_cast<uint64_t>(hash32);
        }
        size_t list_size;
        ifs.read((char*)&list_size, sizeof(size_t));
        vector<int> rep_list(list_size);
        ifs.read((char*)rep_list.data(), sizeof(int) * list_size);
        inverted_index[hash64] = std::move(rep_list);
    }
    
    ifs.close();
    std::cerr << "Loaded clustering state from: " << filepath << std::endl;
    std::cerr << "  - " << sketch_count << " genomes" << std::endl;
    std::cerr << "  - " << rep_count << " clusters (representatives)" << std::endl;
    std::cerr << "  - " << index_size << " unique hashes in inverted index" << std::endl;
    return true;
}


vector<vector<int>> KssdIncrementalCluster(
    KssdClusterState& state,
    vector<KssdSketchInfo>& new_sketches,
    int threads)
{
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "Incremental Clustering with Inverted Index" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Existing clusters: " << state.representatives.size() << std::endl;
    std::cerr << "New genomes: " << new_sketches.size() << std::endl;
    std::cerr << "Inverted index size: " << state.inverted_index.size() << " unique hashes" << std::endl;
    
    int old_genome_count = state.all_sketches.size();
    double radio = calculateMaxSizeRatio(state.threshold, state.kmer_size);
    
    // Pre-calculate minimum jaccard threshold
    double x = std::exp(-state.threshold * state.kmer_size);
    double jaccard_min = x / (2.0 - x);
    

    int new_genome_start_idx = old_genome_count;
    for (auto& new_sketch : new_sketches) {
        state.all_sketches.push_back(new_sketch);
    }
    

    int new_clusters = 0;
    int assigned_to_existing = 0;
    uint64_t total_candidates = 0;
    uint64_t total_distance_calcs = 0;
    
    for (size_t new_idx = 0; new_idx < new_sketches.size(); new_idx++) {
        int genome_idx = new_genome_start_idx + new_idx;  // 在 all_sketches 中的索引
        auto& new_sketch = state.all_sketches[genome_idx];
        int sizeQry = (int)kssd_hash_count(new_sketch);
        
        phmap::flat_hash_map<int, int> candidate_counts;  // rep_idx → intersection count
        candidate_counts.reserve(1000); 
        
        if(new_sketch.use64) {
            for(uint64_t h : new_sketch.hash64_arr) {
                auto it = state.inverted_index.find(h);
                if(it != state.inverted_index.end()) {
                    for(int rep_idx : it->second) {
                        candidate_counts[rep_idx]++;
                    }
                }
            }
        } else {
            for(uint32_t h : new_sketch.hash32_arr) {
                auto it = state.inverted_index.find(static_cast<uint64_t>(h));
                if(it != state.inverted_index.end()) {
                    for(int rep_idx : it->second) {
                        candidate_counts[rep_idx]++;
                    }
                }
            }
        }
        
        total_candidates += candidate_counts.size();
        

        double best_dist = std::numeric_limits<double>::max();
        int best_rep_idx = -1;
        

        std::vector<std::pair<int, int>> candidates;  // (rep_idx, intersection_count)
        candidates.reserve(candidate_counts.size());
        for (const auto& pair : candidate_counts) {
            int rep_idx = pair.first;
            int cnt = pair.second;
            candidates.emplace_back(rep_idx, cnt);
        }
        
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for (size_t i = 0; i < candidates.size(); i++) {
            int rep_idx = candidates[i].first;
            int common = candidates[i].second;

            if (rep_idx < 0 || rep_idx >= (int)state.representatives.size()) continue;

            const auto& rep = state.representatives[rep_idx];
            int sizeRef = (int)kssd_hash_count(rep);
            
            // Size ratio filtering
            double ratio = (double)sizeQry / sizeRef;
            if (ratio > radio || ratio < 1.0 / radio) {
                continue;
            }
            
            // Common filtering: min_jaccard = common / (sizeQry + sizeRef - common)
            // => common >= jaccard_min * (sizeQry + sizeRef) / (1 + jaccard_min)
            int min_common_needed = (int)(jaccard_min * (sizeQry + sizeRef) / (1.0 + jaccard_min));
            if (common < min_common_needed) {
                continue;
            }
            
            // Calculate precise distance
            #pragma omp atomic
            total_distance_calcs++;
            
            double dist;
            if(rep.use64)
                dist = calculate_mash_distance(rep.hash64_arr, new_sketch.hash64_arr, state.kmer_size);
            else
                dist = calculate_mash_distance(rep.hash32_arr, new_sketch.hash32_arr, state.kmer_size);
            
            if (dist <= state.threshold) {
                #pragma omp critical
                {
                    if (dist < best_dist) {
                        best_dist = dist;
                        best_rep_idx = rep_idx;
                    }
                }
            }
        }
        
        
        if (best_rep_idx != -1) {

            state.clusters[best_rep_idx].push_back(genome_idx);
            assigned_to_existing++;
        } else {
  
            int new_rep_idx = state.representatives.size();
            state.representative_ids.push_back(genome_idx);
            state.representatives.push_back(new_sketch);
            state.clusters.push_back(vector<int>()); 
            new_clusters++;
            

            if(new_sketch.use64) {
                for(uint64_t h : new_sketch.hash64_arr) {
                    state.inverted_index[h].push_back(new_rep_idx);
                }
            } else {
                for(uint32_t h : new_sketch.hash32_arr) {
                    state.inverted_index[static_cast<uint64_t>(h)].push_back(new_rep_idx);
                }
            }
        }
        
        if ((new_idx + 1) % 10000 == 0) {
            std::cerr << "---finished clustering: " << (new_idx + 1) << " new genomes" << std::endl;
        }
    }
    
    std::cerr << "\n===== Incremental Clustering Results =====" << std::endl;
    std::cerr << "Assigned to existing clusters: " << assigned_to_existing << std::endl;
    std::cerr << "New clusters created: " << new_clusters << std::endl;
    std::cerr << "Total clusters now: " << state.clusters.size() << std::endl;
    if (!new_sketches.empty()) {
        std::cerr << "Average candidates per query: " << (total_candidates / (double)new_sketches.size()) << std::endl;
    }
    std::cerr << "Total distance calculations: " << total_distance_calcs << std::endl;
    if (total_distance_calcs > 0) {
        std::cerr << "Speedup ratio: " << (state.representatives.size() * new_sketches.size() / (double)total_distance_calcs) << "x" << std::endl;
    } else {
        std::cerr << "Speedup ratio: N/A (no distance calculations performed)" << std::endl;
    }
    std::cerr << "=========================================\n" << std::endl;
    
    return state.clusters;
}

// ===== MinHash Incremental Clustering Implementation =====

MinHashClusterState MinHashInitialClusterWithState(
    vector<SketchInfo>& sketches,
    double threshold,
    int threads,
    int kmer_size,
    int sketch_size,
    bool is_containment)
{
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "MinHash Initial Clustering with State Saving" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    // Call existing clustering function
    vector<vector<int>> clusters = MinHashGreedyClusterWithInvertedIndex(
        sketches, 0, threshold, threads, kmer_size);
    
    // Build complete state
    MinHashClusterState state;
    state.all_sketches = sketches;
    state.clusters = clusters;
    state.threshold = threshold;
    state.kmer_size = kmer_size;
    state.sketch_size = sketch_size;
    state.is_containment = is_containment;
    
    // Extract representatives
    state.representative_ids.reserve(clusters.size());
    state.representatives.reserve(clusters.size());
    
    for (const auto& cluster : clusters) {
        if (!cluster.empty()) {
            int rep_id = cluster[0];  // First member is representative
            state.representative_ids.push_back(rep_id);
            state.representatives.push_back(sketches[rep_id]);
        }
    }
    
    // Build inverted index from representatives
    std::cerr << "\nBuilding inverted index from " << state.representatives.size() 
             << " representatives..." << std::endl;
    
    state.inverted_index.clear();
    state.inverted_index.reserve(10000000);  // Estimate: ~50K reps × ~1000 hashes
    
    for (size_t rep_idx = 0; rep_idx < state.representatives.size(); rep_idx++) {
        const auto& rep = state.representatives[rep_idx];
        std::vector<uint64_t> hashes = rep.minHash->storeMinHashes();
        for (uint64_t hash : hashes) {
            state.inverted_index[hash].push_back(rep_idx);
        }
    }
    
    std::cerr << "Inverted index built: " << state.inverted_index.size() << " unique hashes" << std::endl;
    std::cerr << "State ready for incremental updates!" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    return state;
}

vector<vector<int>> MinHashIncrementalCluster(
    MinHashClusterState& state,
    vector<SketchInfo>& new_sketches,
    int threads)
{
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "MinHash Incremental Clustering with Inverted Index" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Existing clusters: " << state.representatives.size() << std::endl;
    std::cerr << "New genomes: " << new_sketches.size() << std::endl;
    std::cerr << "Inverted index size: " << state.inverted_index.size() << " unique hashes" << std::endl;
    
    if (new_sketches.empty()) {
        std::cerr << "ERROR: No new sketches to process" << std::endl;
        return state.clusters;
    }
    
    if (state.representatives.empty()) {
        std::cerr << "ERROR: No existing representatives" << std::endl;
        return state.clusters;
    }
    
    int old_genome_count = state.all_sketches.size();
    
    int new_genome_start_idx = old_genome_count;
    for (auto& new_sketch : new_sketches) {
        state.all_sketches.push_back(new_sketch);
    }

    if (state.rep_hash_arrays.empty()) {
        state.rep_hash_arrays.resize(state.representatives.size());
        for (size_t i = 0; i < state.representatives.size(); i++) {
            if (state.representatives[i].minHash)
                state.rep_hash_arrays[i] = state.representatives[i].minHash->storeMinHashes();
        }
    }

    double x = std::exp(-state.threshold * state.kmer_size);
    double jaccard_min = x / (2.0 - x);
    
    int new_clusters = 0;
    int assigned_to_existing = 0;
    uint64_t total_candidates = 0;
    uint64_t total_distance_calcs = 0;
    
    for (size_t new_idx = 0; new_idx < new_sketches.size(); new_idx++) {
        int genome_idx = new_genome_start_idx + new_idx;
        auto& new_sketch = state.all_sketches[genome_idx];
        
        if (!new_sketch.minHash) {
            std::cerr << "ERROR: New sketch at index " << new_idx << " has null minHash, skipping" << std::endl;
            continue;
        }
        
        std::vector<uint64_t> new_hashes = new_sketch.minHash->storeMinHashes();
        int sizeQry = new_hashes.size();
        
        phmap::flat_hash_map<int, int> candidate_counts;
        candidate_counts.reserve(1000);
        
        for (uint64_t hash : new_hashes) {
            auto it = state.inverted_index.find(hash);
            if (it != state.inverted_index.end()) {
                for (int rep_idx : it->second) {
                    candidate_counts[rep_idx]++;
                }
            }
        }
        
        total_candidates += candidate_counts.size();
        
        double best_dist = std::numeric_limits<double>::max();
        int best_rep_idx = -1;
        
        std::vector<std::pair<int, int>> candidates;
        candidates.reserve(candidate_counts.size());
        for (const auto& pair : candidate_counts)
            candidates.emplace_back(pair.first, pair.second);
        
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for (size_t i = 0; i < candidates.size(); i++) {
            int rep_idx = candidates[i].first;
            int common = candidates[i].second;
            
            if (rep_idx < 0 || rep_idx >= (int)state.rep_hash_arrays.size()) continue;
            
            int sizeRef = state.rep_hash_arrays[rep_idx].size();
            if (sizeRef == 0) continue;

            int min_common_needed;
            if (state.is_containment) {
                int minSz = std::min(sizeQry, sizeRef);
                min_common_needed = (int)(jaccard_min * minSz);
            } else {
                min_common_needed = (int)(jaccard_min * (sizeQry + sizeRef) / (1.0 + jaccard_min));
            }
            if (common < min_common_needed) continue;

            #pragma omp atomic
            total_distance_calcs++;

            double jaccard;
            if (state.is_containment) {
                int minSz = std::min(sizeQry, sizeRef);
                jaccard = (double)common / minSz;
            } else {
                int denom = sizeQry + sizeRef - common;
                if (denom == 0) continue;
                jaccard = (double)common / denom;
            }

            double dist;
            if (jaccard >= 1.0) {
                dist = 0.0;
            } else if (jaccard <= 0.0) {
                dist = 1.0;
            } else {
                dist = -std::log(2.0 * jaccard / (1.0 + jaccard)) / state.kmer_size;
                if (dist > 1.0) dist = 1.0;
            }

            if (dist <= state.threshold) {
                #pragma omp critical
                {
                    if (dist < best_dist) {
                        best_dist = dist;
                        best_rep_idx = rep_idx;
                    }
                }
            }
        }
        
        if (best_rep_idx != -1) {
            state.clusters[best_rep_idx].push_back(genome_idx);
            assigned_to_existing++;
        } else {
            int new_rep_idx = state.representatives.size();
            state.representative_ids.push_back(genome_idx);
            state.representatives.push_back(new_sketch);
            state.rep_hash_arrays.push_back(new_hashes);
            state.clusters.push_back(vector<int>());
            new_clusters++;
            
            for (uint64_t hash : new_hashes) {
                state.inverted_index[hash].push_back(new_rep_idx);
            }
        }
        
        if ((new_idx + 1) % 10000 == 0) {
            std::cerr << "---finished clustering: " << (new_idx + 1) << " new genomes" << std::endl;
        }
    }
    
    std::cerr << "\n===== Incremental Clustering Results =====" << std::endl;
    std::cerr << "Assigned to existing clusters: " << assigned_to_existing << std::endl;
    std::cerr << "New clusters created: " << new_clusters << std::endl;
    std::cerr << "Total clusters now: " << state.clusters.size() << std::endl;
    if (new_sketches.size() > 0) {
        std::cerr << "Average candidates per query: " << (total_candidates / (double)new_sketches.size()) << std::endl;
    }
    std::cerr << "Total distance calculations: " << total_distance_calcs << std::endl;
    if (total_distance_calcs > 0) {
        std::cerr << "Speedup ratio: " << (state.representatives.size() * new_sketches.size() / (double)total_distance_calcs) << "x" << std::endl;
    } else {
        std::cerr << "Speedup ratio: N/A (no distance calculations performed)" << std::endl;
    }
    std::cerr << "=========================================\n" << std::endl;
    
    return state.clusters;
}

bool MinHashClusterState::save(const string& filepath) const {
    std::ofstream ofs(filepath, std::ios::binary);
    if (!ofs) {
        std::cerr << "ERROR: Cannot open file for writing: " << filepath << std::endl;
        return false;
    }
    
    // Write header
    char magic[8] = "MINHASH";
    ofs.write(magic, 8);
    
    // Write parameters
    ofs.write((char*)&threshold, sizeof(double));
    ofs.write((char*)&kmer_size, sizeof(int));
    ofs.write((char*)&sketch_size, sizeof(int));
    ofs.write((char*)&is_containment, sizeof(bool));
    
    // Write representative IDs
    size_t rep_count = representative_ids.size();
    ofs.write((char*)&rep_count, sizeof(size_t));
    ofs.write((char*)representative_ids.data(), sizeof(int) * rep_count);
    
    // Write all sketches
    size_t sketch_count = all_sketches.size();
    ofs.write((char*)&sketch_count, sizeof(size_t));
    for (const auto& sketch : all_sketches) {
        // Save basic info
        ofs.write((char*)&sketch.id, sizeof(int));
        ofs.write((char*)&sketch.totalSeqLength, sizeof(uint64_t));
        
        // Save MinHash hashes
        std::vector<uint64_t> hashes = sketch.minHash->storeMinHashes();
        size_t hash_count = hashes.size();
        ofs.write((char*)&hash_count, sizeof(size_t));
        if (hash_count > 0) {
            ofs.write((char*)hashes.data(), sizeof(uint64_t) * hash_count);
        }
        
        // Save file name
        size_t name_len = sketch.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sketch.fileName.c_str(), name_len);
    }
    
    // Write clusters
    size_t cluster_count = clusters.size();
    ofs.write((char*)&cluster_count, sizeof(size_t));
    for (const auto& cluster : clusters) {
        size_t member_count = cluster.size();
        ofs.write((char*)&member_count, sizeof(size_t));
        ofs.write((char*)cluster.data(), sizeof(int) * member_count);
    }
    
    // Write inverted index
    size_t index_size = inverted_index.size();
    ofs.write((char*)&index_size, sizeof(size_t));
    std::cerr << "Saving inverted index: " << index_size << " unique hashes..." << std::endl;
    
    for (const auto& pair : inverted_index) {
        uint64_t hash = pair.first;
        const std::vector<int>& rep_list = pair.second;
        ofs.write((char*)&hash, sizeof(uint64_t));
        size_t list_size = rep_list.size();
        ofs.write((char*)&list_size, sizeof(size_t));
        ofs.write((char*)rep_list.data(), sizeof(int) * list_size);
    }
    
    ofs.close();
    std::cerr << "Saved clustering state to: " << filepath << std::endl;
    std::cerr << "  - " << sketch_count << " genomes" << std::endl;
    std::cerr << "  - " << rep_count << " clusters (representatives)" << std::endl;
    std::cerr << "  - " << index_size << " unique hashes in inverted index" << std::endl;
    return true;
}

bool MinHashClusterState::load(const string& filepath) {
    std::ifstream ifs(filepath, std::ios::binary);
    if (!ifs) {
        std::cerr << "ERROR: Cannot open file for reading: " << filepath << std::endl;
        return false;
    }
    
    // Read header
    char magic[8];
    ifs.read(magic, 8);
    if (strncmp(magic, "MINHASH", 8) != 0) {
        std::cerr << "ERROR: Invalid file format (not a MinHash cluster state)" << std::endl;
        return false;
    }
    
    // Read parameters
    ifs.read((char*)&threshold, sizeof(double));
    ifs.read((char*)&kmer_size, sizeof(int));
    ifs.read((char*)&sketch_size, sizeof(int));
    ifs.read((char*)&is_containment, sizeof(bool));
    
    // Read representative IDs
    size_t rep_count;
    ifs.read((char*)&rep_count, sizeof(size_t));
    representative_ids.resize(rep_count);
    ifs.read((char*)representative_ids.data(), sizeof(int) * rep_count);
    
    // Read all sketches (but we'll reload from files, so just skip the data)
    size_t sketch_count;
    ifs.read((char*)&sketch_count, sizeof(size_t));
    // Skip sketch data - we'll reload from files
    for (size_t i = 0; i < sketch_count; i++) {
        int id;
        uint64_t totalSeqLength;
        ifs.read((char*)&id, sizeof(int));
        ifs.read((char*)&totalSeqLength, sizeof(uint64_t));
        
        // Skip MinHash hashes
        size_t hash_count;
        ifs.read((char*)&hash_count, sizeof(size_t));
        if (hash_count > 0) {
            ifs.seekg(sizeof(uint64_t) * hash_count, std::ios::cur);
        }
        
        // Skip file name
        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        ifs.seekg(name_len, std::ios::cur);
    }
    
    // Representatives will be rebuilt from loaded sketches
    representatives.clear();
    
    // Read clusters
    size_t cluster_count;
    ifs.read((char*)&cluster_count, sizeof(size_t));
    clusters.resize(cluster_count);
    for (auto& cluster : clusters) {
        size_t member_count;
        ifs.read((char*)&member_count, sizeof(size_t));
        cluster.resize(member_count);
        ifs.read((char*)cluster.data(), sizeof(int) * member_count);
    }
    
    // Read inverted index
    size_t index_size;
    ifs.read((char*)&index_size, sizeof(size_t));
    std::cerr << "Loading inverted index: " << index_size << " unique hashes..." << std::endl;
    
    inverted_index.clear();
    inverted_index.reserve(index_size);
    
    for (size_t i = 0; i < index_size; i++) {
        uint64_t hash;
        ifs.read((char*)&hash, sizeof(uint64_t));
        
        size_t list_size;
        ifs.read((char*)&list_size, sizeof(size_t));
        
        vector<int> rep_list(list_size);
        ifs.read((char*)rep_list.data(), sizeof(int) * list_size);
        
        inverted_index[hash] = std::move(rep_list);
    }
    
    ifs.close();
    std::cerr << "Loaded clustering state from: " << filepath << std::endl;
    std::cerr << "  - " << sketch_count << " genomes" << std::endl;
    std::cerr << "  - " << rep_count << " clusters (representatives)" << std::endl;
    std::cerr << "  - " << index_size << " unique hashes in inverted index" << std::endl;
    return true;
}

void MinHashClusterState::build_inverted_index() {
    inverted_index.clear();
    inverted_index.reserve(10000000);
    
    for (size_t rep_idx = 0; rep_idx < representatives.size(); rep_idx++) {
        const auto& rep = representatives[rep_idx];
        if (!rep.minHash) {
            std::cerr << "ERROR: Representative " << rep_idx << " has null minHash, cannot build inverted index" << std::endl;
            exit(1);
        }
        try {
            std::vector<uint64_t> hashes = rep.minHash->storeMinHashes();
            if (hashes.empty()) {
                std::cerr << "WARNING: Representative " << rep_idx << " has empty hash array" << std::endl;
                continue;
            }
            for (uint64_t hash : hashes) {
                inverted_index[hash].push_back(rep_idx);
            }
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Exception while processing representative " << rep_idx << ": " << e.what() << std::endl;
            exit(1);
        } catch (...) {
            std::cerr << "ERROR: Unknown exception while processing representative " << rep_idx << std::endl;
            exit(1);
        }
    }
}

void MinHashClusterState::update_inverted_index(int rep_idx) {
    if (rep_idx < 0 || rep_idx >= representatives.size()) {
        std::cerr << "ERROR: update_inverted_index: rep_idx " << rep_idx << " out of range" << std::endl;
        return;
    }
    const auto& rep = representatives[rep_idx];
    if (!rep.minHash) {
        std::cerr << "ERROR: update_inverted_index: representative " << rep_idx << " has null minHash" << std::endl;
        return;
    }
    std::vector<uint64_t> hashes = rep.minHash->storeMinHashes();
    for (uint64_t hash : hashes) {
        inverted_index[hash].push_back(rep_idx);
    }
}

// =====================================================================
// RepDB: Representative Database methods for KssdClusterState
// =====================================================================

bool KssdClusterState::save_repdb(const string& filepath) const {
    std::ofstream ofs(filepath, std::ios::binary);
    if (!ofs) {
        std::cerr << "ERROR: Cannot open RepDB file for writing: " << filepath << std::endl;
        return false;
    }

    const char magic[] = "REPDB002";  // v2: inverted index uses 64-bit hash keys
    ofs.write(magic, 8);

    ofs.write((char*)&threshold, sizeof(double));
    ofs.write((char*)&kmer_size, sizeof(int));
    ofs.write((char*)&params.half_k, sizeof(int));
    ofs.write((char*)&params.half_subk, sizeof(int));
    ofs.write((char*)&params.drlevel, sizeof(int));
    ofs.write((char*)&params.genomeNumber, sizeof(int));

    size_t rep_count = representatives.size();
    ofs.write((char*)&rep_count, sizeof(size_t));

    for (size_t i = 0; i < rep_count; i++) {
        ofs.write((char*)&representative_ids[i], sizeof(int));

        const auto& sk = representatives[i];
        ofs.write((char*)&sk.id, sizeof(int));
        ofs.write((char*)&sk.totalSeqLength, sizeof(uint64_t));
        ofs.write((char*)&sk.use64, sizeof(bool));
        ofs.write((char*)&sk.sketchsize, sizeof(uint32_t));

        size_t h32sz = sk.hash32_arr.size();
        size_t h64sz = sk.hash64_arr.size();
        ofs.write((char*)&h32sz, sizeof(size_t));
        ofs.write((char*)&h64sz, sizeof(size_t));
        if (h32sz > 0)
            ofs.write((char*)sk.hash32_arr.data(), sizeof(uint32_t) * h32sz);
        if (h64sz > 0)
            ofs.write((char*)sk.hash64_arr.data(), sizeof(uint64_t) * h64sz);

        size_t name_len = sk.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sk.fileName.c_str(), name_len);
    }

    size_t cluster_count = clusters.size();
    ofs.write((char*)&cluster_count, sizeof(size_t));
    for (const auto& cl : clusters) {
        size_t member_count = cl.size();
        ofs.write((char*)&member_count, sizeof(size_t));
        ofs.write((char*)cl.data(), sizeof(int) * member_count);
    }

    size_t all_count = all_sketches.size();
    ofs.write((char*)&all_count, sizeof(size_t));
    for (const auto& sk : all_sketches) {
        size_t name_len = sk.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sk.fileName.c_str(), name_len);
        ofs.write((char*)&sk.totalSeqLength, sizeof(uint64_t));
    }

    size_t index_size = inverted_index.size();
    ofs.write((char*)&index_size, sizeof(size_t));
    for (const auto& pair : inverted_index) {
        uint64_t hash = pair.first;
        const std::vector<int>& rep_list = pair.second;
        ofs.write((char*)&hash, sizeof(uint64_t));
        size_t list_size = rep_list.size();
        ofs.write((char*)&list_size, sizeof(size_t));
        ofs.write((char*)rep_list.data(), sizeof(int) * list_size);
    }

    ofs.close();
    std::cerr << "RepDB saved to: " << filepath << std::endl;
    std::cerr << "  Representatives: " << rep_count << std::endl;
    std::cerr << "  Total genomes:   " << all_count << std::endl;
    std::cerr << "  Inverted index:  " << index_size << " unique hashes" << std::endl;
    return true;
}

bool KssdClusterState::load_repdb(const string& filepath) {
    std::ifstream ifs(filepath, std::ios::binary);
    if (!ifs) {
        std::cerr << "ERROR: Cannot open RepDB file for reading: " << filepath << std::endl;
        return false;
    }

    char magic[8];
    ifs.read(magic, 8);
    std::string mag(magic, 8);
    bool repdb_v2 = (mag == "REPDB002");
    if (!repdb_v2 && mag != "REPDB001") {
        std::cerr << "ERROR: Invalid RepDB file (bad magic): " << filepath << std::endl;
        return false;
    }

    ifs.read((char*)&threshold, sizeof(double));
    ifs.read((char*)&kmer_size, sizeof(int));
    ifs.read((char*)&params.half_k, sizeof(int));
    ifs.read((char*)&params.half_subk, sizeof(int));
    ifs.read((char*)&params.drlevel, sizeof(int));
    ifs.read((char*)&params.genomeNumber, sizeof(int));

    size_t rep_count;
    ifs.read((char*)&rep_count, sizeof(size_t));
    representative_ids.resize(rep_count);
    representatives.resize(rep_count);

    for (size_t i = 0; i < rep_count; i++) {
        ifs.read((char*)&representative_ids[i], sizeof(int));

        auto& sk = representatives[i];
        ifs.read((char*)&sk.id, sizeof(int));
        ifs.read((char*)&sk.totalSeqLength, sizeof(uint64_t));
        ifs.read((char*)&sk.use64, sizeof(bool));
        ifs.read((char*)&sk.sketchsize, sizeof(uint32_t));

        size_t h32sz, h64sz;
        ifs.read((char*)&h32sz, sizeof(size_t));
        ifs.read((char*)&h64sz, sizeof(size_t));
        if (h32sz > 0) {
            sk.hash32_arr.resize(h32sz);
            ifs.read((char*)sk.hash32_arr.data(), sizeof(uint32_t) * h32sz);
        }
        if (h64sz > 0) {
            sk.hash64_arr.resize(h64sz);
            ifs.read((char*)sk.hash64_arr.data(), sizeof(uint64_t) * h64sz);
        }

        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        sk.fileName.resize(name_len);
        ifs.read(&sk.fileName[0], name_len);
    }

    size_t cluster_count;
    ifs.read((char*)&cluster_count, sizeof(size_t));
    clusters.resize(cluster_count);
    for (auto& cl : clusters) {
        size_t member_count;
        ifs.read((char*)&member_count, sizeof(size_t));
        cl.resize(member_count);
        ifs.read((char*)cl.data(), sizeof(int) * member_count);
    }

    size_t all_count;
    ifs.read((char*)&all_count, sizeof(size_t));
    all_sketches.resize(all_count);
    for (auto& sk : all_sketches) {
        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        sk.fileName.resize(name_len);
        ifs.read(&sk.fileName[0], name_len);
        ifs.read((char*)&sk.totalSeqLength, sizeof(uint64_t));
    }

    size_t index_size;
    ifs.read((char*)&index_size, sizeof(size_t));
    inverted_index.clear();
    inverted_index.reserve(index_size);
    for (size_t i = 0; i < index_size; i++) {
        uint64_t hash64;
        if(repdb_v2) {
            ifs.read((char*)&hash64, sizeof(uint64_t));
        } else {
            uint32_t h32;
            ifs.read((char*)&h32, sizeof(uint32_t));
            hash64 = static_cast<uint64_t>(h32);
        }
        size_t list_size;
        ifs.read((char*)&list_size, sizeof(size_t));
        vector<int> rep_list(list_size);
        ifs.read((char*)rep_list.data(), sizeof(int) * list_size);
        inverted_index[hash64] = std::move(rep_list);
    }

    ifs.close();
    std::cerr << "RepDB loaded from: " << filepath << std::endl;
    std::cerr << "  Representatives: " << rep_count << std::endl;
    std::cerr << "  Total genomes:   " << all_count << std::endl;
    std::cerr << "  Inverted index:  " << index_size << " unique hashes" << std::endl;
    std::cerr << "  Threshold:       " << threshold << std::endl;
    std::cerr << "  Kmer size:       " << kmer_size << std::endl;
    return true;
}

vector<RepDBQueryResult> KssdClusterState::query_topk(
    const KssdSketchInfo& query, int topk, int threads) const
{
    double radio = calculateMaxSizeRatio(threshold, kmer_size);
    double x = std::exp(-threshold * kmer_size);
    double jaccard_min = x / (2.0 - x);
    int sizeQry = (int)kssd_hash_count(query);

    phmap::flat_hash_map<int, int> candidate_counts;
    candidate_counts.reserve(1000);

    if(query.use64) {
        for(uint64_t h : query.hash64_arr) {
            auto it = inverted_index.find(h);
            if(it != inverted_index.end()) {
                for(int rep_idx : it->second) {
                    candidate_counts[rep_idx]++;
                }
            }
        }
    } else {
        for(uint32_t h : query.hash32_arr) {
            auto it = inverted_index.find(static_cast<uint64_t>(h));
            if(it != inverted_index.end()) {
                for(int rep_idx : it->second) {
                    candidate_counts[rep_idx]++;
                }
            }
        }
    }

    struct DistPair {
        int rep_idx;
        double distance;
    };
    std::vector<DistPair> scored;
    scored.reserve(candidate_counts.size());

    std::vector<std::pair<int, int>> candidates;
    candidates.reserve(candidate_counts.size());
    for (const auto& p : candidate_counts) {
        candidates.emplace_back(p.first, p.second);
    }

    std::vector<DistPair> thread_results;
    std::mutex mtx;

    #pragma omp parallel num_threads(threads)
    {
        std::vector<DistPair> local;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < candidates.size(); i++) {
            int rep_idx = candidates[i].first;
            int common = candidates[i].second;
            const auto& rep = representatives[rep_idx];
            int sizeRef = (int)kssd_hash_count(rep);

            double ratio = (double)sizeQry / sizeRef;
            if (ratio > radio || ratio < 1.0 / radio)
                continue;

            int min_common = (int)(jaccard_min * (sizeQry + sizeRef) / (1.0 + jaccard_min));
            if (common < min_common)
                continue;

            double dist;
            if(rep.use64)
                dist = calculate_mash_distance(rep.hash64_arr, query.hash64_arr, kmer_size);
            else
                dist = calculate_mash_distance(rep.hash32_arr, query.hash32_arr, kmer_size);
            local.push_back({rep_idx, dist});
        }
        std::lock_guard<std::mutex> lock(mtx);
        scored.insert(scored.end(), local.begin(), local.end());
    }

    std::sort(scored.begin(), scored.end(),
        [](const DistPair& a, const DistPair& b) { return a.distance < b.distance; });

    int k = std::min(topk, (int)scored.size());
    vector<RepDBQueryResult> results;
    results.reserve(k);

    for (int i = 0; i < k; i++) {
        RepDBQueryResult r;
        r.rep_idx = scored[i].rep_idx;
        r.distance = scored[i].distance;

        if (r.rep_idx >= 0 && r.rep_idx < (int)representatives.size()) {
            r.genome_id = representative_ids[r.rep_idx];
            r.genome_name = representatives[r.rep_idx].fileName;
        } else {
            r.genome_id = -1;
            r.genome_name = "N/A";
        }
        r.cluster_id = r.rep_idx;
        r.cluster_size = (r.rep_idx >= 0 && r.rep_idx < (int)clusters.size())
                         ? (int)clusters[r.rep_idx].size() : 0;
        results.push_back(r);
    }
    return results;
}

RepDBQueryResult KssdClusterState::assign(
    const KssdSketchInfo& query, int threads) const
{
    auto results = query_topk(query, 1, threads);
    if (!results.empty() && results[0].distance <= threshold) {
        return results[0];
    }
    RepDBQueryResult none;
    none.rep_idx = -1;
    none.genome_id = -1;
    none.genome_name = "unassigned";
    none.distance = -1.0;
    none.cluster_id = -1;
    none.cluster_size = 0;
    return none;
}

void KssdClusterState::print_stats(std::ostream& out) const {
    size_t total_genomes = 0;
    for (const auto& cl : clusters) total_genomes += cl.size();

    out << "========================================" << std::endl;
    out << "        RepDB Statistics Report" << std::endl;
    out << "========================================" << std::endl;
    out << std::endl;

    out << "[Basic Info]" << std::endl;
    out << "  Threshold:              " << threshold << std::endl;
    out << "  Kmer size:              " << kmer_size << std::endl;
    out << "  KSSD half_k:            " << params.half_k << std::endl;
    out << "  KSSD half_subk:         " << params.half_subk << std::endl;
    out << "  KSSD drlevel:           " << params.drlevel << std::endl;
    out << std::endl;

    out << "[Scale]" << std::endl;
    out << "  Total genomes:          " << total_genomes << std::endl;
    out << "  Representatives:        " << representatives.size() << std::endl;
    out << "  Clusters:               " << clusters.size() << std::endl;
    double compression = total_genomes > 0
        ? (1.0 - (double)representatives.size() / total_genomes) * 100.0 : 0.0;
    out << "  Compression ratio:      " << std::fixed << std::setprecision(2)
        << compression << "%" << std::endl;
    out << std::endl;

    out << "[Inverted Index]" << std::endl;
    out << "  Unique hashes:          " << inverted_index.size() << std::endl;
    size_t total_postings = 0;
    size_t max_posting = 0;
    for (const auto& p : inverted_index) {
        total_postings += p.second.size();
        if (p.second.size() > max_posting) max_posting = p.second.size();
    }
    double avg_posting = inverted_index.size() > 0
        ? (double)total_postings / inverted_index.size() : 0;
    out << "  Total postings:         " << total_postings << std::endl;
    out << "  Avg posting length:     " << std::fixed << std::setprecision(2)
        << avg_posting << std::endl;
    out << "  Max posting length:     " << max_posting << std::endl;
    out << std::endl;

    out << "[Cluster Size Distribution]" << std::endl;
    if (!clusters.empty()) {
        std::vector<int> sizes;
        sizes.reserve(clusters.size());
        size_t singleton = 0;
        size_t max_sz = 0;
        for (const auto& cl : clusters) {
            sizes.push_back(cl.size());
            if (cl.size() <= 1) singleton++;
            if (cl.size() > max_sz) max_sz = cl.size();
        }
        std::sort(sizes.begin(), sizes.end());
        double mean = (double)total_genomes / clusters.size();
        int median = sizes[sizes.size() / 2];

        out << "  Min cluster size:       " << sizes.front() << std::endl;
        out << "  Max cluster size:       " << sizes.back() << std::endl;
        out << "  Mean cluster size:      " << std::fixed << std::setprecision(2)
            << mean << std::endl;
        out << "  Median cluster size:    " << median << std::endl;
        out << "  Singletons:             " << singleton << " ("
            << std::fixed << std::setprecision(1)
            << (100.0 * singleton / clusters.size()) << "%)" << std::endl;

        size_t p90 = sizes[(size_t)(sizes.size() * 0.9)];
        size_t p95 = sizes[(size_t)(sizes.size() * 0.95)];
        size_t p99 = sizes[(size_t)(sizes.size() * 0.99)];
        out << "  P90 cluster size:       " << p90 << std::endl;
        out << "  P95 cluster size:       " << p95 << std::endl;
        out << "  P99 cluster size:       " << p99 << std::endl;
    }
    out << std::endl;

    out << "[Representative Sketch Sizes]" << std::endl;
    if (!representatives.empty()) {
        size_t min_sk = UINT64_MAX, max_sk = 0, sum_sk = 0;
        for (const auto& r : representatives) {
            size_t sz = kssd_hash_count(r);
            if (sz < min_sk) min_sk = sz;
            if (sz > max_sk) max_sk = sz;
            sum_sk += sz;
        }
        out << "  Min sketch size:        " << min_sk << std::endl;
        out << "  Max sketch size:        " << max_sk << std::endl;
        out << "  Mean sketch size:       " << std::fixed << std::setprecision(1)
            << (double)sum_sk / representatives.size() << std::endl;
    }

    uint64_t total_seq_len = 0;
    for (const auto& sk : all_sketches) {
        total_seq_len += sk.totalSeqLength;
    }
    if (total_seq_len > 0) {
        out << std::endl;
        out << "[Genome Coverage]" << std::endl;
        uint64_t rep_seq_len = 0;
        for (const auto& r : representatives) {
            rep_seq_len += r.totalSeqLength;
        }
        out << "  Total sequence length:  " << total_seq_len << " bp" << std::endl;
        out << "  Representative seq len: " << rep_seq_len << " bp" << std::endl;
        out << "  Coverage ratio:         " << std::fixed << std::setprecision(2)
            << (100.0 * rep_seq_len / total_seq_len) << "%" << std::endl;
    }

    out << "========================================" << std::endl;
}

// =====================================================================
// RepDB: Representative Database methods for MinHashClusterState
// =====================================================================

static inline double minhash_mash_distance(int common, int sizeQry, int sizeRef, int kmer_size, bool is_containment) {
    if (common <= 0) return 1.0;
    double jaccard;
    if (is_containment) {
        int minSize = std::min(sizeQry, sizeRef);
        if (minSize == 0) return 1.0;
        jaccard = (double)common / minSize;
    } else {
        int denom = sizeQry + sizeRef - common;
        if (denom == 0) return 0.0;
        jaccard = (double)common / denom;
    }
    if (jaccard >= 1.0) return 0.0;
    if (jaccard <= 0.0) return 1.0;
    double dist = -std::log(2.0 * jaccard / (1.0 + jaccard)) / kmer_size;
    return dist > 1.0 ? 1.0 : dist;
}

bool MinHashClusterState::save_repdb(const string& filepath) const {
    std::ofstream ofs(filepath, std::ios::binary);
    if (!ofs) {
        std::cerr << "ERROR: Cannot open RepDB file for writing: " << filepath << std::endl;
        return false;
    }

    const char magic[] = "MHREPDB1";
    ofs.write(magic, 8);

    ofs.write((char*)&threshold, sizeof(double));
    ofs.write((char*)&kmer_size, sizeof(int));
    ofs.write((char*)&sketch_size, sizeof(int));
    ofs.write((char*)&is_containment, sizeof(bool));

    size_t rep_count = representatives.size();
    ofs.write((char*)&rep_count, sizeof(size_t));

    for (size_t i = 0; i < rep_count; i++) {
        ofs.write((char*)&representative_ids[i], sizeof(int));

        const auto& sk = representatives[i];
        ofs.write((char*)&sk.id, sizeof(int));
        ofs.write((char*)&sk.totalSeqLength, sizeof(uint64_t));
        ofs.write((char*)&sk.isContainment, sizeof(bool));

        const std::vector<uint64_t>& hashes = (i < rep_hash_arrays.size())
            ? rep_hash_arrays[i]
            : (sk.minHash ? sk.minHash->storeMinHashes() : std::vector<uint64_t>());
        size_t hash_count = hashes.size();
        ofs.write((char*)&hash_count, sizeof(size_t));
        if (hash_count > 0)
            ofs.write((char*)hashes.data(), sizeof(uint64_t) * hash_count);

        size_t name_len = sk.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sk.fileName.c_str(), name_len);
    }

    size_t cluster_count = clusters.size();
    ofs.write((char*)&cluster_count, sizeof(size_t));
    for (const auto& cl : clusters) {
        size_t member_count = cl.size();
        ofs.write((char*)&member_count, sizeof(size_t));
        ofs.write((char*)cl.data(), sizeof(int) * member_count);
    }

    size_t all_count = all_sketches.size();
    ofs.write((char*)&all_count, sizeof(size_t));
    for (const auto& sk : all_sketches) {
        size_t name_len = sk.fileName.length();
        ofs.write((char*)&name_len, sizeof(size_t));
        ofs.write(sk.fileName.c_str(), name_len);
        ofs.write((char*)&sk.totalSeqLength, sizeof(uint64_t));
    }

    size_t index_size = inverted_index.size();
    ofs.write((char*)&index_size, sizeof(size_t));
    for (const auto& pair : inverted_index) {
        uint64_t hash = pair.first;
        const std::vector<int>& rep_list = pair.second;
        ofs.write((char*)&hash, sizeof(uint64_t));
        size_t list_size = rep_list.size();
        ofs.write((char*)&list_size, sizeof(size_t));
        ofs.write((char*)rep_list.data(), sizeof(int) * list_size);
    }

    ofs.close();
    std::cerr << "MinHash RepDB saved to: " << filepath << std::endl;
    std::cerr << "  Representatives: " << rep_count << std::endl;
    std::cerr << "  Total genomes:   " << all_count << std::endl;
    std::cerr << "  Inverted index:  " << index_size << " unique hashes" << std::endl;
    return true;
}

bool MinHashClusterState::load_repdb(const string& filepath) {
    std::ifstream ifs(filepath, std::ios::binary);
    if (!ifs) {
        std::cerr << "ERROR: Cannot open RepDB file for reading: " << filepath << std::endl;
        return false;
    }

    char magic[8];
    ifs.read(magic, 8);
    if (std::string(magic, 8) != "MHREPDB1") {
        std::cerr << "ERROR: Invalid MinHash RepDB file (bad magic): " << filepath << std::endl;
        return false;
    }

    ifs.read((char*)&threshold, sizeof(double));
    ifs.read((char*)&kmer_size, sizeof(int));
    ifs.read((char*)&sketch_size, sizeof(int));
    ifs.read((char*)&is_containment, sizeof(bool));

    size_t rep_count;
    ifs.read((char*)&rep_count, sizeof(size_t));
    representative_ids.resize(rep_count);
    representatives.resize(rep_count);
    rep_hash_arrays.resize(rep_count);

    for (size_t i = 0; i < rep_count; i++) {
        ifs.read((char*)&representative_ids[i], sizeof(int));

        auto& sk = representatives[i];
        ifs.read((char*)&sk.id, sizeof(int));
        ifs.read((char*)&sk.totalSeqLength, sizeof(uint64_t));
        ifs.read((char*)&sk.isContainment, sizeof(bool));

        size_t hash_count;
        ifs.read((char*)&hash_count, sizeof(size_t));
        rep_hash_arrays[i].resize(hash_count);
        if (hash_count > 0)
            ifs.read((char*)rep_hash_arrays[i].data(), sizeof(uint64_t) * hash_count);

        sk.minHash = nullptr;

        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        sk.fileName.resize(name_len);
        ifs.read(&sk.fileName[0], name_len);
    }

    size_t cluster_count;
    ifs.read((char*)&cluster_count, sizeof(size_t));
    clusters.resize(cluster_count);
    for (auto& cl : clusters) {
        size_t member_count;
        ifs.read((char*)&member_count, sizeof(size_t));
        cl.resize(member_count);
        ifs.read((char*)cl.data(), sizeof(int) * member_count);
    }

    size_t all_count;
    ifs.read((char*)&all_count, sizeof(size_t));
    all_sketches.resize(all_count);
    for (auto& sk : all_sketches) {
        size_t name_len;
        ifs.read((char*)&name_len, sizeof(size_t));
        sk.fileName.resize(name_len);
        ifs.read(&sk.fileName[0], name_len);
        ifs.read((char*)&sk.totalSeqLength, sizeof(uint64_t));
    }

    size_t index_size;
    ifs.read((char*)&index_size, sizeof(size_t));
    inverted_index.clear();
    inverted_index.reserve(index_size);
    for (size_t i = 0; i < index_size; i++) {
        uint64_t hash;
        ifs.read((char*)&hash, sizeof(uint64_t));
        size_t list_size;
        ifs.read((char*)&list_size, sizeof(size_t));
        vector<int> rep_list(list_size);
        ifs.read((char*)rep_list.data(), sizeof(int) * list_size);
        inverted_index[hash] = std::move(rep_list);
    }

    ifs.close();
    std::cerr << "MinHash RepDB loaded from: " << filepath << std::endl;
    std::cerr << "  Representatives: " << rep_count << std::endl;
    std::cerr << "  Total genomes:   " << all_count << std::endl;
    std::cerr << "  Inverted index:  " << index_size << " unique hashes" << std::endl;
    std::cerr << "  Threshold:       " << threshold << std::endl;
    std::cerr << "  Kmer size:       " << kmer_size << std::endl;
    std::cerr << "  Sketch size:     " << sketch_size << std::endl;
    std::cerr << "  Containment:     " << (is_containment ? "yes" : "no") << std::endl;
    return true;
}

vector<RepDBQueryResult> MinHashClusterState::query_topk(
    const SketchInfo& query, int topk, int threads) const
{
    if (!query.minHash) {
        std::cerr << "ERROR: query has null minHash" << std::endl;
        return {};
    }

    std::vector<uint64_t> query_hashes = query.minHash->storeMinHashes();
    int sizeQry = query_hashes.size();

    phmap::flat_hash_map<int, int> candidate_counts;
    candidate_counts.reserve(1000);

    for (uint64_t hash : query_hashes) {
        auto it = inverted_index.find(hash);
        if (it != inverted_index.end()) {
            for (int rep_idx : it->second) {
                candidate_counts[rep_idx]++;
            }
        }
    }

    struct DistPair {
        int rep_idx;
        double distance;
    };
    std::vector<DistPair> scored;
    scored.reserve(candidate_counts.size());

    std::vector<std::pair<int, int>> candidates;
    candidates.reserve(candidate_counts.size());
    for (const auto& p : candidate_counts)
        candidates.emplace_back(p.first, p.second);

    std::mutex mtx;

    #pragma omp parallel num_threads(threads)
    {
        std::vector<DistPair> local;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < candidates.size(); i++) {
            int rep_idx = candidates[i].first;
            int common = candidates[i].second;
            if (rep_idx < 0 || rep_idx >= (int)rep_hash_arrays.size()) continue;

            int sizeRef = rep_hash_arrays[rep_idx].size();
            double dist = minhash_mash_distance(common, sizeQry, sizeRef, kmer_size, is_containment);

            if (!std::isnan(dist) && !std::isinf(dist))
                local.push_back({rep_idx, dist});
        }
        std::lock_guard<std::mutex> lock(mtx);
        scored.insert(scored.end(), local.begin(), local.end());
    }

    std::sort(scored.begin(), scored.end(),
        [](const DistPair& a, const DistPair& b) { return a.distance < b.distance; });

    int k = std::min(topk, (int)scored.size());
    vector<RepDBQueryResult> results;
    results.reserve(k);

    for (int i = 0; i < k; i++) {
        RepDBQueryResult r;
        r.rep_idx = scored[i].rep_idx;
        r.distance = scored[i].distance;
        if (r.rep_idx >= 0 && r.rep_idx < (int)representatives.size()) {
            r.genome_id = representative_ids[r.rep_idx];
            r.genome_name = representatives[r.rep_idx].fileName;
        } else {
            r.genome_id = -1;
            r.genome_name = "N/A";
        }
        r.cluster_id = r.rep_idx;
        r.cluster_size = (r.rep_idx >= 0 && r.rep_idx < (int)clusters.size())
                         ? (int)clusters[r.rep_idx].size() : 0;
        results.push_back(r);
    }
    return results;
}

RepDBQueryResult MinHashClusterState::assign(
    const SketchInfo& query, int threads) const
{
    auto results = query_topk(query, 1, threads);
    if (!results.empty() && results[0].distance <= threshold) {
        return results[0];
    }
    RepDBQueryResult none;
    none.rep_idx = -1;
    none.genome_id = -1;
    none.genome_name = "unassigned";
    none.distance = -1.0;
    none.cluster_id = -1;
    none.cluster_size = 0;
    return none;
}

void MinHashClusterState::print_stats(std::ostream& out) const {
    size_t total_genomes = 0;
    for (const auto& cl : clusters) total_genomes += cl.size();

    out << "========================================" << std::endl;
    out << "    MinHash RepDB Statistics Report" << std::endl;
    out << "========================================" << std::endl;
    out << std::endl;

    out << "[Basic Info]" << std::endl;
    out << "  Threshold:              " << threshold << std::endl;
    out << "  Kmer size:              " << kmer_size << std::endl;
    out << "  Sketch size:            " << sketch_size << std::endl;
    out << "  Containment mode:       " << (is_containment ? "yes" : "no") << std::endl;
    out << std::endl;

    out << "[Scale]" << std::endl;
    out << "  Total genomes:          " << total_genomes << std::endl;
    out << "  Representatives:        " << representatives.size() << std::endl;
    out << "  Clusters:               " << clusters.size() << std::endl;
    double compression = total_genomes > 0
        ? (1.0 - (double)representatives.size() / total_genomes) * 100.0 : 0.0;
    out << "  Compression ratio:      " << std::fixed << std::setprecision(2)
        << compression << "%" << std::endl;
    out << std::endl;

    out << "[Inverted Index]" << std::endl;
    out << "  Unique hashes:          " << inverted_index.size() << std::endl;
    size_t total_postings = 0;
    size_t max_posting = 0;
    for (const auto& p : inverted_index) {
        total_postings += p.second.size();
        if (p.second.size() > max_posting) max_posting = p.second.size();
    }
    double avg_posting = inverted_index.size() > 0
        ? (double)total_postings / inverted_index.size() : 0;
    out << "  Total postings:         " << total_postings << std::endl;
    out << "  Avg posting length:     " << std::fixed << std::setprecision(2)
        << avg_posting << std::endl;
    out << "  Max posting length:     " << max_posting << std::endl;
    out << std::endl;

    out << "[Cluster Size Distribution]" << std::endl;
    if (!clusters.empty()) {
        std::vector<int> sizes;
        sizes.reserve(clusters.size());
        size_t singleton = 0;
        for (const auto& cl : clusters) {
            sizes.push_back(cl.size());
            if (cl.size() <= 1) singleton++;
        }
        std::sort(sizes.begin(), sizes.end());
        double mean = (double)total_genomes / clusters.size();
        int median = sizes[sizes.size() / 2];

        out << "  Min cluster size:       " << sizes.front() << std::endl;
        out << "  Max cluster size:       " << sizes.back() << std::endl;
        out << "  Mean cluster size:      " << std::fixed << std::setprecision(2)
            << mean << std::endl;
        out << "  Median cluster size:    " << median << std::endl;
        out << "  Singletons:             " << singleton << " ("
            << std::fixed << std::setprecision(1)
            << (100.0 * singleton / clusters.size()) << "%)" << std::endl;

        size_t p90 = sizes[(size_t)(sizes.size() * 0.9)];
        size_t p95 = sizes[(size_t)(sizes.size() * 0.95)];
        size_t p99 = sizes[(size_t)(sizes.size() * 0.99)];
        out << "  P90 cluster size:       " << p90 << std::endl;
        out << "  P95 cluster size:       " << p95 << std::endl;
        out << "  P99 cluster size:       " << p99 << std::endl;
    }

    uint64_t total_seq_len = 0;
    for (const auto& sk : all_sketches)
        total_seq_len += sk.totalSeqLength;
    if (total_seq_len > 0) {
        out << std::endl;
        out << "[Genome Coverage]" << std::endl;
        uint64_t rep_seq_len = 0;
        for (const auto& r : representatives)
            rep_seq_len += r.totalSeqLength;
        out << "  Total sequence length:  " << total_seq_len << " bp" << std::endl;
        out << "  Representative seq len: " << rep_seq_len << " bp" << std::endl;
        out << "  Coverage ratio:         " << std::fixed << std::setprecision(2)
            << (100.0 * rep_seq_len / total_seq_len) << "%" << std::endl;
    }

    out << "========================================" << std::endl;
}


#endif

