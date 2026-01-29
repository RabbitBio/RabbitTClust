#include <iostream>
#include "MST.h"
#include <stdio.h>
#include <queue>
#include "UnionFind.h"
#include <omp.h>
#include <fstream>
#include <sys/stat.h>
#include <climits>
#include <algorithm>
#include "phmap.h"  // Google's Swiss Tables for faster hash maps
using namespace std;

bool cmpEdge(EdgeInfo e1, EdgeInfo e2){
  return e1.dist < e2.dist;
}

bool cmpNeighbor(NeighborNode n1, NeighborNode n2){
  return n1.distance < n2.distance;
}


inline double calr(double D, int k) {
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


struct DSU {
  std::vector<int> p, r;
  DSU(int n) : p(n), r(n, 0) {
    std::iota(p.begin(), p.end(), 0);
  }
  int find(int x) {
    return p[x] == x ? x : p[x] = find(p[x]);
  }
  int unite(int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return a;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) r[a]++;
    return a;
  }
};


vector<EdgeInfo> kruskalAlgorithm(vector<EdgeInfo>graph, int vertices){
  UnionFind uf(vertices);
  vector<EdgeInfo>spanningTree;
  if(graph.size() <= 0) return spanningTree;
  //sort(graph.begin(),graph.end(),comparator);
  spanningTree.push_back(graph[0]);
  uf.merge(graph[0].preNode,graph[0].sufNode);
  for(int i=1;i<graph.size();i++){
    if(!uf.connected(graph[i].preNode,graph[i].sufNode)){
      uf.merge(graph[i].preNode,graph[i].sufNode);
      spanningTree.push_back(graph[i]);
    }
    //if(spanningTree.size() == vertices - 1) break;
  }
  uf.clear();
  return spanningTree;
}

vector<EdgeInfo> generateForest(vector <EdgeInfo> mst, double threshhold){
  vector<EdgeInfo> forest;
  for(int i = 0; i < mst.size(); i++){
    if(mst[i].dist <= threshhold){
      forest.push_back(mst[i]);
    }
  }
  return forest;
}

vector<EdgeInfo> modifyForest(vector<EdgeInfo> forest, vector<int> noiseArr, int threads){
  vector<EdgeInfo> newForest;
  bool * removeTag = new bool[forest.size()];
#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < forest.size(); i++){
    bool remove = false;
    for(int j = 0; j < noiseArr.size(); j++){
      if(forest[i].preNode == noiseArr[j] || forest[i].sufNode == noiseArr[j]){
        remove = true;
        break;
      }
    }
    removeTag[i] = remove;
  }
  for(int i = 0; i < forest.size(); i++){
    if(!removeTag[i])
      newForest.push_back(forest[i]);
  }
  return newForest;
}


vector<vector<int>> generateClusterWithBfs(vector<EdgeInfo> forest, int vertices){
  vector<vector<int>>res;

  vector<int>*G=new vector<int>[vertices];
  for(auto it:forest){
    int u=it.preNode;
    int v=it.sufNode;
    G[u].push_back(v);
    G[v].push_back(u);
  }
  bool *visited=new bool[vertices];
  memset(visited,0,sizeof(bool)*vertices);
  for(int i=0;i<vertices;i++){
    if(visited[i]==0){
      visited[i]=1;
      queue<int>Q;
      Q.push(i);
      vector<int>tmpCluster;
      tmpCluster.push_back(i);
      while(!Q.empty()){
        int k=Q.front();
        Q.pop();
        for(auto v:G[k]){
          if(visited[v])continue;
          visited[v]=1;
          Q.push(v);
          tmpCluster.push_back(v);
        }
      }
      res.push_back(tmpCluster);
    }
  }
  return res;
}

vector<vector<int>> generateCluster(vector<EdgeInfo> forest, int vertices){
  UnionFind uf(vertices);
  vector<vector<int>> cluster;
  for(int i = 0; i < forest.size(); i++){
    uf.merge(forest[i].preNode, forest[i].sufNode);
  }


  vector<int> elements;
  vector<int> remainElements;
  vector<int> tmpCluster;
  for(int i = 0; i < vertices; i++){
    elements.push_back(i);
  }

  while(!elements.empty()){
    tmpCluster.push_back(elements[0]);
    for(int i = 1; i < elements.size(); i++){
      if(uf.connected(elements[0], elements[i])){
        tmpCluster.push_back(elements[i]);
      }
      else{
        remainElements.push_back(elements[i]);
      }

    }
    if(elements.size() == 1){//the last element which has not been clustered.
      cluster.push_back(elements);
      elements.clear();
      break;
    }
    cluster.push_back(tmpCluster);
    tmpCluster.clear();
    elements = remainElements;
    remainElements.clear();
  }

  return cluster;
}

//typedef pair<int, int> PairInt;
bool cmpPair(PairInt p1, PairInt p2){
  return p1.second < p2.second;
}

vector<int> getNoiseNode(vector<PairInt> densePairArr, int alpha){
  int clusterSize = densePairArr.size();
  vector<int> noiseArr;
  std::sort(densePairArr.begin(), densePairArr.end(), cmpPair);

  int sizeDense = densePairArr.size();
  int indexQ1 = sizeDense / 4;
  int denseQ1 = densePairArr[indexQ1].second;
  //int thresholdDense = std::min((int)round(clusterSize*alpha), 2); 
  int thresholdDense = alpha;
  thresholdDense = std::min(denseQ1-1, thresholdDense);
  thresholdDense = std::max(thresholdDense, 0);

  for(int i = 0; i < densePairArr.size(); i++){
    if(densePairArr[i].second <= thresholdDense){
      noiseArr.push_back(densePairArr[i].first);
    }
    else{
      break;
    }
  }
  return noiseArr;
}



  // init from memory inverted index or dictFile and indexFile
vector<EdgeInfo> compute_kssd_mst(vector<KssdSketchInfo>& sketches, KssdParameters info, const string folder_path, int start_index, bool no_dense, bool isContainment, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr, double threshold, KssdInvertedIndex* inverted_index){  int half_k = info.half_k;
  int drlevel = info.drlevel;
  bool use64 = half_k - drlevel > 8 ? true : false;
  int kmer_size = half_k * 2;
  phmap::flat_hash_map<uint64_t, vector<uint32_t>> hash_map_arr;
  uint32_t* sketchSizeArr = NULL;
  size_t* offset = NULL;
  uint32_t* indexArr = NULL;
  int radio = calr(threshold, kmer_size-1);
  
  // Precompute constants for distance calculation
  const double inv_kmer_size = 1.0 / kmer_size;
  const double log_2 = log(2.0); 
  
  // Use memory index if available, otherwise load from file
  const phmap::flat_hash_map<uint64_t, vector<uint32_t>>* hash_map_ptr = nullptr;
  if(inverted_index != nullptr && inverted_index->use64 == use64) {
    // Use memory index directly
    if(use64) {
      hash_map_ptr = &(inverted_index->hash_map_64);  // Use pointer to avoid copy
    } else {
      // For non-use64, need to build sketchSizeArr, offset, and indexArr from memory index
      size_t hashSize = inverted_index->hashSize;
      sketchSizeArr = new uint32_t[hashSize];
      offset = new size_t[hashSize];
      uint64_t totalHashNumber = 0;
      
      for(size_t i = 0; i < hashSize; i++){
        sketchSizeArr[i] = inverted_index->hash_map_32[i].size();
        totalHashNumber += sketchSizeArr[i];
        offset[i] = sketchSizeArr[i];
        if(i > 0) offset[i] += offset[i-1];
      }
      
      indexArr = new uint32_t[totalHashNumber];
      size_t idx = 0;
      for(size_t i = 0; i < hashSize; i++){
        for(size_t j = 0; j < inverted_index->hash_map_32[i].size(); j++){
          indexArr[idx++] = inverted_index->hash_map_32[i][j];
        }
      }
    }
    cerr << "-----using memory inverted index in compute_kssd_mst()" << endl;
  } else if(use64) {
    // Load from file (original code)
    size_t hash_number;
    string cur_index_file = folder_path + '/' + "kssd.sketch.index";
    FILE* fp_index = fopen(cur_index_file.c_str(), "rb");
    if(!fp_index){
      cerr << "ERROR: compute_kssd_mst(), cannot open index file: " << cur_index_file << endl;
      exit(1);
    }
    int read_hash_num = fread(&hash_number, sizeof(size_t), 1, fp_index);
    uint64_t * hash_arr = new uint64_t[hash_number];
    uint32_t * hash_size_arr = new uint32_t[hash_number];
    size_t read_hash_arr = fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
    size_t read_hash_size_arr = fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
    if(read_hash_num != 1 || read_hash_arr != hash_number || read_hash_size_arr != hash_number){
      cerr << "ERROR: compute_kssd_mst(), error read hash_number, hash_arr, and hash_size_arr" << endl;
      exit(1);
    }

    fclose(fp_index);

    string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
    FILE* fp_dict = fopen(cur_dict_file.c_str(), "rb");
    if(!fp_dict){
      cerr << "ERROR: compute_kssd_mst(), cannot open dict file: " << cur_dict_file << endl;
      exit(1);
    }
    uint32_t max_hash_size = 1LLU << 24;
    uint32_t* cur_point = new uint32_t[max_hash_size];
    for(size_t i = 0; i < hash_number; i++){
      uint32_t cur_hash_size = hash_size_arr[i];
      if(cur_hash_size > max_hash_size){
        max_hash_size = cur_hash_size;
        cur_point = new uint32_t[max_hash_size];
      }
      uint32_t hash_size = fread(cur_point, sizeof(uint32_t), cur_hash_size, fp_dict);
      if(hash_size != cur_hash_size){
        cerr << "ERROR: compute_kssd_mst(), the read hash number is not equal to the saved hash number information" << endl;
        exit(1);
      }
      vector<uint32_t> cur_genome_arr(cur_point, cur_point + cur_hash_size);
      uint64_t cur_hash = hash_arr[i];
      hash_map_arr.insert({cur_hash, cur_genome_arr});
    }
    delete [] cur_point;
    delete [] hash_arr;
    delete [] hash_size_arr;
    fclose(fp_dict);
  }
  else{
    // Load from file for non-use64 (original code)
    cerr << "-----not use hash64 in compute_kssd_mst() " << endl;
    size_t hashSize;
    uint64_t totalIndex;
    string cur_index_file = folder_path + '/' + "kssd.sketch.index";
    FILE * fp_index = fopen(cur_index_file.c_str(), "rb");
    if(!fp_index){
      cerr << "ERROR: compute_kssd_mst(), cannot open the index sketch file: " << cur_index_file << endl;
      exit(1);
    }
    int read_hash_size = fread(&hashSize, sizeof(size_t), 1, fp_index);
    int read_total_index = fread(&totalIndex, sizeof(uint64_t), 1, fp_index);
    //sketchSizeArr = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
    sketchSizeArr = new uint32_t[hashSize];
    size_t read_sketch_size_arr = fread(sketchSizeArr, sizeof(uint32_t), hashSize, fp_index);

    //offset = (size_t*)malloc(hashSize * sizeof(size_t));
    offset = new size_t[hashSize];
    uint64_t totalHashNumber = 0;
    for(size_t i = 0; i < hashSize; i++){
      totalHashNumber += sketchSizeArr[i];
      offset[i] = sketchSizeArr[i];
      if(i > 0) offset[i] += offset[i-1];
    }
    if(totalHashNumber != totalIndex){
      cerr << "ERROR: compute_kssd_mst(), mismatched total hash number" << endl;
      exit(1);
    }
    fclose(fp_index);

    //cerr << "the hashSize is: " << hashSize << endl;
    //cerr << "totalIndex is: " << totalIndex << endl;
    //cerr << "totalHashNumber is: " << totalHashNumber << endl;
    //cerr << "offset[n-1] is: " << offset[hashSize-1] << endl;;

    //indexArr = (uint32_t*)malloc(totalHashNumber * sizeof(uint32_t));
    indexArr = new uint32_t[totalHashNumber];
    string cur_dict_file = folder_path + '/' + "kssd.sketch.dict";
    FILE * fp_dict = fopen(cur_dict_file.c_str(), "rb");
    if(!fp_dict){
      cerr << "ERROR: compute_kssd_mst(), cannot open the dictionary sketch file: " << cur_dict_file << endl;
      exit(1);
    }
    size_t read_index_arr = fread(indexArr, sizeof(uint32_t), totalHashNumber, fp_dict);
    if(read_hash_size != 1 || read_total_index != 1 || read_sketch_size_arr != hashSize || read_index_arr != totalHashNumber){
      cerr << "ERROR: compute_kssd_mst(), error read hash_size, total_index, sketch_size_arr, index_arr" << endl;
      exit(1);
    }
    fclose(fp_dict);
  }

  //int denseSpan = 10;
  double step = 1.0 / denseSpan;

  //double step = threshold / denseSpan;
  //cerr << "the threshold is: " << threshold << endl;
  //cerr << "the step is: " << step << endl;
  int N = sketches.size();
  denseArr = new int*[denseSpan];
  int** denseLocalArr = new int*[denseSpan * threads];
  double distRadius[denseSpan];
  for(int i = 0; i < denseSpan; i++){
    distRadius[i] = step * i;
    denseArr[i] = new int[N];
    memset(denseArr[i], 0, N * sizeof(int));
    for(int j = 0; j < threads; j++){
      denseLocalArr[i * threads + j] = new int[N];
      memset(denseLocalArr[i*threads+j], 0, N * sizeof(int));
    }
  }
  //for ANI distribution calculation.
  aniArr = new uint64_t[101];
  memset(aniArr, 0, 101 * sizeof(uint64_t));
  uint64_t** threadsANI = new uint64_t*[threads];
  for(int i = 0; i < threads; i++){
    threadsANI[i] = new uint64_t[101];
    memset(threadsANI[i], 0, 101 * sizeof(uint64_t));
  }

  // start to generate the sub_mst
  std::vector<std::vector<EdgeInfo>> mstArr(threads);
  int** intersectionArr = new int*[threads];
  for(int i = 0; i < threads; i++){
    intersectionArr[i] = new int[sketches.size()];
  }
  int** seenStamp = new int*[threads];
  int*  epoch     = new int[threads];
  std::vector<int>* touched = new std::vector<int>[threads];

  for(int t = 0; t < threads; ++t){
    seenStamp[t] = new int[N];
    memset(seenStamp[t], 0, N * sizeof(int));
    epoch[t] = 1;
    touched[t].reserve(4096);
  }

  int subSize = 8;
  int id = 0;
  //int tailNum = sketches.size() % subSize;
  int tailNum = (N - start_index) % subSize;
  //int N = sketches.size();
  //uint64_t totalCompNum = (uint64_t)N * (uint64_t)(N-1)/2;
  uint64_t totalCompNum = (uint64_t)(N-start_index) * (uint64_t)(N+start_index)/2;
  uint64_t percentNum = totalCompNum / 100;
  cerr << "---the percentNum is: " << percentNum << endl;
  cerr << "---the start_index is: " << start_index << endl;
  uint64_t percentId = 0;
#pragma omp parallel for num_threads(threads) schedule (dynamic)
  for(id = start_index; id < sketches.size() - tailNum; id+=subSize){
    int thread_id = omp_get_thread_num();
    for(int i = id; i < id + subSize; i++){
      // ===== NEW: per-thread alias =====
      int* inter = intersectionArr[thread_id];
      int* stamp = seenStamp[thread_id];
      int& ep     = epoch[thread_id];
      auto& cand  = touched[thread_id];
    
      // ===== NEW: start a new epoch, no memset(inter,0,N) =====
      cand.clear();
      ++ep;
      if(__builtin_expect(ep == INT_MAX, 0)){ 
        memset(stamp, 0, N * sizeof(int));
        ep = 1;
      }
    
      // Cache sketch size for current genome
      int size0;
      if(use64){
        size0 = (int)sketches[i].hash64_arr.size();
      }else{
        size0 = (int)sketches[i].hash32_arr.size();
      }
      
      // Early skip if sketch is empty
      if(__builtin_expect(size0 == 0, 0)) continue;
    
      // ===== NEW: accumulate intersection counts AND collect candidates =====
      if(use64){
        const auto& hash_arr_i = sketches[i].hash64_arr;
        const size_t hash_size = hash_arr_i.size();
        for(size_t j = 0; j < hash_size; j++){
          uint64_t hash64 = hash_arr_i[j];
          
          // Prefetch next hash for better cache performance
          if(__builtin_expect(j + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[j + 1], 0, 1);
          }
    
          const phmap::flat_hash_map<uint64_t, vector<uint32_t>>* map_to_use = hash_map_ptr ? hash_map_ptr : &hash_map_arr;
          auto it = map_to_use->find(hash64);
          if(__builtin_expect(it == map_to_use->end(), 0)) continue;
    
          const auto& vec = it->second;
          const size_t vec_size = vec.size();
          for(size_t k = 0; k < vec_size; k++){
            int cur = (int)vec[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      }
      else{
        const auto& hash_arr_i = sketches[i].hash32_arr;
        const size_t hash_size = hash_arr_i.size();
        for(size_t j = 0; j < hash_size; j++){
          uint32_t hash = hash_arr_i[j];
          
          // Prefetch next hash
          if(__builtin_expect(j + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[j + 1], 0, 1);
          }
          
          if(__builtin_expect(sketchSizeArr[hash] <= 1, 0)) continue;
    
          size_t start = (hash > 0) ? offset[hash-1] : 0;
          size_t end   = offset[hash];
          
          // Prefetch index array region
          if(__builtin_expect(end > start, 1)){
            __builtin_prefetch(&indexArr[start], 0, 1);
          }
    
          for(size_t k = start; k < end; k++){
            int cur = (int)indexArr[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      }
    
      // ===== NEW: candidate-only loop (replace for(j=0;j<i;j++)) =====
      const int cand_size = (int)cand.size();
      for(int idx = 0; idx < cand_size; ++idx){
        int j = cand[idx];
        if(__builtin_expect(j >= i, 0)) continue; 
    
        double tmpDist;
        int common = inter[j];
    
        int size1;
        if(use64){
          size1 = (int)sketches[j].hash64_arr.size();
        }else{
          size1 = (int)sketches[j].hash32_arr.size();
        }
    
        // Fast path: skip if size ratio check fails
        if(__builtin_expect(size0 > 0 && size1 > 0, 1)){
          int min_size = size0 < size1 ? size0 : size1;
          int max_size = size0 > size1 ? size0 : size1;
          if(__builtin_expect(max_size > radio * min_size, 0)) continue;
        }else{
          continue;
        }
    
        if(!isContainment){
          int denom = size0 + size1 - common;
          double jaccard;
          if(__builtin_expect(denom == 0, 0)) jaccard = 0.0;
          else jaccard = (double)common / denom;
    
          double mashD;
          if(__builtin_expect(jaccard == 1.0, 0)) mashD = 0.0;
          else if(__builtin_expect(jaccard == 0.0, 0)) mashD = 1.0;
          else {
            // Optimized: use precomputed constants
            double ratio = (2.0 * jaccard) / (1.0 + jaccard);
            mashD = -inv_kmer_size * log(ratio);
          }
          tmpDist = mashD;
        }else{
          int denom = size0 < size1 ? size0 : size1;
          double containment;
          if(__builtin_expect(denom == 0, 0)) containment = 0.0;
          else containment = (double)common / denom;
    
          double AafD;
          if(__builtin_expect(containment == 1.0, 0)) AafD = 0.0;
          else if(__builtin_expect(containment == 0.0, 0)) AafD = 1.0;
          else AafD = -inv_kmer_size * log(containment);
          tmpDist = AafD;
        }
    
        if(!no_dense){
          // Dense semantics are cumulative over thresholds:
          // base: for each t, if(tmpDist <= distRadius[t]) ++
          // Optimized: store only the first satisfied bucket t0, then prefix-sum later.
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(__builtin_expect(t0 < denseSpan, 1)){
            denseLocalArr[t0 * threads + thread_id][i]++;
            denseLocalArr[t0 * threads + thread_id][j]++;
          }

          // Calculate ANI (optimized: use multiplication instead of division)
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI * 100.0);  // 100.0 = 1/0.01
          if(__builtin_expect(ANI >= 101, 0)) ANI = 100;
          threadsANI[thread_id][ANI]++;
        }
    
        mstArr[thread_id].push_back(EdgeInfo{i, j, tmpDist});
      }
    }    
    if(thread_id == 0){
      //uint64_t computedNum = (uint64_t)(N - id) * (uint64_t)id + (uint64_t)id * (uint64_t)id / 2;
      uint64_t computedNum = (uint64_t)(id-start_index) * (uint64_t)(id+start_index) / 2;
      if(computedNum >= percentId * percentNum){
        fprintf(stderr, "---finish MST generation %d %\n", percentId);
        percentId++;
      }
    }

    sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
    vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
    mstArr[thread_id].swap(tmpMst);
    vector<EdgeInfo>().swap(tmpMst);
  }
  cerr << "-----finish the 100 % multiThreads mst generate" << endl;

  if(!no_dense){
    for(int i = 0; i < 101; i++){
      for(int j = 0; j < threads; j++){
        aniArr[i] += threadsANI[j][i];
      }
    }

    // Merge per-thread start-bucket counts (optimized: cache-friendly order + parallel k).
    for(int t = 0; t < denseSpan; t++){
      #pragma omp parallel for schedule(static) num_threads(threads)
      for(int k = 0; k < N; k++){
        int sum = 0;
        for(int th = 0; th < threads; th++){
          sum += denseLocalArr[t*threads + th][k];
        }
        denseArr[t][k] = sum;
      }
    }
  }
  for(int i = 0; i < threads; i++){
    delete [] threadsANI[i];
  }
  delete [] threadsANI;
  
  for(int i = 0; i < denseSpan * threads; i++){
    delete [] denseLocalArr[i];
  }
  delete [] denseLocalArr;

  if(tailNum != 0){
    for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
      int* inter = intersectionArr[0];
      int* stamp = seenStamp[0];
      int& ep     = epoch[0];
      auto& cand  = touched[0];
      
      cand.clear();
      ++ep;
      if(__builtin_expect(ep == INT_MAX, 0)){
        memset(stamp, 0, N * sizeof(int));
        ep = 1;
      }
      
      // Cache sketch size
      int size0;
      if(use64){
        size0 = (int)sketches[i].hash64_arr.size();
      }else{
        size0 = (int)sketches[i].hash32_arr.size();
      }
      
      if(__builtin_expect(size0 == 0, 0)) continue;
      
      if(use64){
        const auto& hash_arr_i = sketches[i].hash64_arr;
        const size_t hash_size = hash_arr_i.size();
        for(size_t jj = 0; jj < hash_size; jj++){
          uint64_t hash64 = hash_arr_i[jj];
          
          if(__builtin_expect(jj + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[jj + 1], 0, 1);
          }
          
          const phmap::flat_hash_map<uint64_t, vector<uint32_t>>* map_to_use = hash_map_ptr ? hash_map_ptr : &hash_map_arr;
          auto it = map_to_use->find(hash64);
          if(__builtin_expect(it == map_to_use->end(), 0)) continue;
      
          const auto& vec = it->second;
          const size_t vec_size = vec.size();
          for(size_t k = 0; k < vec_size; k++){
            int cur = (int)vec[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      }
      else{
        const auto& hash_arr_i = sketches[i].hash32_arr;
        const size_t hash_size = hash_arr_i.size();
        for(size_t jj = 0; jj < hash_size; jj++){
          uint32_t hash = hash_arr_i[jj];
          
          if(__builtin_expect(jj + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[jj + 1], 0, 1);
          }
          
          if(__builtin_expect(sketchSizeArr[hash] <= 1, 0)) continue;
      
          size_t start = (hash > 0) ? offset[hash-1] : 0;
          size_t end   = offset[hash];
          
          if(__builtin_expect(end > start, 1)){
            __builtin_prefetch(&indexArr[start], 0, 1);
          }
      
          for(size_t k = start; k < end; k++){
            int cur = (int)indexArr[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      }
      
      // candidate-only
      const int cand_size = (int)cand.size();
      for(int idx = 0; idx < cand_size; ++idx){
        int j = cand[idx];
        if(__builtin_expect(j >= i, 0)) continue;
      
        double tmpDist;
        int common = inter[j];
      
        int size1;
        if(use64){
          size1 = (int)sketches[j].hash64_arr.size();
        }else{
          size1 = (int)sketches[j].hash32_arr.size();
        }
      
        // Fast path: skip if size ratio check fails
        if(__builtin_expect(size0 > 0 && size1 > 0, 1)){
          int min_size = size0 < size1 ? size0 : size1;
          int max_size = size0 > size1 ? size0 : size1;
          if(__builtin_expect(max_size > radio * min_size, 0)) continue;
        }else{
          continue;
        }
      
        if(!isContainment){
          int denom = size0 + size1 - common;
          double jaccard;
          if(__builtin_expect(denom == 0, 0)) jaccard = 0.0;
          else jaccard = (double)common / denom;
      
          double mashD;
          if(__builtin_expect(jaccard == 1.0, 0)) mashD = 0.0;
          else if(__builtin_expect(jaccard == 0.0, 0)) mashD = 1.0;
          else {
            double ratio = (2.0 * jaccard) / (1.0 + jaccard);
            mashD = -inv_kmer_size * log(ratio);
          }
          tmpDist = mashD;
        }else{
          int denom = size0 < size1 ? size0 : size1;
          double containment;
          if(__builtin_expect(denom == 0, 0)) containment = 0.0;
          else containment = (double)common / denom;
      
          double AafD;
          if(__builtin_expect(containment == 1.0, 0)) AafD = 0.0;
          else if(__builtin_expect(containment == 0.0, 0)) AafD = 1.0;
          else AafD = -inv_kmer_size * log(containment);
          tmpDist = AafD;
        }
      
        if(!no_dense){
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(__builtin_expect(t0 < denseSpan, 1)){
            denseArr[t0][i]++;
            denseArr[t0][j]++;
          }
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI * 100.0);
          if(__builtin_expect(ANI >= 101, 0)) ANI = 100;
          aniArr[ANI]++;
        }
      
        mstArr[0].push_back(EdgeInfo{i, j, tmpDist});
      }
      
    }
    if(mstArr[0].size() != 0){
      sort(mstArr[0].begin(), mstArr[0].end(), cmpEdge);
      vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[0], sketches.size());
      mstArr[0].swap(tmpMst);
    }

  }

  // Convert start-bucket counts to cumulative density counts to match base semantics.
  if(!no_dense){
    #pragma omp parallel for schedule(static) num_threads(threads)
    for(int k = 0; k < N; k++){
      int acc = 0;
      for(int t = 0; t < denseSpan; t++){
        acc += denseArr[t][k];
        denseArr[t][k] = acc;
      }
    }
  }

  vector<EdgeInfo> finalGraph;
  for(int i = 0; i < threads; i++){
    finalGraph.insert(finalGraph.end(), mstArr[i].begin(), mstArr[i].end());
    vector<EdgeInfo>().swap(mstArr[i]);
  }

  sort(finalGraph.begin(), finalGraph.end(), cmpEdge);

  vector<EdgeInfo> mst = kruskalAlgorithm(finalGraph, sketches.size());
  vector<EdgeInfo>().swap(finalGraph);
  for(int t = 0; t < threads; ++t){
    delete [] seenStamp[t];
  }
  delete [] seenStamp;
  delete [] epoch;
  delete [] touched;
  for(int t = 0; t < threads; ++t){
    delete [] intersectionArr[t];
  }
  delete [] intersectionArr;
  
  if(!use64){
    delete [] sketchSizeArr;
    delete [] offset;
    delete [] indexArr;
  }

  return mst;
}

vector<EdgeInfo> modifyMST(vector<SketchInfo>& sketches, int start_index, int sketch_func_id, int threads, bool no_dense, int** &denseArr, int denseSpan, uint64_t* &aniArr){
  //int denseSpan = 10;
  double step = 1.0 / denseSpan;

  //double step = threshold / denseSpan;
  //cerr << "the threshold is: " << threshold << endl;
  //cerr << "the step is: " << step << endl;
  int N = sketches.size();
  denseArr = new int*[denseSpan];
  int** denseLocalArr = new int*[denseSpan * threads];
  double distRadius[denseSpan];
  for(int i = 0; i < denseSpan; i++){
    //distRadius[i] = 0.05 - 0.005 * i;
    //distRadius[i] = threshold - step * i;
    distRadius[i] = step * i;
    denseArr[i] = new int[N];
    memset(denseArr[i], 0, N * sizeof(int));
    for(int j = 0; j < threads; j++){
      denseLocalArr[i * threads + j] = new int[N];
      memset(denseLocalArr[i*threads+j], 0, N * sizeof(int));
    }
  }
  //for ANI distribution calculation.
  aniArr = new uint64_t[101];
  memset(aniArr, 0, 101 * sizeof(uint64_t));
  uint64_t** threadsANI = new uint64_t*[threads];
  //uint64_t threadsANI[threads][101];
  for(int i = 0; i < threads; i++){
    threadsANI[i] = new uint64_t[101];
    memset(threadsANI[i], 0, 101 * sizeof(uint64_t));
  }

  vector<EdgeInfo> mstArr[threads];
  int subSize = 8;
  int id = 0;
  int tailNum = sketches.size() % subSize;
  //int N = sketches.size();
  uint64_t totalCompNum = (uint64_t)N * (uint64_t)(N-1)/2;
  uint64_t percentNum = totalCompNum / 100;
  cerr << "---the percentNum is: " << percentNum << endl;
  cerr << "---the start_index is: " << start_index << endl;
  uint64_t percentId = 0;
#pragma omp parallel for num_threads(threads) schedule (dynamic)
  for(id = 0; id < sketches.size() - tailNum; id+=subSize){
    int thread_id = omp_get_thread_num();
    for(int i = id; i < id+subSize; i++){
      for(int j = max(i+1, start_index); j < sketches.size(); j++){
        double tmpDist;
        if(sketch_func_id == 0)
        {
          if(sketches[i].isContainment)
          {
            //tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
            tmpDist = sketches[i].minHash->containDistance(sketches[j].minHash);
          }
          else
          {
            tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
          }
        }
        else if(sketch_func_id == 1){
          tmpDist = sketches[i].KSSD->distance(sketches[j].KSSD);
        }
        else if(sketch_func_id == 2){
          tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
        }
        else if(sketch_func_id == 3){
          tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
        }
        else if(sketch_func_id == 4){
          tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
        }
        else	
          break;

        if(!no_dense){
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(t0 < denseSpan){
            denseLocalArr[t0 * threads + thread_id][i]++;
            denseLocalArr[t0 * threads + thread_id][j]++;
          }
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI / 0.01);
          assert(ANI < 101);
          threadsANI[thread_id][ANI]++;
        }

        EdgeInfo tmpE{i, j, tmpDist};
        mstArr[thread_id].push_back(tmpE);
      }
    }
    if(thread_id == 0){
      uint64_t computedNum = (uint64_t)(N - id) * (uint64_t)id + (uint64_t)id * (uint64_t)id / 2;
      if(computedNum >= percentId * percentNum){
        fprintf(stderr, "---finish MST generation %d %\n", percentId);
        percentId++;
      }
    }

    sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
    vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
    mstArr[thread_id].swap(tmpMst);
    vector<EdgeInfo>().swap(tmpMst);
  }
  cerr << "-----finish the 100 % multiThreads mst generate" << endl;

  if(!no_dense){
    for(int i = 0; i < 101; i++){
      for(int j = 0; j < threads; j++){
        aniArr[i] += threadsANI[j][i];
      }
    }
    // Merge per-thread start-bucket counts (optimized: cache-friendly order + parallel k).
    for(int t = 0; t < denseSpan; t++){
      #pragma omp parallel for schedule(static) num_threads(threads)
      for(int k = 0; k < N; k++){
        int sum = 0;
        for(int th = 0; th < threads; th++){
          sum += denseLocalArr[t*threads + th][k];
        }
        denseArr[t][k] = sum;
      }
    }
  }
  for(int i = 0; i < threads; i++){
    delete [] threadsANI[i];
  }
  delete [] threadsANI;

  for(int i = 0; i < denseSpan * threads; i++){
    delete [] denseLocalArr[i];
  }
  delete [] denseLocalArr;

  if(tailNum != 0){
    for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
      for(int j = i+1; j < sketches.size(); j++){
        //double tmpDist = 1.0 - minHashes[i].minHash->jaccard(minHashes[j].minHash);
        double tmpDist;
        if(sketch_func_id == 0){
          if(sketches[i].isContainment)
            //tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
            tmpDist = sketches[i].minHash->containDistance(sketches[j].minHash);
          else
          {
            tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
          }
        }
        else if(sketch_func_id == 1)
          tmpDist = sketches[i].KSSD->distance(sketches[j].KSSD);
        else if(sketch_func_id == 2) 
          tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
        else if(sketch_func_id == 3)
          tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
        else if(sketch_func_id == 4)
          tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
        else	
          break;

        if(!no_dense){
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(t0 < denseSpan){
            denseArr[t0][i]++;
            denseArr[t0][j]++;
          }
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI/0.01);
          assert(ANI < 101);
          aniArr[ANI]++;
        }

        EdgeInfo tmpE{i, j, tmpDist};
        mstArr[0].push_back(tmpE);
      }
    }
    if(mstArr[0].size() != 0){
      sort(mstArr[0].begin(), mstArr[0].end(), cmpEdge);
      vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[0], sketches.size());
      mstArr[0].swap(tmpMst);
    }

  }

  // Convert start-bucket counts to cumulative density counts to match base semantics.
  if(!no_dense){
    #pragma omp parallel for schedule(static) num_threads(threads)
    for(int k = 0; k < N; k++){
      int acc = 0;
      for(int t = 0; t < denseSpan; t++){
        acc += denseArr[t][k];
        denseArr[t][k] = acc;
      }
    }
  }

  vector<EdgeInfo> finalGraph;
  for(int i = 0; i < threads; i++){
    finalGraph.insert(finalGraph.end(), mstArr[i].begin(), mstArr[i].end());
    vector<EdgeInfo>().swap(mstArr[i]);
  }

  sort(finalGraph.begin(), finalGraph.end(), cmpEdge);

  vector<EdgeInfo> mst = kruskalAlgorithm(finalGraph, sketches.size());
  vector<EdgeInfo>().swap(finalGraph);

  return mst;
}

// string build_kssd_newick_tree(vector<pair<int, double>>* G, bool* visited, const vector<KssdSketchInfo>& sketches, bool sketch_by_file, int node){
// 	visited[node] = true;
// 	//if(G[node].size() == 0)	return to_string(node);
// 	string name;
// 	if(sketch_by_file) name = sketches[node].fileName;
// 	else name = sketches[node].seqInfo.name;
// 	if(G[node].size() == 0){
// 		return name;
// 	}
// 	string children("");
// 	for(auto x : G[node]){
// 		int cur_node = x.first;
// 		double dist = x.second;
// 		if(!visited[cur_node]){
// 			string child = build_kssd_newick_tree(G, visited, sketches, sketch_by_file, cur_node);
// 			children += child + ':' + to_string(dist) + ',';
// 		}
// 	}
// 	//if(children.length() == 0) return to_string(node);
// 	if(children.length() == 0) return name;
// 	//children = '(' + children.substr(0, children.length()-1) + ')' + to_string(node);
// 	children = '(' + children.substr(0, children.length()-1) + ')' + name;
// 	return children;
// }
static string build_newick_tree_recursive(
    int node,
    const vector<vector<pair<int,double>>>& children,
    const vector<SketchInfo>& sketches,
    bool sketch_by_file)
{
  if (children[node].empty()) {
    string name = sketch_by_file ? sketches[node].fileName
      : sketches[node].seqInfo.name;
    return name;
  }

  string s = "(";
  for (size_t i = 0; i < children[node].size(); ++i) {
    int child = children[node][i].first;
    double bl  = children[node][i].second;
    if (i > 0) s += ",";
    s += build_newick_tree_recursive(child, children, sketches, sketch_by_file);
    s += ":" + to_string(bl);
  }
  s += ")";  
  return s;
}


// string get_kssd_newick_tree(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file){
// 	int vertice_number = sketches.size();
// 	vector<pair<int, double>>* G = new vector<pair<int, double>> [vertice_number];
// 	for(auto f : mst){
// 		int u = f.preNode;
// 		int v = f.sufNode;
// 		double dist = f.dist;
// 		G[u].push_back(make_pair(v, dist));
// 		G[v].push_back(make_pair(u, dist));
// 	}
// 	bool* visited = new bool[vertice_number];
// 	memset(visited, 0, vertice_number*sizeof(bool));

// 	string newick_tree = build_kssd_newick_tree(G, visited, sketches, sketch_by_file, 0);
// 	newick_tree += ';';

// 	return newick_tree;
// }

string get_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file){
  int N = sketches.size();
  if (N == 0) return ";";
  if (N == 1) {
    string name = sketch_by_file ? sketches[0].fileName
      : sketches[0].seqInfo.name;
    return name + ";";
  }

  vector<EdgeInfo> edges = mst;
  sort(edges.begin(), edges.end(),
      [](const EdgeInfo& a, const EdgeInfo& b){
      return a.dist < b.dist;
      });

  int maxNodes = 2 * N - 1;
  vector<vector<pair<int,double>>> children(maxNodes);
  vector<double> height(maxNodes, 0.0);  
  vector<int> repNode(maxNodes, -1);     
  for (int i = 0; i < N; ++i) {
    repNode[i] = i;  
  }

  DSU dsu(N);
  int nextNode = N;  

  for (const auto& e : edges) {
    int u = e.preNode;
    int v = e.sufNode;
    double w = e.dist;

    int ru = dsu.find(u);
    int rv = dsu.find(v);
    if (ru == rv) continue; 

    int nodeU = repNode[ru];
    int nodeV = repNode[rv];

    double hU = height[nodeU];
    double hV = height[nodeV];
    double h  = w;  

    double blU = max(0.0, h - hU);
    double blV = max(0.0, h - hV);

    int newNode = nextNode++;
    children[newNode].push_back({nodeU, blU});
    children[newNode].push_back({nodeV, blV});
    height[newNode] = h;

    int rNew = dsu.unite(ru, rv);
    repNode[rNew] = newNode;
  }

  int rootRep = dsu.find(0);
  int root    = repNode[rootRep];

  string newick_tree = build_newick_tree_recursive(root, children, sketches, sketch_by_file);
  newick_tree += ";";
  return newick_tree;
}

static std::string build_kssd_newick_tree(
    int node,
    const std::vector<std::vector<std::pair<int,double>>>& children,
    const std::vector<KssdSketchInfo>& sketches,
    bool sketch_by_file)
{
  if (children[node].empty()) {
    std::string name = sketch_by_file ? sketches[node].fileName
      : sketches[node].seqInfo.name;
    return name;
  }

  std::string s = "(";
  for (size_t i = 0; i < children[node].size(); ++i) {
    int child = children[node][i].first;
    double bl  = children[node][i].second;
    if (i > 0) s += ",";
    s += build_kssd_newick_tree(child, children, sketches, sketch_by_file);
    s += ":" + std::to_string(bl);
  }
  s += ")"; 
  return s;
}

std::string get_kssd_newick_tree(const std::vector<KssdSketchInfo>& sketches,
    const std::vector<EdgeInfo>& mst,
    bool sketch_by_file)
{
  int N = sketches.size();
  if (N == 0) return ";";
  if (N == 1) {
    std::string name = sketch_by_file ? sketches[0].fileName
      : sketches[0].seqInfo.name;
    return name + ";";
  }

  std::vector<EdgeInfo> edges = mst;
  std::sort(edges.begin(), edges.end(),
      [](const EdgeInfo& a, const EdgeInfo& b){
      return a.dist < b.dist;
      });


  int maxNodes = 2 * N - 1;
  std::vector<std::vector<std::pair<int,double>>> children(maxNodes);
  std::vector<double> height(maxNodes, 0.0);     
  std::vector<int> repNode(maxNodes, -1);       
  for (int i = 0; i < N; ++i) {
    repNode[i] = i;   
  }

  DSU dsu(N);
  int nextNode = N;      

  for (const auto& e : edges) {
    int u = e.preNode;
    int v = e.sufNode;
    double w = e.dist;

    int ru = dsu.find(u);
    int rv = dsu.find(v);
    if (ru == rv) continue; 

    int nodeU = repNode[ru];
    int nodeV = repNode[rv];

    double hU = height[nodeU];
    double hV = height[nodeV];
    double h  = w;  

    double blU = std::max(0.0, h - hU);
    double blV = std::max(0.0, h - hV);

    int newNode = nextNode++;
    children[newNode].push_back({nodeU, blU});
    children[newNode].push_back({nodeV, blV});
    height[newNode] = h;

    int rNew = dsu.unite(ru, rv);
    repNode[rNew] = newNode;
  }

  int rootRep = dsu.find(0);
  int root    = repNode[rootRep];

  std::string newick_tree = build_kssd_newick_tree(root, children, sketches, sketch_by_file);
  newick_tree += ";";
  return newick_tree;
}

// 从 MST 构造 single-linkage linkage matrix
vector<LinkageRow> get_linkage_from_mst(int N, const vector<EdgeInfo>& mst){
  vector<LinkageRow> Z;
  if (N <= 1) return Z;

  vector<EdgeInfo> edges = mst;
  sort(edges.begin(), edges.end(),
      [](const EdgeInfo& a, const EdgeInfo& b){
      return a.dist < b.dist;
      });

  DSU dsu(N);
  int next_cluster_id = N;
  vector<int> cluster_id(N);
  vector<int> cluster_size(2 * N - 1);

  for (int i = 0; i < N; ++i) {
    cluster_id[i]   = i;
    cluster_size[i] = 1;
  }

  Z.reserve(N - 1);

  for (const auto& e : edges) {
    int ru = dsu.find(e.preNode);
    int rv = dsu.find(e.sufNode);
    if (ru == rv) continue;

    int id_u = cluster_id[ru];
    int id_v = cluster_id[rv];

    int new_id   = next_cluster_id++;
    int new_size = cluster_size[id_u] + cluster_size[id_v];

    LinkageRow row;
    row.c1   = id_u;
    row.c2   = id_v;
    row.dist = e.dist;
    row.size = new_size;
    Z.push_back(row);

    int rNew = dsu.unite(ru, rv);
    cluster_id[rNew] = new_id;
    cluster_size[new_id] = new_size;
  }

  return Z;
}

// MinHash version of MST computation using inverted index (similar to compute_kssd_mst)
vector<EdgeInfo> compute_minhash_mst(vector<SketchInfo>& sketches, int start_index, bool no_dense, bool isContainment, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr, double threshold, int kmerSize, MinHashInvertedIndex* inverted_index){
  int kmer_size = kmerSize;
  int radio = calr(threshold, kmer_size-1);
  
  // Precompute constants for distance calculation
  const double inv_kmer_size = 1.0 / kmer_size;
  
  // Use memory index if available
  const phmap::flat_hash_map<uint64_t, vector<uint32_t>>* hash_map_ptr = nullptr;
  if(inverted_index != nullptr) {
    hash_map_ptr = &(inverted_index->hash_map);  // Use pointer to avoid copy
    cerr << "-----using memory inverted index in compute_minhash_mst()" << endl;
  } else {
    cerr << "WARNING: compute_minhash_mst() called without inverted index. Consider using modifyMST() for traditional approach." << endl;
    // Without inverted index, we can still work but need to build index on-the-fly (slower)
    // For better performance, use modifyMST() instead
  }

  // CRITICAL OPTIMIZATION: Pre-cache all hash arrays to avoid repeated storeMinHashes() calls
  // This is similar to how KSSD directly accesses hash32_arr/hash64_arr
  cerr << "-----pre-caching MinHash sketches..." << endl;
  vector<vector<uint64_t>> cached_hashes(sketches.size());
  vector<int> cached_sizes(sketches.size());
  
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 100)
  for(size_t i = 0; i < sketches.size(); i++){
    if(sketches[i].minHash != nullptr) {
      cached_hashes[i] = sketches[i].minHash->storeMinHashes();
      cached_sizes[i] = (int)cached_hashes[i].size();
    } else {
      cached_sizes[i] = 0;
    }
  }
  cerr << "-----MinHash sketches cached" << endl;

  //int denseSpan = 10;
  double step = 1.0 / denseSpan;

  int N = sketches.size();
  denseArr = new int*[denseSpan];
  int** denseLocalArr = new int*[denseSpan * threads];
  double distRadius[denseSpan];
  for(int i = 0; i < denseSpan; i++){
    distRadius[i] = step * i;
    denseArr[i] = new int[N];
    memset(denseArr[i], 0, N * sizeof(int));
    for(int j = 0; j < threads; j++){
      denseLocalArr[i * threads + j] = new int[N];
      memset(denseLocalArr[i*threads+j], 0, N * sizeof(int));
    }
  }
  //for ANI distribution calculation.
  aniArr = new uint64_t[101];
  memset(aniArr, 0, 101 * sizeof(uint64_t));
  uint64_t** threadsANI = new uint64_t*[threads];
  for(int i = 0; i < threads; i++){
    threadsANI[i] = new uint64_t[101];
    memset(threadsANI[i], 0, 101 * sizeof(uint64_t));
  }

  // start to generate the sub_mst
  std::vector<std::vector<EdgeInfo>> mstArr(threads);
  int** intersectionArr = new int*[threads];
  for(int i = 0; i < threads; i++){
    intersectionArr[i] = new int[sketches.size()];
  }
  int** seenStamp = new int*[threads];
  int*  epoch     = new int[threads];
  std::vector<int>* touched = new std::vector<int>[threads];

  for(int t = 0; t < threads; ++t){
    seenStamp[t] = new int[N];
    memset(seenStamp[t], 0, N * sizeof(int));
    epoch[t] = 1;
    touched[t].reserve(4096);
  }

  int subSize = 8;
  int id = 0;
  int tailNum = (N - start_index) % subSize;
  uint64_t totalCompNum = (uint64_t)(N-start_index) * (uint64_t)(N+start_index)/2;
  uint64_t percentNum = totalCompNum / 100;
  cerr << "---the percentNum is: " << percentNum << endl;
  cerr << "---the start_index is: " << start_index << endl;
  uint64_t percentId = 0;
  
#pragma omp parallel for num_threads(threads) schedule (dynamic)
  for(id = start_index; id < sketches.size() - tailNum; id+=subSize){
    int thread_id = omp_get_thread_num();
    for(int i = id; i < id + subSize; i++){
      // ===== per-thread alias =====
      int* inter = intersectionArr[thread_id];
      int* stamp = seenStamp[thread_id];
      int& ep     = epoch[thread_id];
      auto& cand  = touched[thread_id];
    
      // ===== start a new epoch =====
      cand.clear();
      ++ep;
      if(__builtin_expect(ep == INT_MAX, 0)){ 
        memset(stamp, 0, N * sizeof(int));
        ep = 1;
      }
    
      // Get MinHash sketch for current genome (use cached version)
      int size0 = cached_sizes[i];
      
      // Early skip if sketch is empty
      if(__builtin_expect(size0 == 0, 0)) continue;
      
      const auto& hash_arr_i = cached_hashes[i];
    
      // ===== accumulate intersection counts AND collect candidates using inverted index =====
      if(hash_map_ptr != nullptr) {
        // Use memory inverted index
        const size_t hash_size = hash_arr_i.size();
        for(size_t j = 0; j < hash_size; j++){
          uint64_t hash64 = hash_arr_i[j];
          
          // Prefetch next hash for better cache performance
          if(__builtin_expect(j + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[j + 1], 0, 1);
          }
    
          auto it = hash_map_ptr->find(hash64);
          if(__builtin_expect(it == hash_map_ptr->end(), 0)) continue;
    
          const auto& vec = it->second;
          const size_t vec_size = vec.size();
          for(size_t k = 0; k < vec_size; k++){
            int cur = (int)vec[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      } else {
        // Fallback: without pre-built index, use cached hashes to compute intersections
        // This is slower but works. For better performance, build inverted index during sketch generation.
        const size_t hash_size = hash_arr_i.size();
        for(size_t j = 0; j < sketches.size(); j++){
          if(j == i || cached_sizes[j] == 0) continue;
          const auto& hash_arr_j = cached_hashes[j];
          
          // Count intersection using cached arrays (much faster than calling storeMinHashes)
          int common = 0;
          // Use hash set for faster lookup
          phmap::flat_hash_set<uint64_t> hash_set_j(hash_arr_j.begin(), hash_arr_j.end());
          for(uint64_t hash_i : hash_arr_i) {
            if(hash_set_j.count(hash_i) > 0) {
              common++;
            }
          }
          
          if(common > 0) {
            if(__builtin_expect(stamp[j] != ep, 1)){
              stamp[j] = ep;
              inter[j] = common;
              cand.push_back(j);
            }else{
              inter[j] += common;
            }
          }
        }
      }
    
      // ===== candidate-only loop (replace for(j=0;j<i;j++)) =====
      const int cand_size = (int)cand.size();
      for(int idx = 0; idx < cand_size; ++idx){
        int j = cand[idx];
        if(__builtin_expect(j >= i, 0)) continue; 
    
        double tmpDist;
        int common = inter[j];
    
        // Use cached hash array (avoid expensive storeMinHashes() call)
        int size1 = cached_sizes[j];
        if(__builtin_expect(size1 == 0, 0)) continue;
        const auto& hash_arr_j = cached_hashes[j];
    
        // Fast path: skip if size ratio check fails
        if(__builtin_expect(size0 > 0 && size1 > 0, 1)){
          int min_size = size0 < size1 ? size0 : size1;
          int max_size = size0 > size1 ? size0 : size1;
          if(__builtin_expect(max_size > radio * min_size, 0)) continue;
        }else{
          continue;
        }
    
        if(!isContainment){
          int denom = size0 + size1 - common;
          double jaccard;
          if(__builtin_expect(denom == 0, 0)) jaccard = 0.0;
          else jaccard = (double)common / denom;
    
          double mashD;
          if(__builtin_expect(jaccard == 1.0, 0)) mashD = 0.0;
          else if(__builtin_expect(jaccard == 0.0, 0)) mashD = 1.0;
          else {
            // Optimized: use precomputed constants
            double ratio = (2.0 * jaccard) / (1.0 + jaccard);
            mashD = -inv_kmer_size * log(ratio);
          }
          tmpDist = mashD;
        }else{
          int denom = size0 < size1 ? size0 : size1;
          double containment;
          if(__builtin_expect(denom == 0, 0)) containment = 0.0;
          else containment = (double)common / denom;
    
          double AafD;
          if(__builtin_expect(containment == 1.0, 0)) AafD = 0.0;
          else if(__builtin_expect(containment == 0.0, 0)) AafD = 1.0;
          else AafD = -inv_kmer_size * log(containment);
          tmpDist = AafD;
        }
    
        if(!no_dense){
          // Dense semantics are cumulative over thresholds
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(__builtin_expect(t0 < denseSpan, 1)){
            denseLocalArr[t0 * threads + thread_id][i]++;
            denseLocalArr[t0 * threads + thread_id][j]++;
          }

          // Calculate ANI
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI * 100.0);
          if(__builtin_expect(ANI >= 101, 0)) ANI = 100;
          threadsANI[thread_id][ANI]++;
        }
    
        mstArr[thread_id].push_back(EdgeInfo{i, j, tmpDist});
      }
    }    
    if(thread_id == 0){
      uint64_t computedNum = (uint64_t)(id-start_index) * (uint64_t)(id+start_index) / 2;
      if(computedNum >= percentId * percentNum){
        fprintf(stderr, "---finish MST generation %d %\n", percentId);
        percentId++;
      }
    }

    sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
    vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
    mstArr[thread_id].swap(tmpMst);
    vector<EdgeInfo>().swap(tmpMst);
  }
  cerr << "-----finish the 100 % multiThreads mst generate" << endl;

  if(!no_dense){
    for(int i = 0; i < 101; i++){
      for(int j = 0; j < threads; j++){
        aniArr[i] += threadsANI[j][i];
      }
    }

    // Merge per-thread start-bucket counts
    for(int t = 0; t < denseSpan; t++){
      #pragma omp parallel for schedule(static) num_threads(threads)
      for(int k = 0; k < N; k++){
        int sum = 0;
        for(int th = 0; th < threads; th++){
          sum += denseLocalArr[t*threads + th][k];
        }
        denseArr[t][k] = sum;
      }
    }
  }
  for(int i = 0; i < threads; i++){
    delete [] threadsANI[i];
  }
  delete [] threadsANI;
  
  for(int i = 0; i < denseSpan * threads; i++){
    delete [] denseLocalArr[i];
  }
  delete [] denseLocalArr;

  if(tailNum != 0){
    for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
      int* inter = intersectionArr[0];
      int* stamp = seenStamp[0];
      int& ep     = epoch[0];
      auto& cand  = touched[0];
      
      cand.clear();
      ++ep;
      if(__builtin_expect(ep == INT_MAX, 0)){
        memset(stamp, 0, N * sizeof(int));
        ep = 1;
      }
      
      // Get MinHash sketch (use cached version)
      int size0 = cached_sizes[i];
      
      if(__builtin_expect(size0 == 0, 0)) continue;
      
      const auto& hash_arr_i = cached_hashes[i];
      
      if(hash_map_ptr != nullptr) {
        const size_t hash_size = hash_arr_i.size();
        for(size_t jj = 0; jj < hash_size; jj++){
          uint64_t hash64 = hash_arr_i[jj];
          
          if(__builtin_expect(jj + 1 < hash_size, 1)){
            __builtin_prefetch(&hash_arr_i[jj + 1], 0, 1);
          }
          
          auto it = hash_map_ptr->find(hash64);
          if(__builtin_expect(it == hash_map_ptr->end(), 0)) continue;
      
          const auto& vec = it->second;
          const size_t vec_size = vec.size();
          for(size_t k = 0; k < vec_size; k++){
            int cur = (int)vec[k];
            if(__builtin_expect(stamp[cur] != ep, 1)){
              stamp[cur] = ep;
              inter[cur] = 1;
              cand.push_back(cur);
            }else{
              inter[cur] += 1;
            }
          }
        }
      }
      
      // candidate-only
      const int cand_size = (int)cand.size();
      for(int idx = 0; idx < cand_size; ++idx){
        int j = cand[idx];
        if(__builtin_expect(j >= i, 0)) continue;
      
        double tmpDist;
        int common = inter[j];
      
        // Use cached hash array (avoid expensive storeMinHashes() call)
        int size1 = cached_sizes[j];
        if(__builtin_expect(size1 == 0, 0)) continue;
        const auto& hash_arr_j = cached_hashes[j];
      
        // Fast path: skip if size ratio check fails
        if(__builtin_expect(size0 > 0 && size1 > 0, 1)){
          int min_size = size0 < size1 ? size0 : size1;
          int max_size = size0 > size1 ? size0 : size1;
          if(__builtin_expect(max_size > radio * min_size, 0)) continue;
        }else{
          continue;
        }
      
        if(!isContainment){
          int denom = size0 + size1 - common;
          double jaccard;
          if(__builtin_expect(denom == 0, 0)) jaccard = 0.0;
          else jaccard = (double)common / denom;
      
          double mashD;
          if(__builtin_expect(jaccard == 1.0, 0)) mashD = 0.0;
          else if(__builtin_expect(jaccard == 0.0, 0)) mashD = 1.0;
          else {
            double ratio = (2.0 * jaccard) / (1.0 + jaccard);
            mashD = -inv_kmer_size * log(ratio);
          }
          tmpDist = mashD;
        }else{
          int denom = size0 < size1 ? size0 : size1;
          double containment;
          if(__builtin_expect(denom == 0, 0)) containment = 0.0;
          else containment = (double)common / denom;
      
          double AafD;
          if(__builtin_expect(containment == 1.0, 0)) AafD = 0.0;
          else if(__builtin_expect(containment == 0.0, 0)) AafD = 1.0;
          else AafD = -inv_kmer_size * log(containment);
          tmpDist = AafD;
        }
      
        if(!no_dense){
          int t0 = (int)(std::lower_bound(distRadius, distRadius + denseSpan, tmpDist) - distRadius);
          if(__builtin_expect(t0 < denseSpan, 1)){
            denseArr[t0][i]++;
            denseArr[t0][j]++;
          }
          double tmpANI = 1.0 - tmpDist;
          int ANI = (int)(tmpANI * 100.0);
          if(__builtin_expect(ANI >= 101, 0)) ANI = 100;
          aniArr[ANI]++;
        }
      
        mstArr[0].push_back(EdgeInfo{i, j, tmpDist});
      }
      
    }
    if(mstArr[0].size() != 0){
      sort(mstArr[0].begin(), mstArr[0].end(), cmpEdge);
      vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[0], sketches.size());
      mstArr[0].swap(tmpMst);
    }

  }

  // Convert start-bucket counts to cumulative density counts
  if(!no_dense){
    #pragma omp parallel for schedule(static) num_threads(threads)
    for(int k = 0; k < N; k++){
      int acc = 0;
      for(int t = 0; t < denseSpan; t++){
        acc += denseArr[t][k];
        denseArr[t][k] = acc;
      }
    }
  }

  vector<EdgeInfo> finalGraph;
  for(int i = 0; i < threads; i++){
    finalGraph.insert(finalGraph.end(), mstArr[i].begin(), mstArr[i].end());
    vector<EdgeInfo>().swap(mstArr[i]);
  }

  sort(finalGraph.begin(), finalGraph.end(), cmpEdge);

  vector<EdgeInfo> mst = kruskalAlgorithm(finalGraph, sketches.size());
  vector<EdgeInfo>().swap(finalGraph);
  for(int t = 0; t < threads; ++t){
    delete [] seenStamp[t];
  }
  delete [] seenStamp;
  delete [] epoch;
  delete [] touched;
  for(int t = 0; t < threads; ++t){
    delete [] intersectionArr[t];
  }
  delete [] intersectionArr;

  return mst;
}

