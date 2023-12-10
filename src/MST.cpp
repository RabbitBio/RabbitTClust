#include <iostream>
#include "MST.h"
#include <stdio.h>
#include <queue>
#include "UnionFind.h"
#include <omp.h>
#include <fstream>

using namespace std;

bool cmpEdge(EdgeInfo e1, EdgeInfo e2){
	return e1.dist < e2.dist;
}

bool cmpNeighbor(NeighborNode n1, NeighborNode n2){
	return n1.distance < n2.distance;
}


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

vector<EdgeInfo> compute_kssd_mst(vector<KssdSketchInfo>& sketches, KssdParameters info, const string folder_path, const string dictFile, const string indexFile,  int start_index, bool isContainment, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr){
	// init from the dictFile and indexFile
	int half_k = info.half_k;
	int drlevel = info.drlevel;
	bool use64 = half_k - drlevel > 8 ? true : false;
	int kmer_size = half_k * 2;
	robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
	uint32_t* sketchSizeArr = NULL;
	size_t* offset = NULL;
	uint32_t* indexArr = NULL;
	//uint64_t* dict;
	if(use64)
	{
		cerr << "-----use hash64 in compute_kssd_mst() " << endl;
		size_t hash_number;
		FILE* fp_index = fopen((folder_path+'/'+indexFile).c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: compute_kssd_mst(), cannot open index file: " << indexFile << endl;
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

		FILE* fp_dict = fopen((folder_path+'/'+dictFile).c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: compute_kssd_mst(), cannot open dict file: " << dictFile << endl;
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
		cerr << "-----not use hash64 in compute_kssd_mst() " << endl;
		size_t hashSize;
		uint64_t totalIndex;
		FILE * fp_index = fopen(indexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: compute_kssd_mst(), cannot open the index sketch file: " << indexFile << endl;
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
		FILE * fp_dict = fopen(dictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: compute_kssd_mst(), cannot open the dictionary sketch file: " << dictFile << endl;
			exit(1);
		}
		size_t read_index_arr = fread(indexArr, sizeof(uint32_t), totalHashNumber, fp_dict);
		if(read_hash_size != 1 || read_total_index != 1 || read_sketch_size_arr != hashSize || read_index_arr != totalHashNumber){
			cerr << "ERROR: compute_kssd_mst(), error read hash_size, total_index, sketch_size_arr, index_arr" << endl;
			exit(1);
		}
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
	vector<EdgeInfo> mstArr[threads];
	int** intersectionArr = new int*[threads];
	for(int i = 0; i < threads; i++){
		intersectionArr[i] = new int[sketches.size()];
	}
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
			memset(intersectionArr[thread_id], 0, sketches.size() * sizeof(int));
			if(use64){
				for(size_t j = 0; j < sketches[i].hash64_arr.size(); j++){
					uint64_t hash64 = sketches[i].hash64_arr[j];
					//if(!(dict[hash64/64] & (0x8000000000000000LLU >> (hash64 % 64))))	continue;
					if(hash_map_arr.count(hash64) == 0) continue;
					//for(auto x : hash_map_arr[hash64])
					for(size_t k = 0; k < hash_map_arr[hash64].size(); k++){
						size_t cur_index = hash_map_arr[hash64][k];
						intersectionArr[thread_id][cur_index]++;
						//cerr << hash64 << '\t' << cur_index << endl;
					}
				}
			}
			else{
				for(size_t j = 0; j < sketches[i].hash32_arr.size(); j++){
					uint32_t hash = sketches[i].hash32_arr[j];
					if(sketchSizeArr[hash] == 0) continue;
					size_t start = hash > 0 ? offset[hash-1] : 0;
					size_t end = offset[hash];
					for(size_t k = start; k < end; k++){
						size_t curIndex = indexArr[k];
						intersectionArr[thread_id][curIndex]++;
					}
				}
			}

			for(int j = max(i+1, start_index); j < sketches.size(); j++){
				double tmpDist;
				int common = intersectionArr[thread_id][j];
				int size0, size1;
				if(use64){
					size0 = sketches[i].hash64_arr.size();
					size1 = sketches[j].hash64_arr.size();
				}
				else{
					size0 = sketches[i].hash32_arr.size();
					size1 = sketches[j].hash32_arr.size();
				}
				if(!isContainment){
					int denom = size0 + size1 - common;
					double jaccard;
					if(size0 == 0 || size1 == 0)
						jaccard = 0.0;
					else
						jaccard = (double)common / denom;
					double mashD;
					if(jaccard == 1.0)
						mashD = 0.0;
					else if(jaccard == 0.0)
						mashD = 1.0;
					else
						mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
					tmpDist = mashD;
				}
				else{
					int denom = std::min(size0, size1);
					double containment;
					if(size0 == 0 || size1 == 0)
						containment = 0.0;
					else
						containment = (double)common / denom;
					double AafD;
					if(containment == 1.0)
						AafD = 0.0;
					else if(containment == 0.0)
						AafD = 1.0;
					else
						AafD = (double)-1.0 / kmer_size * log(containment);
					tmpDist = AafD;
				}

				for(int t = 0; t < denseSpan; t++){
					if(tmpDist <= distRadius[t]){
						denseLocalArr[t * threads + thread_id][i]++;
						denseLocalArr[t * threads + thread_id][j]++;
					}
				}
				double tmpANI = 1.0 - tmpDist;
				int ANI = (int)(tmpANI / 0.01);
				assert(ANI < 101);
				threadsANI[thread_id][ANI]++;
					
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

	for(int i = 0; i < 101; i++){
		for(int j = 0; j < threads; j++){
			aniArr[i] += threadsANI[j][i];
		}
	}
	for(int i = 0; i < threads; i++){
		delete(threadsANI[i]);
	}

	for(int i = 0; i < denseSpan; i++){
		for(int j = 0; j < threads; j++){
			for(int k = 0; k < N; k++){
				denseArr[i][k] += denseLocalArr[i*threads+j][k];
			}
		}
	}
	for(int i = 0; i < denseSpan * threads; i++){
		delete(denseLocalArr[i]);
	}

	if(tailNum != 0){
		for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
			if(use64){
				for(size_t j = 0; j < sketches[i].hash64_arr.size(); j++){
					uint64_t hash64 = sketches[i].hash64_arr[j];
					//if(!(dict[hash64/64] & (0x8000000000000000LLU >> (hash64 % 64))))	continue;
					if(hash_map_arr.count(hash64) == 0) continue;
					//for(auto x : hash_map_arr[hash64])
					for(size_t k = 0; k < hash_map_arr[hash64].size(); k++){
						size_t cur_index = hash_map_arr[hash64][k];
						intersectionArr[0][cur_index]++;
						//cerr << hash64 << '\t' << cur_index << endl;
					}
				}
			}
			else{
				for(size_t j = 0; j < sketches[i].hash32_arr.size(); j++){
					uint32_t hash = sketches[i].hash32_arr[j];
					if(sketchSizeArr[hash] == 0) continue;
					size_t start = hash > 0 ? offset[hash-1] : 0;
					size_t end = offset[hash];
					for(size_t k = start; k < end; k++){
						size_t curIndex = indexArr[k];
						intersectionArr[0][curIndex]++;
					}
				}
			}

			for(int j = i+1; j < sketches.size(); j++){
				//double tmpDist = 1.0 - minHashes[i].minHash->jaccard(minHashes[j].minHash);
				double tmpDist;
				int common = intersectionArr[0][j];
				int size0, size1;
				if(use64){
					size0 = sketches[i].hash64_arr.size();
					size1 = sketches[j].hash64_arr.size();
				}
				else{
					size0 = sketches[i].hash32_arr.size();
					size1 = sketches[j].hash32_arr.size();
				}
				if(!isContainment){
					int denom = size0 + size1 - common;
					double jaccard;
					if(size0 == 0 || size1 == 0)
						jaccard = 0.0;
					else
						jaccard = (double)common / denom;
					double mashD;
					if(jaccard == 1.0)
						mashD = 0.0;
					else if(jaccard == 0.0)
						mashD = 1.0;
					else
						mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
					tmpDist = mashD;
				}
				else{
					int denom = std::min(size0, size1);
					double containment;
					if(size0 == 0 || size1 == 0)
						containment = 0.0;
					else
						containment = (double)common / denom;
					double AafD;
					if(containment == 1.0)
						AafD = 0.0;
					else if(containment == 0.0)
						AafD = 1.0;
					else
						AafD = (double)-1.0 / kmer_size * log(containment);
					tmpDist = AafD;
				}

				for(int t = 0; t < denseSpan; t++){
					if(tmpDist <= distRadius[t]){
						denseArr[t][i]++;
						denseArr[t][j]++;
					}
				}
				double tmpANI = 1.0 - tmpDist;
				int ANI = (int)(tmpANI/0.01);
				assert(ANI < 101);
				aniArr[ANI]++;

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


vector<EdgeInfo> modifyMST(vector<SketchInfo>& sketches, int start_index, int sketch_func_id, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr){
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

				for(int t = 0; t < denseSpan; t++){
					if(tmpDist <= distRadius[t]){
						denseLocalArr[t * threads + thread_id][i]++;
						denseLocalArr[t * threads + thread_id][j]++;
					}
				}
				double tmpANI = 1.0 - tmpDist;
				int ANI = (int)(tmpANI / 0.01);
				assert(ANI < 101);
				threadsANI[thread_id][ANI]++;
					
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

	for(int i = 0; i < 101; i++){
		for(int j = 0; j < threads; j++){
			aniArr[i] += threadsANI[j][i];
		}
	}
	for(int i = 0; i < threads; i++){
		delete(threadsANI[i]);
	}

	for(int i = 0; i < denseSpan; i++){
		for(int j = 0; j < threads; j++){
			for(int k = 0; k < N; k++){
				denseArr[i][k] += denseLocalArr[i*threads+j][k];
			}
		}
	}
	for(int i = 0; i < denseSpan * threads; i++){
		delete(denseLocalArr[i]);
	}

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

				for(int t = 0; t < denseSpan; t++){
					if(tmpDist <= distRadius[t]){
						denseArr[t][i]++;
						denseArr[t][j]++;
					}
				}
				double tmpANI = 1.0 - tmpDist;
				int ANI = (int)(tmpANI/0.01);
				assert(ANI < 101);
				aniArr[ANI]++;

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

string build_newick_tree(vector<pair<int, double>>* G, bool* visited, const vector<SketchInfo>& sketches, bool sketch_by_file, int node){
	visited[node] = true;
	//if(G[node].size() == 0)	return to_string(node);
	string name;
	if(sketch_by_file) name = sketches[node].fileName;
	else name = sketches[node].seqInfo.name;
	if(G[node].size() == 0){
		return name;
	}
	string children("");
	for(auto x : G[node]){
		int cur_node = x.first;
		double dist = x.second;
		if(!visited[cur_node]){
			string child = build_newick_tree(G, visited, sketches, sketch_by_file, cur_node);
			children += child + ':' + to_string(dist) + ',';
		}
	}
	//if(children.length() == 0) return to_string(node);
	if(children.length() == 0) return name;
	//children = '(' + children.substr(0, children.length()-1) + ')' + to_string(node);
	children = '(' + children.substr(0, children.length()-1) + ')' + name;
	return children;
}


string get_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file){
	int vertice_number = sketches.size();
	vector<pair<int, double>>* G = new vector<pair<int, double>> [vertice_number];
	for(auto f : mst){
		int u = f.preNode;
		int v = f.sufNode;
		double dist = f.dist;
		G[u].push_back(make_pair(v, dist));
		G[v].push_back(make_pair(u, dist));
	}
	bool* visited = new bool[vertice_number];
	memset(visited, 0, vertice_number*sizeof(bool));

	string newick_tree = build_newick_tree(G, visited, sketches, sketch_by_file, 0);
	newick_tree += ';';

	return newick_tree;
}




