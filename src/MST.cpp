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

vector<EdgeInfo> generateMST(vector<SketchInfo>& sketches, string sketchFunc, int threads){

	//for ANI distribution calculation.
	uint64_t aniArr[101];
	uint64_t  threadsANI[threads][101];
	for(int i = 0; i < 101; i++){
		aniArr[i] = 0;
		for(int j = 0; j < threads; j++){
			threadsANI[j][i] = 0;
		}
	}
	vector<EdgeInfo> mstArr[threads];
	int subSize = 8;
	int id = 0;
	int tailNum = sketches.size() % subSize;
	int N = sketches.size();
	uint64_t totalCompNum = (uint64_t)N * (uint64_t)(N-1)/2;
	uint64_t percentNum = totalCompNum / 100;
	cerr << "the percentNum is: " << percentNum << endl;
	uint64_t percentId = 0;
	#pragma omp parallel for num_threads(threads) schedule (dynamic)
	for(id = 0; id < sketches.size() - tailNum; id+=subSize){
		int thread_id = omp_get_thread_num();
		for(int i = id; i < id+subSize; i++){
			for(int j = i+1; j < sketches.size(); j++){
				double tmpDist;
				if(sketchFunc == "MinHash")
				{
					if(sketches[i].isContainment)
					{
						tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
						//double jaccard = sketches[i].minHash->containJaccard(sketches[j].minHash);
						//tmpDist = -log(2 * jaccard / (1.0 + jaccard)) / 21;
					}
					else
					{
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
						//tmpDist = 1.0 - sketches[i].minHash->jaccard(sketches[j].minHash);
						double tmpANI = 1.0 - tmpDist;
						int ANI = (int)(tmpANI / 0.01);
						assert(ANI < 101);
						threadsANI[thread_id][ANI]++;
					}
				}
				else if(sketchFunc == "WMH"){
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				}
				else if(sketchFunc == "HLL"){
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				}
				else if(sketchFunc == "OMH"){
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				}
				else	
					break;
					
				EdgeInfo tmpE{i, j, tmpDist};
				mstArr[thread_id].push_back(tmpE);
			}
		}
		//if(thread_id == 0 && id % 10000 < 50)	cerr << "finish MST: " << id << endl;
		if(thread_id == 0){
			uint64_t computedNum = (uint64_t)(N - id) * (uint64_t)id + (uint64_t)id * (uint64_t)id / 2;
			if(computedNum >= percentId * percentNum){
				fprintf(stderr, "finish MST generation %d %\n", percentId);
				percentId++;
			}
		}

		sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
		vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
		mstArr[thread_id].swap(tmpMst);
		vector<EdgeInfo>().swap(tmpMst);
	}

	for(int i = 0; i < 101; i++){
		for(int j = 0; j < threads; j++){
			aniArr[i] += threadsANI[j][i];
		}
	}
	if(tailNum != 0){
		for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
			for(int j = i+1; j < sketches.size(); j++){
				//double tmpDist = 1.0 - minHashes[i].minHash->jaccard(minHashes[j].minHash);
				double tmpDist;
				if(sketchFunc == "MinHash"){
					if(sketches[i].isContainment)
						tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
					else
					{
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
					//tmpDist = sketches[i].minHash->jaccard(sketches[j].minHash);
						double tmpANI = 1.0 - tmpDist;
						int ANI = (int)(tmpANI/0.01);
						assert(ANI < 101);
						aniArr[ANI]++;
					}
				}
				else if(sketchFunc == "WMH")
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				else if(sketchFunc == "HLL")
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				else if(sketchFunc == "OMH")
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				else	
					break;

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
	int sketchSize = sketches[0].minHash->getSketchSize();
	string outputFile = "ANI_calculation_" + to_string(sketchSize); 
	ofstream ofs(outputFile);
	for(int i = 0; i < 101; i++){
		ofs << i << '\t' << aniArr[i] << endl;
	}
	ofs.close();

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

vector<EdgeInfo> generateForest(vector <EdgeInfo> mst, double threshhold){
	vector<EdgeInfo> forest;
	for(int i = 0; i < mst.size(); i++){
		if(mst[i].dist <= threshhold){
			forest.push_back(mst[i]);
		}
	}
	return forest;
}

vector<EdgeInfo> modifyForest(vector<EdgeInfo> forest, vector<int> noiseArr){
	vector<EdgeInfo> newForest;
	bool * removeTag = new bool[forest.size()];
	#pragma omp parallel for num_threads(48)
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
		if(removeTag[i] == 0)
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

vector < vector<int> > generateCluster(vector<EdgeInfo> forest, int vertices){
	UnionFind uf(vertices);
	vector< vector<int> > cluster;
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

//vector<int> getNoiseNode(int * densArr, int size, double alpha){
vector<int> getNoiseNode(vector<PairInt> densPairArr, double alpha){
	vector<int> noiseArr;
	std::sort(densPairArr.begin(), densPairArr.end(), cmpPair);
	vector<int> denseSet;
	int curDense = densPairArr[0].second;
	denseSet.push_back(curDense);
	for(int i = 1; i < densPairArr.size(); i++){
		if(densPairArr[i].second - curDense > curDense * 0.01){
			curDense = densPairArr[i].second;
			denseSet.push_back(curDense);
		}
	}

	//for(int i = 0; i < densPairArr.size(); i++){
	//	cout << i << '\t' << densPairArr[i].first << '\t' << densPairArr[i].second << endl;
	//}
	//cout << "end the dence print" << endl;
	//int minValidDense = 0;
	//for(int i = 0; i < size; i++){
	//	if(densPairArr[i].second > 0){
	//		minValidDense = densPairArr[i].second;
	//		break;
	//	}
	//}
	int size = denseSet.size();
	if(size < 4) return noiseArr;
	int Q1_index = size / 4;
	int Q2_index = size / 4 * 2;
	int Q3_index = size / 4 * 3;
	int Q1 = denseSet[Q1_index];
	int Q2 = denseSet[Q2_index];
	int Q3 = denseSet[Q3_index];
	//int Q1 = densPairArr[Q1_index].second;
	//int Q3 = densPairArr[Q3_index].second;
	//int thresholdDense = std::max((int)(Q1 - alpha * (Q3 - Q1)), minValidDense);
	//int thresholdDense = std::max((int)(Q1 - alpha * (Q3 - Q1)), 0);
	int thresholdDense = std::max((int)(Q2 - alpha * (Q3 - Q2)), 0);
	cout << "the thresholdDense is: " << thresholdDense << ", the alpha is: " << alpha << endl;
	for(int i = 0; i < densPairArr.size(); i++){
		if(densPairArr[i].second < thresholdDense){
			noiseArr.push_back(densPairArr[i].first);
		}
		else{
			break;
		}
	}
	return noiseArr;
}

vector<EdgeInfo> modifyMST(vector<SketchInfo>& sketches, string sketchFunc, int threads, int* &dens1Arr, int* &dens2Arr, int* &dens3Arr, int* &dens4Arr, int* &dens5Arr, string prefixName){
	dens1Arr = new int[sketches.size()];
	dens2Arr = new int[sketches.size()];
	dens3Arr = new int[sketches.size()];
	dens4Arr = new int[sketches.size()];
	dens5Arr = new int[sketches.size()];

	memset(dens1Arr, 0, sketches.size() * sizeof(int));
	memset(dens2Arr, 0, sketches.size() * sizeof(int));
	memset(dens3Arr, 0, sketches.size() * sizeof(int));
	memset(dens4Arr, 0, sketches.size() * sizeof(int));
	memset(dens5Arr, 0, sketches.size() * sizeof(int));

	int* dens_local1Arr[threads];
	int* dens_local2Arr[threads];
	int* dens_local3Arr[threads];
	int* dens_local4Arr[threads];
	int* dens_local5Arr[threads];
	for(int i = 0; i < threads; i++){
		dens_local5Arr[i] = new int[sketches.size()];
		dens_local1Arr[i] = new int[sketches.size()];
		dens_local2Arr[i] = new int[sketches.size()];
		dens_local3Arr[i] = new int[sketches.size()];
		dens_local4Arr[i] = new int[sketches.size()];
		memset(dens_local1Arr[i], 0, sketches.size() * sizeof(int));
		memset(dens_local2Arr[i], 0, sketches.size() * sizeof(int));
		memset(dens_local3Arr[i], 0, sketches.size() * sizeof(int));
		memset(dens_local4Arr[i], 0, sketches.size() * sizeof(int));
		memset(dens_local5Arr[i], 0, sketches.size() * sizeof(int));
	}

	//for ANI distribution calculation.
	uint64_t aniArr[101];
	uint64_t  threadsANI[threads][101];
	for(int i = 0; i < 101; i++){
		aniArr[i] = 0;
		for(int j = 0; j < threads; j++){
			threadsANI[j][i] = 0;
		}
	}
	vector<EdgeInfo> mstArr[threads];
	int subSize = 8;
	int id = 0;
	int tailNum = sketches.size() % subSize;
	int N = sketches.size();
	uint64_t totalCompNum = (uint64_t)N * (uint64_t)(N-1)/2;
	uint64_t percentNum = totalCompNum / 100;
	cerr << "the percentNum is: " << percentNum << endl;
	uint64_t percentId = 0;
	#pragma omp parallel for num_threads(threads) schedule (dynamic)
	for(id = 0; id < sketches.size() - tailNum; id+=subSize){
		int thread_id = omp_get_thread_num();
		for(int i = id; i < id+subSize; i++){
			for(int j = i+1; j < sketches.size(); j++){
				double tmpDist;
				if(sketchFunc == "MinHash")
				{
					if(sketches[i].isContainment)
					{
						tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
						//double jaccard = sketches[i].minHash->containJaccard(sketches[j].minHash);
						//tmpDist = -log(2 * jaccard / (1.0 + jaccard)) / 21;
					}
					else
					{
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
						if(tmpDist <= 0.05){
							dens_local5Arr[thread_id][i]++;
							dens_local5Arr[thread_id][j]++;
						}
						if(tmpDist <= 0.04){
							dens_local4Arr[thread_id][i]++;
							dens_local4Arr[thread_id][j]++;
						}
						if(tmpDist <= 0.03){
							dens_local3Arr[thread_id][i]++;
							dens_local3Arr[thread_id][j]++;
						}
						if(tmpDist <= 0.02){
							dens_local2Arr[thread_id][i]++;
							dens_local2Arr[thread_id][j]++;
						}
						if(tmpDist <= 0.01){
							dens_local1Arr[thread_id][i]++;
							dens_local1Arr[thread_id][j]++;
						}
						//tmpDist = 1.0 - sketches[i].minHash->jaccard(sketches[j].minHash);
						double tmpANI = 1.0 - tmpDist;
						int ANI = (int)(tmpANI / 0.01);
						assert(ANI < 101);
						threadsANI[thread_id][ANI]++;
					}
				}
				else if(sketchFunc == "WMH"){
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				}
				else if(sketchFunc == "HLL"){
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				}
				else if(sketchFunc == "OMH"){
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				}
				else	
					break;
					
				EdgeInfo tmpE{i, j, tmpDist};
				mstArr[thread_id].push_back(tmpE);
			}
		}
		//if(thread_id == 0 && id % 10000 < 50)	cerr << "finish MST: " << id << endl;
		if(thread_id == 0){
			uint64_t computedNum = (uint64_t)(N - id) * (uint64_t)id + (uint64_t)id * (uint64_t)id / 2;
			if(computedNum >= percentId * percentNum){
				fprintf(stderr, "finish MST generation %d %\n", percentId);
				percentId++;
			}
		}

		sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
		vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
		mstArr[thread_id].swap(tmpMst);
		vector<EdgeInfo>().swap(tmpMst);
	}

	for(int i = 0; i < 101; i++){
		for(int j = 0; j < threads; j++){
			aniArr[i] += threadsANI[j][i];
		}
	}
	for(int i = 0; i < sketches.size(); i++){
		for(int j = 0; j < threads; j++){
			dens1Arr[i] += dens_local1Arr[j][i];
			dens2Arr[i] += dens_local2Arr[j][i];
			dens3Arr[i] += dens_local3Arr[j][i];
			dens4Arr[i] += dens_local4Arr[j][i];
			dens5Arr[i] += dens_local5Arr[j][i];
		}
	}
	for(int i = 0; i < threads; i++){
		delete(dens_local1Arr[i]);
		delete(dens_local2Arr[i]);
		delete(dens_local3Arr[i]);
		delete(dens_local4Arr[i]);
		delete(dens_local5Arr[i]);
	}

	if(tailNum != 0){
		for(int i = sketches.size()-tailNum; i < sketches.size(); i++){
			for(int j = i+1; j < sketches.size(); j++){
				//double tmpDist = 1.0 - minHashes[i].minHash->jaccard(minHashes[j].minHash);
				double tmpDist;
				if(sketchFunc == "MinHash"){
					if(sketches[i].isContainment)
						tmpDist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
					else
					{
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
						if(tmpDist <= 0.05){
							dens5Arr[i]++;
							dens5Arr[j]++;
						}
						if(tmpDist <= 0.04){
							dens4Arr[i]++;
							dens4Arr[j]++;
						}
						if(tmpDist <= 0.03){
							dens3Arr[i]++;
							dens3Arr[j]++;
						}
						if(tmpDist <= 0.02){
							dens2Arr[i]++;
							dens2Arr[j]++;
						}
						if(tmpDist <= 0.01){
							dens1Arr[i]++;
							dens1Arr[j]++;
						}
					//tmpDist = sketches[i].minHash->jaccard(sketches[j].minHash);
						double tmpANI = 1.0 - tmpDist;
						int ANI = (int)(tmpANI/0.01);
						assert(ANI < 101);
						aniArr[ANI]++;
					}
				}
				else if(sketchFunc == "WMH")
					tmpDist = sketches[i].WMinHash->distance(sketches[j].WMinHash);
				else if(sketchFunc == "HLL")
					tmpDist = sketches[i].HLL->distance(*sketches[j].HLL);
				else if(sketchFunc == "OMH")
					tmpDist = sketches[i].OMH->distance(*sketches[j].OMH);
				else	
					break;

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
	int sketchSize = sketches[0].minHash->getSketchSize();
	string outputFile = "ANI_calculation_" + prefixName + to_string(sketchSize); 
	ofstream ofs(outputFile);
	for(int i = 0; i < 101; i++){
		ofs << i << '\t' << aniArr[i] << endl;
	}
	ofs.close();

	string outputFile5 = "Dense_calculation5_" + prefixName + to_string(sketchSize); 
	string outputFile4 = "Dense_calculation4_" + prefixName + to_string(sketchSize); 
	string outputFile3 = "Dense_calculation3_" + prefixName + to_string(sketchSize); 
	string outputFile2 = "Dense_calculation2_" + prefixName + to_string(sketchSize); 
	string outputFile1 = "Dense_calculation1_" + prefixName + to_string(sketchSize); 
	ofstream ofs5(outputFile5);
	ofstream ofs4(outputFile4);
	ofstream ofs3(outputFile3);
	ofstream ofs2(outputFile2);
	ofstream ofs1(outputFile1);
	for(int i = 0; i < sketches.size(); i++){
		ofs5 << i << '\t' << dens5Arr[i] << endl;
		ofs4 << i << '\t' << dens4Arr[i] << endl;
		ofs3 << i << '\t' << dens3Arr[i] << endl;
		ofs2 << i << '\t' << dens2Arr[i] << endl;
		ofs1 << i << '\t' << dens1Arr[i] << endl;
	}
	ofs5.close();
	ofs4.close();
	ofs3.close();
	ofs2.close();
	ofs1.close();
			

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







