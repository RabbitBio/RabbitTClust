#include <iostream>
#include "MST.h"
#include <stdio.h>
#include <queue>
#include "UnionFind.h"
#include <omp.h>

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
	vector<EdgeInfo> mstArr[threads];
	int subSize = 8;
	int id = 0;
	int tailNum = sketches.size() % subSize;
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
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
						//tmpDist = 1.0 - sketches[i].minHash->jaccard(sketches[j].minHash);
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
		if(thread_id == 0 && id % 10000 < 50)	cerr << "finish MST: " << id << endl;

		sort(mstArr[thread_id].begin(), mstArr[thread_id].end(), cmpEdge);
		vector<EdgeInfo> tmpMst = kruskalAlgorithm(mstArr[thread_id], sketches.size());
		mstArr[thread_id].swap(tmpMst);
		vector<EdgeInfo>().swap(tmpMst);
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
						tmpDist = sketches[i].minHash->distance(sketches[j].minHash);
					//tmpDist = sketches[i].minHash->jaccard(sketches[j].minHash);
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

	vector<EdgeInfo> finalGraph;
	//for(int i = 0; i < mstArr.size(); i++){
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








