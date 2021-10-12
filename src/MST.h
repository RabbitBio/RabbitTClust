#ifndef H_MST_GRAPH
#define H_MST_GRAPH

#include <iostream>
#include <vector>
#include <unordered_set>
#include "SketchInfo.h"

struct NeighborNode{
	int id;
	double distance;
	NeighborNode(int i, double d){
		id = i;
		distance = d;
	}
};

struct EdgeInfo{
	int preNode;
	int sufNode;
	double dist;
};


struct Graph{
	int node;
	int curNeighbor;
	std::vector<NeighborNode> neighbor;

};


struct MST{
	std::unordered_set<int> nodes;
	std::vector<EdgeInfo> edges;

};

bool cmpEdge(EdgeInfo e1, EdgeInfo e2);

bool cmpNeighbor(NeighborNode n1, NeighborNode n2);

std::vector<EdgeInfo> kruskalAlgorithm(std::vector<EdgeInfo>graph, int vertices);

vector<EdgeInfo> generateMST(vector<SketchInfo>& sketches, string sketchFunc, int threads);

std::vector<EdgeInfo> generateForest(std::vector <EdgeInfo> mst, double threshhold);

std::vector <std::vector<int> > generateCluster(std::vector<EdgeInfo> forest, int vertices);

#endif
