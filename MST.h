#ifndef H_MST_GRAPH
#define H_MST_GRAPH

#include <iostream>
#include <vector>
#include <unordered_set>

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
	EdgeInfo(int n1, int n2, double d){
		preNode = n1;
		sufNode = n2;
		dist = d;
	}

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

void kruskalMST(MST &mst, std::vector<EdgeInfo> graph, int numberNode);

void primMST(MST &mst, std::vector<Graph> & graphs, int numberNode);
	

void mst2Graph(MST &mst, std::vector<Graph> &graphs);

void mst2Graph(MST &mst, std::vector<Graph> &graphs, int graphSize);

void printMST(MST mst);
void printMST(std::vector<EdgeInfo> mst);

void printGraph(std::vector<Graph> graphs);

void printGraph(std::vector<EdgeInfo> graph);

void creatForest(MST srcMST, MST & resForest, double threshhold);


//void creatClust(std::vector<Graph> graphs, std::vector< std::vector<int> > & clusters);
void creatClust(std::vector<Graph> graphs, std::vector< std::unordered_set<int> > & clusters);

void mergeMST(MST & resMST, MST tmpMST);

void mergeGraph(std::vector<Graph> & resGraph, std::vector<Graph> tmpGraph);




#endif
