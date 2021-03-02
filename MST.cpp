#include <iostream>
#include "MST.h"
#include <stdio.h>
#include <queue>

using namespace std;

bool cmpNeighbor(NeighborNode n1, NeighborNode n2){
	return n1.distance < n2.distance;
}

void creatMST(MST &mst, vector <Graph> &graphs, int numberNode){
	if(graphs.size() == 0){
		cerr << "the graphs is NULL " << endl;
		exit(1);
	}
	mst.nodes.insert(graphs[0].node);
	int numMSTNode = 1;
	while(numMSTNode < numberNode){
		bool findMinEdge = false;

		double tmpMinDist = 1.0;
		int tmpPreNode = graphs[0].node;
		int tmpSufNode = graphs[0].node;

		//when the loop ends, the shortest edge combining the two vertex sets(MST and G-MST) is added into MST.
		for(int i = 0; i < graphs.size(); i++){
		//for(int i = 0; i < graphs.capacity(); i++){
			//the edge to add into the MST in prime algorithm
			//get the minimum edge from the edge of all the preNode in the MST.
			
			if(mst.nodes.count(graphs[i].node) > 0){//the preNode is in the MST.
				for(int j = 0; j < graphs[i].neighbor.size(); j++){
				//for(int j = graphs[i].curNeighbor; j < graphs[i].neighbor.size(); j++){
					//graphs[i].curNeighbor++;
					if(mst.nodes.count(graphs[i].neighbor[j].id) == 0 && tmpMinDist > graphs[i].neighbor[j].distance){
						tmpPreNode = graphs[i].node;
						tmpSufNode = graphs[i].neighbor[j].id;
						tmpMinDist = graphs[i].neighbor[j].distance;
						findMinEdge = true;
						break;//all the neighbor nodes(suffix nodes of the current node) have been sorted, the first one which is satisfied with "count(neighbor)==0" is the shortest one whose prefix node is the current node.

					}
				}
			}//end if: the preNode is in the MST.

		}//end for: the shortest edge is added.


		if(!findMinEdge){
			break;//end while //there are all circles when adding a new edge in the graph(sub-graph).
		}

		//add the shortest edge into the MST.
		mst.nodes.insert(tmpSufNode);
		//cout << tmpSufNode << endl;
		mst.edges.push_back(EdgeInfo(tmpPreNode, tmpSufNode, tmpMinDist));
		numMSTNode++;

	}//end while


	//if the graph is sub-graph, add other edge are all 1.0 to node(0).
	if(numMSTNode < numberNode){
		for(int i = 0; i < numberNode; i++){
			if(mst.nodes.count(i) == 0){
				mst.nodes.insert(i);
				mst.edges.push_back(EdgeInfo(graphs[0].node, i, 1.0));
			}
		}
	
	}//end if(numMSTNode < numberNode)

}


void mst2Graph(MST &mst, std::vector<Graph> &graphs){
	graphs.resize(mst.nodes.size());
	for(int i = 0; i < mst.nodes.size(); i++){
		graphs[i].node = i;
	}
	for(int i = 0; i < mst.edges.size(); i++){
		if(mst.edges[i].preNode >= mst.nodes.size()){
			printf("the node is out of bound, exit\n");
			return;
		}
		//graphs[mst.edges[i].preNode].node = mst.edges[i].preNode;
		graphs[mst.edges[i].preNode].neighbor.push_back(NeighborNode(mst.edges[i].sufNode, mst.edges[i].dist));
	}
}

void printMST(MST mst){
	for(int i = 0; i < mst.edges.size(); i++){
		printf("<%d, %d, %lf>\n", mst.edges[i].preNode, mst.edges[i].sufNode, mst.edges[i].dist);
	}

}

void printGraph(vector<Graph> graphs){
	for(int i = 0; i < graphs.size(); i++){
		for(int j = 0; j < graphs[i].neighbor.size(); j++){
			printf("<%d, %d, %lf> \n", graphs[i].node, graphs[i].neighbor[j].id, graphs[i].neighbor[j].distance);
		}
		cout << endl;

	}//end for i

}


void creatForest(MST srcMST, MST & resForest, double threshold){	
	resForest.nodes = srcMST.nodes;
	for(int i = 0; i < srcMST.edges.size(); i++){
		if(srcMST.edges[i].dist < threshold){
			resForest.edges.push_back(srcMST.edges[i]);
			//resForest.nodes.insert(srcMST.edges[i].preNode);
			//resForest.nodes.insert(srcMST.edges[i].sufNode);
		}

	}

}

//void creatClust(vector<Graph> graphs, vector< vector<int> > & clusters){
void creatClust(vector<Graph> graphs, vector< unordered_set<int> > & clusters){
	queue<int> q;
	unordered_set<int> traversedNode;//the head node.
	unordered_set<int> tmpCluster;
	//vector<int> tmpCluster;

	if(graphs.size() == 0){
		cerr << "the forest graph is NULL, cannot create cluster" << endl;
		return;
	}
	q.push(graphs[0].node);

	int index = 0;
	while(traversedNode.size() != graphs.size()){
		//cerr << traversedNode.size() << endl;
		if(q.empty()){
			//cerr << "the q is empty!" << endl;
			clusters.push_back(tmpCluster);
			tmpCluster.clear();
			for(int i = 0; i < graphs.size(); i++){
				if(traversedNode.count(graphs[i].node) == 0){
					q.push(graphs[i].node);	
					break;
				}
			}//end for

		}
		else{
			int curNode = q.front();
			//cerr << "the curNode is: " << curNode << endl;
			traversedNode.insert(curNode);
			tmpCluster.insert(curNode);
			//tmpCluster.push_back(curNode);
			q.pop();
			for(int i = 0; i < graphs.size(); i++){
				//if(traversedNode.count(graphs[i].node) == 0 && graphs[i].node == curNode)
				if(graphs[i].node == curNode){
					//tmpCluster.insert(curNode);
					for(int j = 0; j < graphs[i].neighbor.size(); j++){
						tmpCluster.insert(graphs[i].neighbor[j].id);
						//tmpCluster.push_back(graphs[i].neighbor[j].id);
						q.push(graphs[i].neighbor[j].id);
					}

					break;//only find one element in the query.

				}//end the one point and it's neighbor add into cluster.

			}//end for


		}//end else
		
	//	index++;
	//	if(index == 10) break;

	}//end while

}


void mergeGraph(vector<Graph> & resGraph, vector<Graph> tmpGraph){
	for(int i = 0; i < tmpGraph.size(); i++){
		for(int j = 0; j < resGraph.size(); j++){
			if(tmpGraph[i].node == resGraph[j].node){
				for(int k = 0; k < tmpGraph[i].neighbor.size(); k++){
					for(int l = 0; l < resGraph[j].neighbor.size(); l++){
						if(tmpGraph[i].neighbor[k].id == resGraph[j].neighbor[l].id){
							resGraph[j].neighbor[l].distance = tmpGraph[i].neighbor[k].distance < resGraph[j].neighbor[l].distance ? tmpGraph[i].neighbor[k].distance : resGraph[j].neighbor[l].distance;
							break;//end the tmpGraph[i].neighbor[k];
							
						}
						if(l == resGraph[j].neighbor.size()-1){
							resGraph[j].neighbor.push_back(tmpGraph[i].neighbor[k]);
						}
					}
				}//end the tmpGraph[i].neighbor[k];
			}//end if
		}
		
	}


}




		
