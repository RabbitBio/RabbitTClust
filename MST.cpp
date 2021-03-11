#include <iostream>
#include "MST.h"
#include <stdio.h>
#include <queue>
#include "UnionFind.h"

using namespace std;

bool cmpInfo(SimilarityInfo s1, SimilarityInfo s2){
	return s1.id < s2.id;
}	

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
    return spanningTree;
}


vector<EdgeInfo> generateForest(vector <EdgeInfo> mst, double threshhold){
	vector<EdgeInfo> forest;
	for(int i = 0; i < mst.size(); i++){
		if(mst[i].dist < threshhold){
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




void updateForest(vector<MST> &forest, MST mst, int curSufNode){//update forest with the result MST;
	for(int i = 0; i < forest.size(); i++){
		if(forest[i].nodes.count(curSufNode) > 0){
			mst.nodes.insert(forest[i].nodes.begin(), forest[i].nodes.end());
			mst.edges.insert(mst.edges.end(), forest[i].edges.begin(), forest[i].edges.end());
			forest.erase(forest.begin() + i);
			break;
		}
	}

}

void updateForest(vector<MST> &forest, int curIndex, int curSufNode){//inner update of the forest.
	if(curIndex >= forest.size()){
		cerr << "error index of forest update" << endl;
		exit(1);
	}
	for(int i = 0; i < forest.size(); i++){
		if(i == curIndex){
			//current MST in the forest, continue;
			continue;
		}
		if(forest[i].nodes.count(curSufNode) > 0){
			//the MST can be merged into the curMST
			forest[curIndex].nodes.insert(forest[i].nodes.begin(), forest[i].nodes.end());
			forest[curIndex].edges.insert(forest[curIndex].edges.end(), forest[i].edges.begin(), forest[i].edges.end());
			forest.erase(forest.begin()+i);//remove the MST which has been merged.
			break;	//only one MST has the intesection node with the current MST on curSufNode.
					//if more than one MST have the intsection on curSufNode, they should have been merged before.	
		}
	}
	

}

void kruskalMST(MST &mst, vector<EdgeInfo> graph, int numberNode){
	mst.nodes.insert(graph[0].preNode);
	mst.nodes.insert(graph[0].sufNode);
	mst.edges.push_back(graph[0]);

	vector<MST> forest;

	for(int i = 1; i < graph.size(); i++){
		if(mst.nodes.count(graph[i].preNode) > 0 && mst.nodes.count(graph[i].sufNode) > 0){
			//both nodes of the edge have been added into the tree, which can cause circle.
			continue;
		}
		else if(mst.nodes.count(graph[i].preNode) == 0 && mst.nodes.count(graph[i].sufNode) == 0){
			//both nodes of the edge are not in the MST.
			//they should be in another MST;
			bool inForest = false;
			bool causeCircle = false;
			for(int i = 0; i < forest.size(); i++){
				if(forest[i].nodes.count(graph[i].preNode) > 0 && forest[i].nodes.count(graph[i].sufNode) > 0){
					//it can cause circle in this MST.
					//so it can cause circle in the final MST, this edge is omited.
					causeCircle = true;
					break;
				}
				else if(forest[i].nodes.count(graph[i].preNode) == 0 && forest[i].nodes.count(graph[i].sufNode) > 0){
					//the edge can be insert into this MST.
					forest[i].nodes.insert(graph[i].preNode);
					forest[i].edges.push_back(graph[i]);
					inForest = true;
					//update the forest to merged the connected MST.
					updateForest(forest, i, graph[i].preNode);
					break;
				}
				else if(forest[i].nodes.count(graph[i].preNode) > 0 && forest[i].nodes.count(graph[i].sufNode) == 0){
					forest[i].nodes.insert(graph[i].sufNode);
					forest[i].edges.push_back(graph[i]);
					inForest = true;
					updateForest(forest, i, graph[i].sufNode);
					break;
				}
				else{
					//the edge cannot add into this MST,find next MST.
					continue;
				}

			}//end for lookup forest.
			if(!inForest && !causeCircle){
				//the edge cannot cause circle and is not in the forest.
				MST newMST;
				newMST.nodes.insert(graph[i].preNode);
				newMST.nodes.insert(graph[i].sufNode);
				newMST.edges.push_back(graph[i]);
				forest.push_back(newMST);
			}

		}
		else if(mst.nodes.count(graph[i].preNode) == 0 && mst.nodes.count(graph[i].sufNode) > 0){
			//the edge can be insert in final MST.
			mst.nodes.insert(graph[i].preNode);
			mst.edges.push_back(graph[i]);
			//update the forest to merged the connected MST into final MST.
			updateForest(forest, mst, graph[i].preNode);
		}
		else{ //mst.node.count(graph[i].preNode > 0 && mst.nodes.count(graph[i].sufNode) == 0)
			mst.nodes.insert(graph[i].sufNode);
			mst.edges.push_back(graph[i]);
			updateForest(forest, mst, graph[i].sufNode);

		}
		if(mst.nodes.size()== numberNode){
			cerr << "the finished edge index is: " << i << endl;
			break;
		}
		
	}
	cerr << "the graph.size() is: " << graph.size() << endl;
	cerr << "the mst.node.size() is: " << mst.nodes.size() << endl;
	cerr << "the mst.edges.size() is: " << mst.edges.size() << endl;

}

void primMST(MST &mst, vector <Graph> &graphs, int numberNode){
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
		int tmpParent = 0;

		//when the loop ends, the shortest edge combining the two vertex sets(MST and G-MST) is added into MST.
		for(int i = 0; i < graphs.size(); i++){
		//for(int i = 0; i < graphs.capacity(); i++){
			//the edge to add into the MST in prime algorithm
			//get the minimum edge from the edge of all the preNode in the MST.
			
			if(mst.nodes.count(graphs[i].node) > 0){//the preNode is in the MST.
				//for(int j = 0; j < graphs[i].neighbor.size(); j++){
				for(int j = graphs[i].curNeighbor; j < graphs[i].neighbor.size(); j++){
					//graphs[i].curNeighbor++;
					if(mst.nodes.count(graphs[i].neighbor[j].id) == 0 && tmpMinDist > graphs[i].neighbor[j].distance){
						tmpPreNode = graphs[i].node;
						tmpSufNode = graphs[i].neighbor[j].id;
						tmpMinDist = graphs[i].neighbor[j].distance;
						findMinEdge = true;
						tmpParent = i;
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
		graphs[tmpParent].curNeighbor++;
		//cout << tmpSufNode << endl;
		cerr << "numMSTNode is: " << numMSTNode << " add SufNode: " << tmpSufNode << endl;
		EdgeInfo tmpE;
		tmpE.preNode = tmpPreNode;
		tmpE.sufNode = tmpSufNode;
		tmpE.dist = tmpMinDist;
		mst.edges.push_back(tmpE);
		numMSTNode++;

	}//end while


//	//if the graph is sub-graph, add other edge are all 1.0 to node(0).
//	if(numMSTNode < numberNode){
//		for(int i = 0; i < numberNode; i++){
//			if(mst.nodes.count(i) == 0){
//				mst.nodes.insert(i);
//				mst.edges.push_back(EdgeInfo(graphs[0].node, i, 1.0));
//			}
//		}
//	
//	}//end if(numMSTNode < numberNode)

}

void mst2Graph(MST &mst, std::vector<Graph> &graphs, int graphSize){
	graphs.resize(graphSize);
	for(int i = 0; i < graphSize; i++){
		graphs[i].node = i;
	}
	for(int i = 0; i < mst.edges.size(); i++){
		if(mst.edges[i].preNode >= graphSize){
			printf("the node is out of bound, exit\n");
			return;
		}
		//graphs[mst.edges[i].preNode].node = mst.edges[i].preNode;
		graphs[mst.edges[i].preNode].neighbor.push_back(NeighborNode(mst.edges[i].sufNode, mst.edges[i].dist));
	}
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

void printMST(vector<EdgeInfo> mst){
	for(int i = 0; i < mst.size(); i++){
		printf("<%d, %d, %lf>\n", mst[i].preNode, mst[i].sufNode, mst[i].dist);
	}
}

void printMST(MST mst){
	for(int i = 0; i < mst.edges.size(); i++){
		printf("<%d, %d, %lf>\n", mst.edges[i].preNode, mst.edges[i].sufNode, mst.edges[i].dist);
		//printf("<%d, %d, %d>\n", mst.edges[i].preNode, mst.edges[i].sufNode, (int)mst.edges[i].dist);
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

void printGraph(vector<EdgeInfo> graph){
	for(int i = 0; i < graph.size(); i++){
		printf("<%d, %d, %lf> \n", graph[i].preNode, graph[i].sufNode, graph[i].dist);
	}
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

//the mst tobe merged should have the same nodes.
void mergeMST(MST & resMST, MST tmpMST){ 
	for(int i = 0; i < tmpMST.edges.size(); i++){
		resMST.edges.push_back(tmpMST.edges[i]);
	}

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




		
