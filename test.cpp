#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include "MST.h"
#include <algorithm> //std::sort
#include <sys/time.h>

using namespace std;

double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}


int main(int argc, char* argv[]){
	
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	double t0 = get_sec();
	int n = 16;
	int N = n * (n-1) / 2;
	double *distArr = (double *)malloc(N *sizeof(double));

	for(int i = 0; i < N; i++){
		distArr[i] = distribution(generator);
	}
	
	double t1 = get_sec();


	vector<Graph> graphs;
	vector<MST> msts;

	int index = 0;
	for(int i = 0; i < n; i++){
		Graph tmpG;	
		tmpG.node = i;
		tmpG.curNeighbor = 0;
		for(int j = i+1; j < n; j++){

			tmpG.neighbor.push_back(NeighborNode(j, distArr[index]));
			//printf("<%d, %d, %lf> \n", i, j, distArr[index]);
			index++;
		}
		//sort the suffixNodes by weight of the edge in an ascending order.
		std::sort(tmpG.neighbor.begin(), tmpG.neighbor.end(), cmpNeighbor);

		graphs.push_back(tmpG);

	}
	cout << endl;
	cout << endl;
	cout << "print the Graph: " << endl;
	printGraph(graphs);
	cout << endl;
	cout << endl;


	MST resMst;
	//creatMST(resMst,midGraphs, n);
	primMST(resMst, graphs, n);



	cout << "print the MST: " << endl;
	for(int i = 0; i < resMst.edges.size(); i++){
		printf("<%d, %d, %lf>\n", resMst.edges[i].preNode, resMst.edges[i].sufNode, resMst.edges[i].dist);
	}

	cout << "output the graph================================" << endl;	

	double t2 = get_sec();

	cout << "start construct the MST================================" << endl;	
	

	double t3 = get_sec();
	//vector<Graph> resGraph;
	//mst2Graph(mst, resGraph);

	double t4 = get_sec();

	cout << "start print the final graph================================" << endl;	
	//for(int i = 0; i <resGraph.size(); i++){
	//	for(int j = 0; j <resGraph[i].neighbor.size(); j++){
	//		printf("<%d, %d, %lf> \n",resGraph[i].node,resGraph[i].neighbor[j].id,resGraph[i].neighbor[j].distance);
	//	}
	//	cout << endl;

	//}//end for i

		
	cerr << "the time of random array is: "  << t1 - t0 << endl;
	cerr << "the time of generate graph is: "  << t2 - t1 << endl;
	cerr << "the time of generate MST is: "  << t3 - t2 << endl;
	cerr << "the time of mst2Graph is: "  << t4 - t3 << endl;
		

		


			



	return 0;
}
