#ifdef GREEDY_CLUST
#include "greedy.h"
#include <map>
#include <omp.h>
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
				if(sketches[repId].isContainment)
					dist = sketches[repId].minHash->containDistance(sketches[j].minHash);
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

#endif
