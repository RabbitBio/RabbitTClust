#ifdef GREEDY_CLUST
#include "greedy.h"
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
vector<vector<int> >greedyCluster(vector<SketchInfo> sketches, string sketchFunc, double threshold, int threads)
{
	int numGenomes = sketches.size();
	int * clustLabels = new int[numGenomes];
	memset(clustLabels, 0, numGenomes*sizeof(int));
	vector<vector<int> > cluster;

	for(int i = 0; i < numGenomes; i++){
		if(clustLabels[i] == 0){
			vector<int> curClust;
			curClust.push_back(i);
			clustLabels[i] = 1;
			for(int j = i+1; j < numGenomes; j++){
				if(clustLabels[j] == 0){
					double dist;
					if(sketchFunc == "MinHash"){
						if(sketches[i].isContainment)
							dist = 1.0 - sketches[i].minHash->containJaccard(sketches[j].minHash);
						else
							dist = sketches[i].minHash->distance(sketches[j].minHash);
					}
					else{
						cerr << "can only support MinHash with greedy incremental clust" << endl;
						exit(1);
					}
					if(dist < threshold){
						curClust.push_back(j);
						clustLabels[j] = 1;
					}
				}
			}
			cluster.push_back(curClust);
		}//end if
	}

	return cluster;
}

#endif