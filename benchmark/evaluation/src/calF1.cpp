/* Author: Xiaoming Xu 
 * Email: xiaoming.xu@mail.sdu.edu.cn
 * Data: 2022/2/18
 *
 * calF1.cpp is used as preprocessing of the evaluation of precision, recall and F1-score for bacteria, refseq, half-Bacteria, sub-Bacteria datasets.
 * The ground truth labels of genomes are as species_taxid which reveals nomenclature of gene feature.
 * The parameter -i and -l corresponding to the cluster of genomes served as sequences and files.
 * The input cluster files are in the CD-HIT format.
 *
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

using namespace std;
struct LabNum{
	int label;
	int number;
};

struct GlobalLabelInfo{
	int clustId;
	int labelNumber;
};

struct PosNum{
	int startPos;
	int clustSize;
};

struct IdNum{
	int id;
	int number;
};

bool cmpLabNum(LabNum ln1, LabNum ln2){
	return ln1.number > ln2.number;
}

bool cmpIdNum(IdNum in1, IdNum in2){
	return in1.number > in2.number;
}

inline void printInfo()
{
	cerr << "run with: ./calF1 RabbitTClust -l(-i) groundTruth bacteria.out bacteria.f1" << endl;
	cerr << "The second argument (RabbitTClust) is applications, including RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "For the third argument, -l means genomes served as files, -i means genomes served as sequences" << endl;
	cerr << "The fourth argument (groundTruth) is the ground truth from assembly_bacteria.txt of the <assembly_accession genomeName species_taxid> " << endl;
	cerr << "The fifth argument (bacteria.out) is the cluster result from RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "The sixth argument (bacteria.f1) is the output file path" << endl;
}

/* The output result is th resLabelArr with size of cluster number, each element is the label for the cluster.
 */
int updateLabel(vector< vector<LabNum> > &labNumArr, unordered_map<int, GlobalLabelInfo> &globalMap, int clustId, int &badLabel, vector<int> &resLabelArr)//return the new badLabel
{
	bool isBad = true;
	while(labNumArr[clustId].size() != 0 && isBad)
	{
		int curLabel = labNumArr[clustId][0].label;
		int curNumber = labNumArr[clustId][0].number;
		if(globalMap.find(curLabel) == globalMap.end())//new label 
		{
			GlobalLabelInfo glab;
			glab.clustId = clustId;
			glab.labelNumber = curNumber;
			globalMap.insert({curLabel, glab});
			resLabelArr[clustId] = curLabel;
			isBad = false;
		}
		else//label collison with previous cluster
		{
			int preClustId = globalMap[curLabel].clustId;
			int preNumber = globalMap[curLabel].labelNumber;
			if(curNumber > preNumber)//the previous cluster is defeated, need to update the previous cluster.
			{
				resLabelArr[clustId] = curLabel;
				isBad = false;
				globalMap[curLabel].clustId = clustId;
				globalMap[curLabel].labelNumber = curNumber;
				badLabel = updateLabel(labNumArr, globalMap, preClustId, badLabel, resLabelArr);
			}
			else//current cluster can not defeat the previous cluster, just erase the biggest label to try new biggest label
			{}
		}

		labNumArr[clustId].erase(labNumArr[clustId].begin());//erase the biggest label in this cluster
	}//end while
	if(isBad)
	{
		resLabelArr[clustId] = badLabel;
		badLabel--;//update the newBadLabel
	}
	return badLabel;
}

void calF1(string application, string argument, string groundTruth, string inputFile, string outputFile)
{	
	ofstream ofs(outputFile);
	ofstream ofs1(outputFile+".humanReadable");

	fstream fs0(groundTruth);
	string line;

	unordered_map<string, int> groundTruthMapFile;
	unordered_map<string, int> groundTruthMapSeq;

	while(getline(fs0, line))
	{
		string assembly_accession, genomeName, species_taxid;
		stringstream ss;
		ss << line;
		ss >> assembly_accession >> genomeName >> species_taxid;
		groundTruthMapFile.insert({assembly_accession, stoi(species_taxid)});
		groundTruthMapSeq.insert({genomeName, stoi(species_taxid)});

	}

	
	fstream fs(inputFile);

	int curStandardIndex = 0;
	vector<int> ourClust;
	vector<int> standardClust;
	unordered_map<string, int> standardMap;
	unordered_map<int, int> curMap;

	int startPos = 0;
	vector< vector<LabNum> > labNumArr;
	vector<PosNum> posArr;

	int numNotIngroundTruth = 0;

	while(getline(fs, line))
	{
		if(line.length() == 0) continue;
		if(application == "MeshClust2")
		{
			if(line[0] == '>')
			{
				if(curMap.size() != 0)
				{
					int clustSize = 0;
					vector<LabNum> curClustInfo;

					for(auto x : curMap)
					{
						LabNum ln;
						ln.label = x.first;
						ln.number = x.second;
						curClustInfo.push_back(ln);
						clustSize += x.second;
					}
					std::sort(curClustInfo.begin(), curClustInfo.end(), cmpLabNum);
					labNumArr.push_back(curClustInfo);

					PosNum pn;
					pn.startPos = startPos;
					pn.clustSize = clustSize;
					posArr.push_back(pn);

					startPos += clustSize;

					unordered_map<int, int>().swap(curMap);
				}
			}
			else{
				stringstream ss;
				ss << line;
				int curId, genomeId;
				string genomeSize, fileName, genomeName;
				string type0, type1, type2;
				if(argument == "-l")
					ss >> curId >> fileName >>genomeSize >> genomeName >> type0 >> type1 >> type2;
				else if(argument == "-i")
					ss >> curId >>genomeSize >> genomeName >> type0 >> type1 >> type2;
				else
				{
					cerr << "error argument, need -l or -i " << endl;
					printInfo();
					return;
				}

				genomeName = genomeName.substr(1);
				if(groundTruthMapSeq.find(genomeName) == groundTruthMapSeq.end())
				{
					//cerr << "the genomeName: " << genomeName << " is not in the groundTruthMapSeq!" << endl;
					numNotIngroundTruth++;
					continue;
				}
				else
				{
					int curLabel = groundTruthMapSeq[genomeName];
					standardClust.push_back(curLabel);
					curMap.insert({curLabel, 0});
					curMap[curLabel]++;
				}

			}
		}//end MeshClust2
		else //other application(RabbitTClust, gclust, Mothur)
		{
			if(line[0] != '\t')
			{
				if(curMap.size() != 0)
				{
					int clustSize = 0;
					vector<LabNum> curClustInfo;

					for(auto x : curMap)
					{
						LabNum ln;
						ln.label = x.first;
						ln.number = x.second;
						curClustInfo.push_back(ln);
						clustSize += x.second;
					}
					std::sort(curClustInfo.begin(), curClustInfo.end(), cmpLabNum);
					labNumArr.push_back(curClustInfo);

					PosNum pn;
					pn.startPos = startPos;
					pn.clustSize = clustSize;
					posArr.push_back(pn);

					startPos += clustSize;

					unordered_map<int, int>().swap(curMap);
				}
			}
			else{
				stringstream ss;
				ss << line;
				int curId, genomeId;
				string genomeSize, fileName, genomeName;
				string type0, type1, type2;
				if(application == "RabbitTClust")
				{
					if(argument == "-l")
					{
						ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName >> type0 >> type1 >> type2;
						int startIndex = fileName.find_last_of('/');
						int endIndex = fileName.find('_', startIndex + 5);
						string key = fileName.substr(startIndex+1, endIndex -startIndex -1);
						if(groundTruthMapFile.find(key) == groundTruthMapFile.end())
						{
							//cerr << "the key: " << key << " is not in the groundTruth!" << endl;
							numNotIngroundTruth++;
							continue;//skip this label
						}
						else
						{
							int curLabel = groundTruthMapFile[key];
							standardClust.push_back(curLabel);
							curMap.insert({curLabel, 0});
							curMap[curLabel]++;
						}
					}
					else if(argument == "-i")
					{
						ss >> curId >> genomeId >> genomeSize >> genomeName >> type0 >> type1 >> type2;
						if(groundTruthMapSeq.find(genomeName) == groundTruthMapSeq.end())
						{
							//cerr << "the genomeName: " << genomeName << " is not in the groundTruthMapSeq!" << endl;
							numNotIngroundTruth++;
							continue;
						}
						else
						{
							int curLabel = groundTruthMapSeq[genomeName];
							standardClust.push_back(curLabel);
							curMap.insert({curLabel, 0});
							curMap[curLabel]++;
						}
					}
				}
				else if(application == "Mothur")//TODO
				{
					ss >> genomeName >> type0 >> type1 >> type2;
					if(groundTruthMapSeq.find(genomeName) == groundTruthMapSeq.end())
					{
						//cerr << "the genomeName: " << genomeName << " is not in the groundTruthMapSeq!" << endl;
						numNotIngroundTruth++;
						continue;
					}
					else
					{
						int curLabel = groundTruthMapSeq[genomeName];
						standardClust.push_back(curLabel);
						curMap.insert({curLabel, 0});
						curMap[curLabel]++;
					}
				}
				else if(application == "gclust")//TODO
				{
					ss >> curId >> genomeSize >> genomeName >> type0 >> type1 >> type2;
					if(groundTruthMapSeq.find(genomeName) == groundTruthMapSeq.end())
					{
						//cerr << "the genomeName: " << genomeName << " is not in the groundTruthMapSeq!" << endl;
						numNotIngroundTruth++;
						continue;
					}
					else
					{
						int curLabel = groundTruthMapSeq[genomeName];
						standardClust.push_back(curLabel);
						curMap.insert({curLabel, 0});
						curMap[curLabel]++;
					}
				}
				else
				{
					cerr << "error application, need RabbitTClust, Mothur, gclust or MeshClust2" << endl;
					printInfo();
					return;
				}

			}//end a cluster calculation
		}
	}//end while

	if(curMap.size() != 0)
	{
		int clustSize = 0;
		vector<LabNum> curClustInfo;

		for(auto x : curMap)
		{
			LabNum ln;
			ln.label = x.first;
			ln.number = x.second;
			curClustInfo.push_back(ln);
			clustSize += x.second;
		}
		std::sort(curClustInfo.begin(), curClustInfo.end(), cmpLabNum);
		labNumArr.push_back(curClustInfo);

		PosNum pn;
		pn.startPos = startPos;
		pn.clustSize = clustSize;
		posArr.push_back(pn);

		startPos += clustSize;

		unordered_map<int, int>().swap(curMap);
	}

	//update the labels
	unordered_map<int, GlobalLabelInfo> globalMap;
	int badLabel = -1;
	vector<int> resLabelArr;
	int clustNumber = labNumArr.size();
	resLabelArr.resize(clustNumber); 
	for(int i = 0; i < clustNumber; i++)
	{
		badLabel = updateLabel(labNumArr, globalMap, i, badLabel, resLabelArr); 
	}
	
	//generate the result
	for(int i = 0; i < posArr.size(); i++)
	{
		int startPos = posArr[i].startPos;
		int clustSize = posArr[i].clustSize;
		for(int j = 0; j < clustSize; j++)
		{
			ourClust.push_back(resLabelArr[i]);
		}
	}

	cerr << "the number of which not in the groundTruth is: " << numNotIngroundTruth << endl;
	cerr << "the size of ourClust is: " << ourClust.size() << endl;
	cerr << "the size of standardClust is: " << standardClust.size() << endl;
	
	if(ourClust.size() != standardClust.size())
	{
		cerr << "the size of ourClust is not equal to the standardClust, exit()" << endl;
		return;
	}
	for(int i = 0; i < ourClust.size(); i++)
	{
		ofs1 << ourClust[i] << '\t' << standardClust[i] << endl;
	}

	for(int i = 0; i < ourClust.size(); i++)
		ofs << ourClust[i] << ' ';
	ofs << endl;
	
	for(int i = 0; i < standardClust.size(); i++)
		ofs << standardClust[i] << ' ';
	ofs << endl;

}

int main(int argc, char* argv[]){
	if(argc < 6){
		printInfo();
		return 1;
	}
	string application = argv[1];
	string argument = argv[2];
	string groundTruth = argv[3];
	string inputFile = argv[4];
	string outputFile = argv[5];

	calF1(application, argument, groundTruth, inputFile, outputFile);

	return 0;

}
