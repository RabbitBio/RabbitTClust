/* Author: Xiaoming Xu 
 * Email: xiaoming.xu@mail.sdu.edu.cn
 *
 * calF1.cpp is used as preprocessing of the evaluation of precision, recall and F1-score.
 * The ground truth labels of genomes are as the first two keywords of nomenclature of gene feature.
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
	cerr << "run with: ./calF1 RabbitTClust -l(-i) bacteria.out bacteria.f1" << endl;
	cerr << "The second argument (RabbitTClust) is applications, including RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "For the third argument, -l means genomes served as files, -i means genomes served as sequences" << endl;
	cerr << "The fourth argument (bacteria.out) is the cluster result from RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "The fifth argument (bacteria.f1) is the output file path" << endl;
}

/* The output result is th resLabelArr with size of cluster number, each element is the label for the cluster.
 *
 *
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

void calF1(string application, string argument, string inputFile, string outputFile)
{	
	ofstream ofs(outputFile);
	ofstream ofs1(outputFile+".humanReadable");

	fstream fs(inputFile);
	string line;

	int curStandardIndex = 0;
	vector<int> ourClust;
	vector<int> standardClust;
	unordered_map<string, int> standardMap;
	unordered_map<int, int> groundTruthMap;
	unordered_map<int, int> curMap;


	int startPos = 0;
	vector< vector<LabNum> > labNumArr;
	vector<PosNum> posArr;

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

				if(type0.substr(0, 10) == "UNVERIFIED")
				{
					type0 = type1;
					type1 = type2;
				}
				if(type0.back() == ',') type0.pop_back();
				if(type1.back() == ',') type1.pop_back();

				string key = type0 + '\t' + type1;

				//for labels creating
				if(standardMap.find(key) == standardMap.end())//new label comes
				{
					standardMap.insert({key, curStandardIndex});
					standardClust.push_back(curStandardIndex);
					curStandardIndex++;
				}
				else
				{
					standardClust.push_back(standardMap[key]);
				}

				int curLabel = standardMap[key];
				//for ground truth labels
				groundTruthMap.insert({curLabel, 0});
				groundTruthMap[curLabel]++;

				//for prediction labels
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
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
						ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName >> type0 >> type1 >> type2;
					else if(argument == "-i")
						ss >> curId >> genomeId >> genomeSize >> genomeName >> type0 >> type1 >> type2;
				}
				else if(application == "Mothur")
					ss >> genomeName >> type0 >> type1 >> type2;
				else if(application == "gclust")
					ss >> curId >> genomeSize >> genomeName >> type0 >> type1 >> type2;
				else
				{
					cerr << "error application, need RabbitTClust, Mothur, gclust or MeshClust2" << endl;
					printInfo();
					return;
				}
				if(type0.substr(0, 10) == "UNVERIFIED")
				{
					type0 = type1;
					type1 = type2;
				}
				if(type0.back() == ',') type0.pop_back();
				if(type1.back() == ',') type1.pop_back();

				string key = type0 + '\t' + type1;


				//for labels creating
				if(standardMap.find(key) == standardMap.end())//new label comes
				{
					standardMap.insert({key, curStandardIndex});
					standardClust.push_back(curStandardIndex);
					curStandardIndex++;
				}
				else
				{
					standardClust.push_back(standardMap[key]);
				}

				int curLabel = standardMap[key];
				//for groundTruth labels
				groundTruthMap.insert({curLabel, 0});
				groundTruthMap[curLabel]++;

				//for prediction labels
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;

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

	//deal with the badLabels
	//get the remaining labels which are not the cluster labels.
	vector<IdNum> remainLabelArr;
	for(auto x : groundTruthMap)
	{
		int curLabel = x.first;
		if(globalMap.find(curLabel) == globalMap.end())//the label has been used
		{
			IdNum idn;
			idn.id = x.first;
			idn.number = x.second;
			remainLabelArr.push_back(idn);
		}
	}
	std::sort(remainLabelArr.begin(), remainLabelArr.end(), cmpIdNum);
	cerr << "the size of remainLabelArr is: " << remainLabelArr.size() << endl;

	vector<IdNum> badClustArr;
	for(int i = 0; i < resLabelArr.size(); i++)
	{
		IdNum idn;
		if(resLabelArr[i] < 0)
		{
			idn.id = i;
			idn.number = posArr[i].clustSize;
			badClustArr.push_back(idn);
		}
	}
	std::sort(badClustArr.begin(), badClustArr.end(), cmpIdNum);
	int remainLabelNumber = remainLabelArr.size();
	int badLabelNumber = badClustArr.size();
	int fixLabelNumber;
	if(remainLabelNumber < badLabelNumber)
	{
		fixLabelNumber = badLabelNumber - remainLabelNumber;
		for(int i = 0; i < remainLabelNumber; i++)
		{
			int curId = badClustArr[i].id;
			resLabelArr[curId] = remainLabelArr[i].id;
		}
	}
	else
	{
		fixLabelNumber = remainLabelNumber - badLabelNumber;
		for(int i = 0; i < badLabelNumber; i++)
		{
			int curId = badClustArr[i].id;
			resLabelArr[curId] = remainLabelArr[i].id;
		}
	}
	cerr << "the fixed number of label is: " << fixLabelNumber << endl;
	
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

//	cerr << "the badLabelNumber is: " << badLabelNumber << endl;
//	cerr << "the final badLabel is: " << badLabel << endl;
//	cerr << "the size of globalMap is: " << globalMap.size() << endl;
//	cerr << "the clust number is: " << clustNumber << endl;
//	cerr << "the standMap size is: " << standardMap.size() << endl;

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
	if(argc < 5){
		printInfo();
		return 1;
	}
	string application = argv[1];
	string argument = argv[2];
	string inputFile = argv[3];
	string outputFile = argv[4];

	calF1(application, argument, inputFile, outputFile);

	return 0;

}
