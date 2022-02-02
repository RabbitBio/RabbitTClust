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

using namespace std;

struct PosNum{
	int startPos;
	int clustSize;
};
struct IdPosNum{
	int label;
	int startPos;
	int clustSize;
};

inline void printInfo()
{
		cerr << "run with: ./calF1 RabbitTClust -l(-i) bacteria.out bacteria.f1" << endl;
		cerr << "The second argument (RabbitTClust) is applications, including RabbitTClust, MeshClust2, gclust or Mothur " << endl;
		cerr << "For the third argument, -l means genomes served as files, -i means genomes served as sequences" << endl;
		cerr << "The fourth argument (bacteria.out) is the cluster result from RabbitTClust, MeshClust2, gclust or Mothur " << endl;
		cerr << "The fifth argument (bacteria.f1) is the output file path" << endl;
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
	unordered_map<int, int> curMap;
	vector<IdPosNum> clustLabelInfo;
	int startPos = 0;

	while(getline(fs, line))
	{
		if(line.length() == 0) continue;
		if(application == "MeshClust2")
		{
			if(line[0] == '>')
			{
				if(curMap.size() != 0)
				{
					int maxNumber = 0;
					int maxIndex = 0;
					int clustSize = 0;

					for(auto x : curMap)
					{
						if(x.second > maxNumber)
						{
							maxNumber = x.second;
							maxIndex = x.first;
						}
						clustSize += x.second;
					}

					for(int i = 0; i < clustSize; i++)
					{
						ourClust.push_back(maxIndex);
					}

					IdPosNum ipn;
					ipn.label = maxIndex;
					ipn.startPos = startPos;
					ipn.clustSize = clustSize;
					clustLabelInfo.push_back(ipn);
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

				//for ground truth labels
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

				//for prediction labels
				int curLabel = standardMap[key];
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
					int maxNumber = 0;
					int maxIndex = 0;
					int clustSize = 0;

					for(auto x : curMap)
					{
						if(x.second > maxNumber)
						{
							maxNumber = x.second;
							maxIndex = x.first;
						}
						clustSize += x.second;
					}

					for(int i = 0; i < clustSize; i++)
					{
						ourClust.push_back(maxIndex);
					}

					IdPosNum ipn;
					ipn.label = maxIndex;
					ipn.startPos = startPos;
					ipn.clustSize = clustSize;
					clustLabelInfo.push_back(ipn);
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


				//for ground truth labels
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

				//for prediction labels
				int curLabel = standardMap[key];
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;

			}//end a cluster calculation
		}
	}//end while

	if(curMap.size() != 0)
	{
		int maxNumber = 0;
		int maxIndex = 0;
		int clustSize = 0;

		for(auto x : curMap)
		{
			if(x.second > maxNumber)
			{
				maxNumber = x.second;
				maxIndex = x.first;
			}
			clustSize += x.second;
		}

		for(int i = 0; i < clustSize; i++)
		{
			ourClust.push_back(maxIndex);
		}

		IdPosNum ipn;
		ipn.label = maxIndex;
		ipn.startPos = startPos;
		ipn.clustSize = clustSize;
		clustLabelInfo.push_back(ipn);
		startPos += clustSize;

		unordered_map<int, int>().swap(curMap);
	}

	int badLabel = -1;
	unordered_map<int, PosNum> uniqueMap;//used for store the representative cluster infos.
	for(int i = 0; i < clustLabelInfo.size(); i++)
	{
		int curLabel = clustLabelInfo[i].label;
		int curStartPos = clustLabelInfo[i].startPos;
		int curClustSize = clustLabelInfo[i].clustSize;
		if(uniqueMap.find(curLabel) == uniqueMap.end())//new label
		{
			PosNum pn;
			pn.startPos = curStartPos;
			pn.clustSize = curClustSize;
			uniqueMap.insert({curLabel, pn});
		}
		else//existed label
		{
			int preStartPos = uniqueMap[curLabel].startPos;
			int preClustSize = uniqueMap[curLabel].clustSize;
			if(curClustSize <= preClustSize)//curClust is less than preClust, change curClust, do not update uniqueMap
			{
				for(int j = curStartPos; j < curStartPos + curClustSize; j++)
					ourClust[j] = badLabel;
			}
			else//preClust is less than curClust, change preClust, update uniqueMap
			{
				for(int j = preStartPos; j < preStartPos + preClustSize; j++)
					ourClust[j] = badLabel;
				uniqueMap[curLabel].startPos = curStartPos;
				uniqueMap[curLabel].clustSize = curClustSize;
			}
			badLabel--;
		}
		
	}
	
	cerr << "the size of clustLabelInfo is: " << clustLabelInfo.size() << endl;
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
