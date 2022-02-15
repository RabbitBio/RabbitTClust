/* Author: Xiaoming Xu
 * Email: xiaoming.xu@mail.sdu.edu.cn
 * 
 * calViral.cpp is used for evaluating the cluster result of viral datasets in viruSITE.
 * The ground truth is from the taxonomy tree from the viruSITE website (http://www.virusite.org/)
 * The parameter -i and -l corresponding to the cluster of genomes served as sequences and files.
 * The input cluster files are in the CD-HIT format.
 * The output files includes both the F1-score and NMI result files.
 *
 */

#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>

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
	cerr << "run with: ./calViral RabbitTClust -l(-i) groundTruth viral.out result_viral" << endl;
	cerr << "The second argument (RabbitTClust) is applications, including RabbitTClust, MeshClust2, MeshClust3(latest version, different output format from MeshClust2), gclust, uclust or Mothur" << endl;
	cerr << "For the third argument, -l means genomes served as files, -i means genomes served as sequences" << endl;
	cerr << "The fourth argument (groundTruth) is the ground truth file from the taxonomy tree " << endl;
	cerr << "The fifth argument (viral.out) is the cluster result from RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "The sixth argument (result_viral) is the output file path" << endl;
}

void calViral(string application, string argument, string groundTruth, string inputFile, string outputFile)
{
	string outputF1 = outputFile + ".f1";
	string outputNMI = outputFile + ".nmi";
	ofstream ofs0_f1(outputF1);
	ofstream ofs1_f1(outputF1 + ".humanReadable");
	ofstream ofs0_NMI(outputNMI);
	ofstream ofs1_NMI(outputNMI + ".humanReadable");

	fstream fs0(groundTruth);
	fstream fs1(inputFile);

	string line;
	
	unordered_map<string, int> groundTruthMap;
	int labelIndex = 0;
	int lineIndex = 0;

	vector<int> ourClustF1;
	vector<int> standardClust;
	unordered_map<int, int> curMap;
	vector<IdPosNum> clustLabelInfo;
	int startPos = 0;
	
	vector<int> ourClustNMI;
	int curClustIndex = 0;

	//for ground truth 
	while(getline(fs0, line))
	{
		if(line[0] == '\t')
		{
			if(lineIndex % 3 == 0)
			{
				stringstream ss;
				ss << line;
				string genomeName;
				ss >> genomeName;
				groundTruthMap.insert({genomeName, labelIndex});
			}
			lineIndex++;
		}
		else//a new label
		{
			labelIndex++;
		}

	}//end while
	cerr << "the size of groundTruthMap is: " << groundTruthMap.size() << endl;

	//for cluster result
	if(application == "MeshClust3")
	{
		int curId;
		string genomeSize, genomeName, fileName;
		while(getline(fs1, line))
		{
			if(line.length() == 0)//finish a cluster
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
						ourClustF1.push_back(maxIndex);
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
			else
			{
				stringstream ss;
				ss << line;
				ss >> curId >> genomeName;
				ourClustNMI.push_back(curId);
				genomeName = genomeName.substr(1);//remove the '>' in the first char of genomeName
				if(groundTruthMap.find(genomeName) == groundTruthMap.end())
				{
					cerr << "error label which is not in groundTruthMap " << endl;
					cerr << "the error label is: " << genomeName << endl;
					return;
				}
				standardClust.push_back(groundTruthMap[genomeName]);

				int curLabel = groundTruthMap[genomeName];
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
			}
		}//end while
	}
	else//other applications(MeshClust2, RabbitTClust, gclust, Mothur)
	{
		while(getline(fs1, line))
		{
			if(line.length() == 0) continue;
			if(application == "MeshClust2" || application == "gclust" || application == "uclust")
			{
				if(line[0] == '>')
				{
					curClustIndex++;
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
							ourClustF1.push_back(maxIndex);
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
					ourClustNMI.push_back(curClustIndex);
					stringstream ss;
					ss << line;
					int curId, genomeId;
					string genomeSize, fileName, genomeName;
					//string type0, type1, type2;
					if(application == "MeshClust2")
					{
						if(argument == "-l")
						 	ss >> curId >> fileName >> genomeSize >> genomeName;
						else if(argument == "-i")
							ss >> curId >> genomeSize >> genomeName;
						genomeName = genomeName.substr(1);//remove the '>' in the first char of genomeName
					}
					else if(application == "gclust")
					{
						ss >> curId >> genomeSize >> genomeName;
						genomeName = genomeName.substr(1);//remove the '>' in the first char of genomeName
					}
					else if(application == "uclust")
					{
						ss >> genomeName;
					}
					else
					{
						cerr << "error argument, need -l or -i " << endl;
						printInfo();
						return;
					}
					
					if(groundTruthMap.find(genomeName) == groundTruthMap.end())
					{
						cerr << "error label which is not in groundTruthMap " << endl;
						cerr << "the error label is: " << genomeName << endl;
						return;
					}

					//for ground truth labels
					standardClust.push_back(groundTruthMap[genomeName]);

					//for prediction labels
					int curLabel = groundTruthMap[genomeName];
					curMap.insert({curLabel, 0});
					curMap[curLabel]++;
						
				}
			}//end MeshClust2
			else//other applications(RabbitTClust, Mothur)
			{
				if(line[0] != '\t')
				{
					curClustIndex++;
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
							ourClustF1.push_back(maxIndex);
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
				else
				{
					ourClustNMI.push_back(curClustIndex);
					stringstream ss;
					ss << line;
					int curId, genomeId;
					string genomeSize, fileName, genomeName;
					//string type0, type1, type2;
					if(application == "RabbitTClust")
					{
						if(argument == "-l")
							ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName;
						else if(argument == "-i")
							ss >> curId >> genomeId >> genomeSize >> genomeName;
						else
						{
							cerr << "error argument, need -l or -i " << endl;
							printInfo();
							return;
						}
					}
					else if(application == "Mothur")
						ss >> genomeName;
					//else if(application == "gclust")
					//	ss >> curId >> genomeSize >> genomeName;
					else
					{
						cerr << "error application, need RabbitTClust, Mothur, glcust or MeshClust2" << endl;
						printInfo();
						return;
					}
					if(groundTruthMap.find(genomeName) == groundTruthMap.end())
					{
						cerr << "error label which is not in groundTruth " << endl;
						cerr << "the error label is: " << genomeName << endl;
						return;
					}

					//for ground truth labels
					standardClust.push_back(groundTruthMap[genomeName]);

					//for prediction labels
					int curLabel = groundTruthMap[genomeName];
					curMap.insert({curLabel, 0});
					curMap[curLabel]++;

				}
			}
		}//end while
	}

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
			ourClustF1.push_back(maxIndex);
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
					ourClustF1[j] = badLabel;
			}
			else//preClust is less than curClust, change preClust, update uniqueMap
			{
				for(int j = preStartPos; j < preStartPos + preClustSize; j++)
					ourClustF1[j] = badLabel;
				uniqueMap[curLabel].startPos = curStartPos;
				uniqueMap[curLabel].clustSize = curClustSize;
			}
			badLabel--;
		}
		
	}
	
	cerr << "the size of clustLabelInfo is: " << clustLabelInfo.size() << endl;
	cerr << "the size of ourClustF1 is: " << ourClustF1.size() << endl;
	cerr << "the size of standardClust is: " << standardClust.size() << endl;

	if(ourClustF1.size() != standardClust.size())
	{
		cerr << "the size of ourClustF1 is not equal to the standardClust, exit()" << endl;
		return;
	}
	for(int i = 0; i < ourClustF1.size(); i++)
	{
		ofs1_f1 << ourClustF1[i] << '\t' << standardClust[i] << endl;
		ofs1_NMI << ourClustNMI[i] << '\t' << standardClust[i] << endl;
	}

	for(int i = 0; i < ourClustF1.size(); i++)
	{
		ofs0_f1 << ourClustF1[i] << ' ';
		ofs0_NMI << ourClustNMI[i] << ' ';
	}
	ofs0_f1 << endl;
	ofs0_NMI << endl;
	
	for(int i = 0; i < standardClust.size(); i++)
	{
		ofs0_f1 << standardClust[i] << ' ';
		ofs0_NMI << standardClust[i] << ' ';
	}
	ofs0_f1 << endl;
	ofs0_NMI << endl;

}

int main(int argc, char * argv[]){
	if(argc < 6){
		printInfo();	
		return 1;
	}

	string application = argv[1];
	string argument = argv[2];
	string groundTruth = argv[3];
	string inputFile = argv[4];
	string outputFile = argv[5];

	calViral(application, argument, groundTruth, inputFile, outputFile);

	return 0;
}
