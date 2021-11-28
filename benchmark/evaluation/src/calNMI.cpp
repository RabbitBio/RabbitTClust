/* The evaluation of cluster is based on the precision and recall.
 * calResult2.cpp is used to calculate the evaluation of precision and recall.
 * The label of genome is as the first two keywords of nomenclature of gene feature.
 * The parameter -i and -l corresponding to the cluster of sequences and genomes.
 * The input cluster file is the CD-HIT format.
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

void calNMI(string application, string argument, string inputFile, string outputFile)
{	
	ofstream ofs(outputFile);
	fstream fs(inputFile);
	string line;
	unordered_map<string, int> mapClust;

	int curClustIndex = 0;
	int curStandardIndex = 0;
	vector<int> ourClust;
	vector<int> standardClust;

	while(getline(fs, line))
	{
		if(line.length() == 0) continue;
		if(application == "MeshClust2")
		{
			if(line[0] == '>')
			{
				curClustIndex++;
			}
			else{
				ourClust.push_back(curClustIndex);
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
				if(mapClust.find(key) == mapClust.end())//new standard cluster
				{
					mapClust.insert({key, curStandardIndex});
					standardClust.push_back(curStandardIndex);
					curStandardIndex++;
				}
				else
				{
					standardClust.push_back(mapClust[key]);
				}
			}
		}//end MeshClust2
		else //other application(RabbitTClust, gclust, Mothur)
		{
			if(line[0] != '\t')
			{
				curClustIndex++;
			}
			else{
				ourClust.push_back(curClustIndex);
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
					cerr << "err input application: " << application <<endl;
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
				if(mapClust.find(key) == mapClust.end())//new standard cluster
				{
					mapClust.insert({key, curStandardIndex});
					standardClust.push_back(curStandardIndex);
					curStandardIndex++;
				}
				else
				{
					standardClust.push_back(mapClust[key]);
				}
			}
		}
	}//end while

	cerr << "the size of ourClust is: " << ourClust.size() << endl;
	cerr << "the size of standardClust is: " << standardClust.size() << endl;

	for(int i = 0; i < ourClust.size(); i++)
		ofs << ourClust[i] << ' ';
	ofs << endl;
	
	for(int i = 0; i < standardClust.size(); i++)
		ofs << standardClust[i] << ' ';
	ofs << endl;

}


int main(int argc, char* argv[]){
	if(argc < 5){
		cerr << "run with: ./nmi RabbitTClust -l(-i) bacteria.out bacteria.nmi" << endl;
		cerr << "where -l means cluster by genomes, -i means cluster by sequences" << endl;
		return 1;
	}
	string application = argv[1];
	string argument = argv[2];
	string inputFile = argv[3];
	string outputFile = argv[4];

	calNMI(application, argument, inputFile, outputFile);

	return 0;
}
