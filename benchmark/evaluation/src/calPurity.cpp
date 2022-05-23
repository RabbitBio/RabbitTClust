/* Author: Xiaoming Xu
 * Email: xiaoming.xu@mail.sdu.edu.cn
 * Data: 2022/5/23
 *
 * calPurity.cpp is used for get the purity of the clustering.
 * The groundTruth labels of genomes are as species_taxid.
 * The parameter -i and -l corresponding to the cluster of genomes served as sequences and files.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <unordered_set>
#include <algorithm>

using namespace std;

struct PurityInfo{
	int totalNumber;
	int dominateNumber;
	double purity;
	int dominateSpeciesId;
	string dominateOrganism;
	PurityInfo(int a, int b, double c, int d, string e):totalNumber(a), dominateNumber(b), purity(c), dominateSpeciesId(d), dominateOrganism(e){}
};

struct SpeciesIdNumInfo{
	int species_id;
	int number;
	SpeciesIdNumInfo(int a, int b):species_id(a), number(b){}
};

struct IdNameNum{
	int id;
	string name;
	int number;
	IdNameNum(int a, string b, int c):id(a), name(b), number(c){}
};


bool cmpPurity(PurityInfo p1, PurityInfo p2){
	return p1.totalNumber > p2.totalNumber;
}
bool cmpIdNum(IdNameNum id1, IdNameNum id2){
	return id1.number > id2.number;
}

inline void printInfo(string args)
{
	cerr << "run with: " << args << " -l(-i) groundTruth bacteria.clust partPurity.out" << endl;
	cerr << "The second argument, -l means genomes served as files, -i means genomes served as sequences" << endl;
	cerr << "The third argument (groundTruth) is the ground truth from assembly_bacteria.txt of the <assembly_accession genomeName species_taxid> " << endl;
	cerr << "The fourth argument (bacteria.clust) is the cluster result from RabbitTClust, MeshClust2, gclust or Mothur " << endl;
	cerr << "The fifth argument (partPurity.out) is the output file path for detailed purity of each cluster" << endl;
	
}
double calPurity(string args, string argument, string groundTruth, string inputFile, string outputFile){
	if(argument != "-l" && argument != "-i"){
		printInfo(args);
		return 0.0;
	}
	fstream fs0(groundTruth);
	string line;

	unordered_map<string, int> groundTruthMapFile;
	unordered_map<string, int> groundTruthMapSeq;
	unordered_map<int, string> groundTruthMapIdName;
	unordered_map<int, int> groundTruthMapIdNumber;

	while(getline(fs0, line))
	{
		string assembly_accession, genomeName, tmpName, organism_name("");
		int species_taxid;
		stringstream ss;
		ss << line;
		ss >> assembly_accession >> genomeName >> species_taxid;
		while(ss >> tmpName){
			organism_name += tmpName + ' ';
		}
		organism_name = organism_name.substr(0, organism_name.length()-1);

		groundTruthMapFile.insert({assembly_accession, species_taxid});
		groundTruthMapSeq.insert({genomeName, species_taxid});
		groundTruthMapIdName.insert({species_taxid, organism_name});
		groundTruthMapIdNumber.insert({species_taxid, 0});
		groundTruthMapIdNumber[species_taxid]++;

	}
	fs0.close();
	cerr << "the groundTruthMapSpeciesIdName size is: " << groundTruthMapIdName.size() << endl;


	int outIndex = outputFile.find_first_of('.');
	string groundTruthOutput = "default_groundTruthIdNumber";
	if(outIndex != -1)	
		groundTruthOutput = outputFile.substr(0, outIndex) + "_groundTruthIdNumber";
	ofstream ofsT(groundTruthOutput);
	vector<IdNameNum> idNameNumArr;
	for(auto x : groundTruthMapIdNumber){
		int id = x.first;
		int number = x.second;
		string name = groundTruthMapIdName[id];
		IdNameNum curIdNameNum(id, name, number);
		idNameNumArr.push_back(curIdNameNum);
	}
	std::sort(idNameNumArr.begin(), idNameNumArr.end(), cmpIdNum);
	ofsT << "Species_taxid_number\tSpecies_taxid\tOrganismName\n";
	for(int i = 0; i <idNameNumArr.size(); i++){
		ofsT << idNameNumArr[i].number << '\t' << idNameNumArr[i].id << '\t' << idNameNumArr[i].name << endl;
	}

	ofsT.close();

	
	
	int numNotIngroundTruth = 0;
	unordered_map<int, int> curMap;
	vector<int> totalArr;
	vector<SpeciesIdNumInfo> dominantArr;

	fstream fs(inputFile);
	while(getline(fs, line))
	{
		if(line.length() == 0) continue;
		if(line[0] != '\t' && curMap.size() != 0){
			int curTotalNum = 0;
			int maxNum = 0;
			int maxSpeciesId = 0;
			for(auto x : curMap){
				curTotalNum += x.second;
				if(maxNum < x.second){
					maxNum = x.second;
					maxSpeciesId = x.first;
				}
			}
			totalArr.push_back(curTotalNum);
			SpeciesIdNumInfo id_num(maxSpeciesId, maxNum);
			dominantArr.push_back(id_num);
			unordered_map<int, int>().swap(curMap);
		}
		else{
			stringstream ss;
			ss << line;
			int curId, genomeId;
			string genomeSize, fileName, genomeName;
			string type0, type1, type2;
			if(argument == "-l")
			{
				ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName >> type0 >> type1 >> type2;
				int startIndex = fileName.find_last_of('/');
				int endIndex = fileName.find('_', startIndex + 5);
				if(fileName.find('_', startIndex+5) == -1)
					endIndex = fileName.find('.', startIndex+5);
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
					curMap.insert({curLabel, 0});
					curMap[curLabel]++;
				}
			}
		}

	}
	fs.close();
	if(curMap.size() != 0){
		int curTotalNum = 0;
		int maxNum = 0;
		int maxSpeciesId = 0;
		for(auto x : curMap){
			curTotalNum += x.second;
			if(maxNum < x.second){
				maxNum = x.second;
				maxSpeciesId = x.first;
			}
		}
		totalArr.push_back(curTotalNum);
		SpeciesIdNumInfo id_num(maxSpeciesId, maxNum);
		dominantArr.push_back(id_num);
		unordered_map<int, int>().swap(curMap);
	}

	//vector<double> partPurity;
	vector<PurityInfo> partPurity;
	assert(totalArr.size() == dominantArr.size());
	int totalGenomeNumber = 0;
	int totalDominantNumber = 0;
	for(int i = 0; i < totalArr.size(); i++)
	{
		totalGenomeNumber += totalArr[i];
		totalDominantNumber += dominantArr[i].number;
		double curPurity = (double)dominantArr[i].number / (double)totalArr[i];
		int curDominantSpeciesId = dominantArr[i].species_id;
		string curDominantOrganism = groundTruthMapIdName[curDominantSpeciesId];
		PurityInfo curPur(totalArr[i], dominantArr[i].number, curPurity, curDominantSpeciesId, curDominantOrganism);
		partPurity.push_back(curPur);
	}
	std::sort(partPurity.begin(), partPurity.end(), cmpPurity);

	assert(totalArr.size() == partPurity.size());
	ofstream ofs(outputFile);
	ofs << "Purity\ttotalNumber\tdominateNumber\tdominateSpeciesId\tdominateOriganism\n";
	for(int i = 0; i < totalArr.size(); i++)
	{
		ofs << partPurity[i].purity << '\t' << partPurity[i].totalNumber << '\t' << partPurity[i].dominateNumber << '\t' << partPurity[i].dominateSpeciesId << '\t' << partPurity[i].dominateOrganism << endl;
	}
	ofs.close();

	double totalPurity = (double)totalDominantNumber / (double)totalGenomeNumber;
	cerr << "the total genome number of " << inputFile << " is: " << totalGenomeNumber << endl;
	cerr << "the total dominant genome number of " << inputFile << " is: " << totalDominantNumber << endl;
	return totalPurity;

}

int main(int argc, char* argv[]){
	if(argc != 5){
		printInfo(argv[0]);
		return 1;
	}
	string args = argv[0];
	string argument = argv[1];
	string groundTruth = argv[2];
	string inputFile = argv[3];
	string outputFile = argv[4];

	double purity = calPurity(args, argument, groundTruth, inputFile, outputFile);
	cout << "the final purity of " << inputFile << " is: " << purity << endl;
	cout << "===============================================================" << endl;

	return 0;

}





