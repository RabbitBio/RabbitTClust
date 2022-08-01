/* Author: Xiaoming Xu
 * Data: 2022/7/22
 * 
 */
#include <iostream>
#include <stdlib.h>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <sys/sysinfo.h>
#include <omp.h>
#include <set>
#include <math.h>
#include <zlib.h>
#include "groundTruth.h"

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


bool cmpPurityNumber(PurityInfo p1, PurityInfo p2){
	return p1.totalNumber > p2.totalNumber;
}
bool cmpPurityPurity(PurityInfo p1, PurityInfo p2){
	if(p1.purity < p2.purity) return true;
	else if(p1.purity > p2.purity) return false;
	else return p1.totalNumber > p2.totalNumber;

}
bool cmpIdNum(IdNameNum id1, IdNameNum id2){
	return id1.number > id2.number;
}

void calPurityFile(string groundTruth, string clustFile, string outputFile);
void calPuritySequence(string groundTruth, string clustFile, string outputFile);


void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions);

int main(int argc , char *argv[]){
	string application = argv[0];
	vector<string> args, descriptions;
	args.push_back(application);
	descriptions.push_back("the application name");

	//========= parameters need changing ========
	//The example is with parameters of specific numbers.
	//args is the tutorial names.
	string pwd = "RabbitTClust/benchmark/evaluation/src";
	string dependency = "None";
	string example = application + " -l bacteria.groundTruth bacteria.clust partPurity.out";
	args.push_back("options(-l, -i)");
	args.push_back("groundTruth");
	args.push_back("clustFile");
	args.push_back("outputPurityFile");
	descriptions.push_back("input option, sketch option for clust, -l or -i");
	descriptions.push_back("input file, the groundTruth file, <assembly_accession, species_taxid, genomeName> per line");
	descriptions.push_back("input file, cluster result file from RabbitTClust");
	descriptions.push_back("output file, output purity info file");

	//-------- no changing -----------
	assert(args.size() == descriptions.size());
  if(argc != args.size()) {
		printInfo(pwd, dependency, example, args, descriptions);
    return 1;
  }
	else if(argc == 2 && (argv[1] == "-h" || argv[1] == "--help"))
	{
		printInfo(pwd, dependency, example, args, descriptions);
		return 1;
	}

	//======== specific implement ========
	string argument = argv[1];
	string groundTruth = argv[2];
	string clustFile = argv[3];
	string outputFile = argv[4];

	bool sketchByFile;
	if(argument == "-l")	sketchByFile = true;
	else if(argument == "-i") sketchByFile = false;
	else{
		cerr << "error option: " << argument << ", need -l or -i " << endl;
		return 1;
	}
	if(sketchByFile){
		calPurityFile(groundTruth, clustFile, outputFile);
	}
	else{
		calPuritySequence(groundTruth, clustFile, outputFile);
		cerr << "have not implemented the calPurity by sequence" << endl;
	}

  return 0;
}
void calPuritySequence(string groundTruth, string clustFile, string outputFile){
	//--------for groundTruth--------------
	unordered_map<string, int> seqName_taxid_map;
	unordered_map<int, string> taxid_organismName_map;
	getGroundTruthBySequence(groundTruth, seqName_taxid_map, taxid_organismName_map);

	//--------for cluster file--------------------------
	unordered_map<int, int> curMap;
	vector<int> totalArr;
	vector<SpeciesIdNumInfo> dominantArr;
	int numNotIngroundTruth = 0;
	string line;

	ifstream ifs1(clustFile);
	while(getline(ifs1, line)){
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
			ss >> curId >> genomeId >> genomeSize >> genomeName >> type0 >> type1 >> type2;
			string key = genomeName;
			if(seqName_taxid_map.find(key) == seqName_taxid_map.end()){
				numNotIngroundTruth++;
				continue;
			}
			else{
				int curLabel = seqName_taxid_map[key];
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
			}
		}
	}//end while getline
	ifs1.close();
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

	vector<PurityInfo> partPurity;
	assert(totalArr.size() == dominantArr.size());
	int totalGenomeNumber = 0;
	int totalDominantNumber = 0;
	int totalCoverageNumber = 0;
	for(int i = 0; i < totalArr.size(); i++)
	{
		if(totalArr[i] > 1){
			totalCoverageNumber += totalArr[i];
		}
		totalGenomeNumber += totalArr[i];
		totalDominantNumber += dominantArr[i].number;
		double curPurity = (double)dominantArr[i].number / (double)totalArr[i];
		int curDominantSpeciesId = dominantArr[i].species_id;
		string curDominantOrganism = taxid_organismName_map[curDominantSpeciesId];
		PurityInfo curPur(totalArr[i], dominantArr[i].number, curPurity, curDominantSpeciesId, curDominantOrganism);
		partPurity.push_back(curPur);
	}
	std::sort(partPurity.begin(), partPurity.end(), cmpPurityNumber);

	assert(totalArr.size() == partPurity.size());
	ofstream ofs(outputFile);
	FILE* fp = fopen(outputFile.c_str(), "w");
	fprintf(fp, "Purity\ttotalNumber\tdominateNumber\tdominateSpeciesId\tdominateOriganism\n");
	for(int i = 0; i < totalArr.size(); i++)
	{
		fprintf(fp, "%8lf\t%8d\t%8d\t\t%8d\t%s\n", partPurity[i].purity, partPurity[i].totalNumber, partPurity[i].dominateNumber, partPurity[i].dominateSpeciesId, partPurity[i].dominateOrganism.c_str());
	}
	fclose(fp);

	double totalPurity = (double)totalDominantNumber / (double)totalGenomeNumber;
	double totalCoverageRate = (double)totalCoverageNumber / (double)totalGenomeNumber;
	cerr << "the coverage is: " << totalCoverageRate << endl;
	cerr << "the final purity is: " << totalPurity << endl;
	cerr << "the total genome number of " << clustFile << " is: " << totalGenomeNumber << endl;
	cerr << "the total dominant genome number of " << clustFile << " is: " << totalDominantNumber << endl;

}

void calPurityFile(string groundTruth, string clustFile, string outputFile){
	//--------for groundTruth--------------
	unordered_map<string, int> accession_taxid_map;
	unordered_map<int, string> taxid_organismName_map;
	getGroundTruthByFile(groundTruth, accession_taxid_map, taxid_organismName_map);

	//--------for cluster file--------------------------
	unordered_map<int, int> curMap;
	vector<int> totalArr;
	vector<SpeciesIdNumInfo> dominantArr;
	int numNotIngroundTruth = 0;
	string line;

	ifstream ifs1(clustFile);
	while(getline(ifs1, line)){
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
			ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName >> type0 >> type1 >> type2;
			int startIndex = fileName.find_last_of('/');
			int endIndex = fileName.find('_', startIndex + 5);
			if(fileName.find('_', startIndex+5) == -1)
				endIndex = fileName.find('.', startIndex+5);
			string key = fileName.substr(startIndex+1, endIndex -startIndex -1);
			if(accession_taxid_map.find(key) == accession_taxid_map.end()){
				numNotIngroundTruth++;
				continue;
			}
			else{
				int curLabel = accession_taxid_map[key];
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
			}
		}
	}//end while getline
	ifs1.close();
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

	vector<PurityInfo> partPurity;
	assert(totalArr.size() == dominantArr.size());
	int totalGenomeNumber = 0;
	int totalDominantNumber = 0;
	int totalCoverageNumber = 0;
	for(int i = 0; i < totalArr.size(); i++)
	{
		if(totalArr[i] > 1){
			totalCoverageNumber += totalArr[i];
		}
		totalGenomeNumber += totalArr[i];
		totalDominantNumber += dominantArr[i].number;
		double curPurity = (double)dominantArr[i].number / (double)totalArr[i];
		int curDominantSpeciesId = dominantArr[i].species_id;
		string curDominantOrganism = taxid_organismName_map[curDominantSpeciesId];
		PurityInfo curPur(totalArr[i], dominantArr[i].number, curPurity, curDominantSpeciesId, curDominantOrganism);
		partPurity.push_back(curPur);
	}
	std::sort(partPurity.begin(), partPurity.end(), cmpPurityNumber);
	//std::sort(partPurity.begin(), partPurity.end(), cmpPurityPurity);

	double minPurity = 1.0;
	assert(totalArr.size() == partPurity.size());
	ofstream ofs(outputFile);
	FILE* fp = fopen(outputFile.c_str(), "w");
	fprintf(fp, "Purity\ttotalNumber\tdominateNumber\tdominateSpeciesId\tdominateOriganism\n");
	for(int i = 0; i < totalArr.size(); i++)
	{
		minPurity = std::min(minPurity, partPurity[i].purity);
		fprintf(fp, "%8lf\t%8d\t%8d\t\t%8d\t%s\n", partPurity[i].purity, partPurity[i].totalNumber, partPurity[i].dominateNumber, partPurity[i].dominateSpeciesId, partPurity[i].dominateOrganism.c_str());
	}
	fclose(fp);

	double totalPurity = (double)totalDominantNumber / (double)totalGenomeNumber;
	double totalCoverageRate = (double)totalCoverageNumber / (double)totalGenomeNumber;
	cerr << "the coverage is: " << totalCoverageRate << endl;
	cerr << "the final purity is: " << totalPurity << endl;
	cerr << "the total genome number of " << clustFile << " is: " << totalGenomeNumber << endl;
	cerr << "the total dominant genome number of " << clustFile << " is: " << totalDominantNumber << endl;
	cerr << "the minimum purity of cluster is: " << minPurity << endl;

}





void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions){
	assert(args.size() == descriptions.size());
	cerr << endl;
	cerr << "example: " << example << endl;
	cerr << endl;
	cerr << "source file path: " << pwd << endl;
	cerr << endl;
	cerr << "dependency: " << dependency << endl;
	cerr << endl;
	cerr << "run as: ";
	for(int i = 0; i < args.size(); i++){
		cerr << args[i] << ' ';
	}
	cerr << endl;
	for(int i = 0; i < args.size(); i++){
		fprintf(stderr, "\tThe %d paramter(%s) is %s\n", i, args[i].c_str(), descriptions[i].c_str());
	}
}
