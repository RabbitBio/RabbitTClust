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
#include <utility>

using namespace std;

//bool reserveAllRepresentativeGenome = true;
bool reserveAllRepresentativeGenome = false;

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

bool cmpVecStr(pair<int, vector<string>> pvs1, pair<int, vector<string>> pvs2){
	return pvs1.second.size() > pvs2.second.size();
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
	string pwd = "RabbitTClust/benchmark/evaluation/src/calPurity.cpp";
	string dependency = "None";
	string example = application + " -l bacteria.groundTruth bacteria.clust bacteria.purity";
	args.push_back("options(-l, -i)");
	args.push_back("groundTruth");
	args.push_back("clustFile");
	args.push_back("bacteria.purity");
	descriptions.push_back("input option, sketch option for clust, -l or -i");
	descriptions.push_back("input file, the groundTruth file, <assembly_accession, species_taxid, genomeName> per line");
	descriptions.push_back("input file, cluster result file from RabbitTClust");
	descriptions.push_back("output file, output purity info file, including total result file(bacteria.purity) and accession file(<accession, taxid> per line)");

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
	unordered_map<int, vector<string>> curIdAccessionsMap;
	vector<vector<pair<string, int>>> badAccessionTaxidArr;
	vector<pair<string, int>> representAccessionTaxidArr;
	vector<vector<vector<pair<string, int>>>> allAccessionTaxidArr;
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
		
			vector<pair<int, vector<string>>> curAccessionsArr;
			for(auto x : curIdAccessionsMap){
				curAccessionsArr.push_back(std::make_pair(x.first, x.second));
			}
			std::sort(curAccessionsArr.begin(), curAccessionsArr.end(), cmpVecStr);

			vector<vector<pair<string, int>>> curClustArr;
			for(int i = 0; i < curAccessionsArr.size(); i++){
				int curId = curAccessionsArr[i].first;
				vector<pair<string, int>> curClustLabelArr;
				for(auto x : curAccessionsArr[i].second){
					curClustLabelArr.push_back(std::make_pair(x, curId));
				}
				curClustArr.push_back(curClustLabelArr);
				vector<pair<string, int>>().swap(curClustLabelArr);
			}
			allAccessionTaxidArr.push_back(curClustArr);
			vector<vector<pair<string, int>>>().swap(curClustArr);
			
			vector<pair<int, vector<string>>>().swap(curAccessionsArr);
			unordered_map<int, vector<string>>().swap(curIdAccessionsMap);
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
				curIdAccessionsMap.insert({curLabel, vector<string>()});
				curIdAccessionsMap[curLabel].push_back(key);
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

		vector<pair<int, vector<string>>> curAccessionsArr;
		for(auto x : curIdAccessionsMap){
			curAccessionsArr.push_back(std::make_pair(x.first, x.second));
		}
		std::sort(curAccessionsArr.begin(), curAccessionsArr.end(), cmpVecStr);

		vector<vector<pair<string, int>>> curClustArr;
		for(int i = 0; i < curAccessionsArr.size(); i++){
			int curId = curAccessionsArr[i].first;
			vector<pair<string, int>> curClustLabelArr;
			for(auto x : curAccessionsArr[i].second){
				curClustLabelArr.push_back(std::make_pair(x, curId));
			}
			curClustArr.push_back(curClustLabelArr);
			vector<pair<string, int>>().swap(curClustLabelArr);
		}
		allAccessionTaxidArr.push_back(curClustArr);
		vector<vector<pair<string, int>>>().swap(curClustArr);
		
		vector<pair<int, vector<string>>>().swap(curAccessionsArr);
		unordered_map<int, vector<string>>().swap(curIdAccessionsMap);
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
	FILE* fp = fopen(outputFile.c_str(), "w");
	fprintf(fp, "Purity\ttotalNumber\tdominateNumber\tdominateSpeciesId\tdominateOriganism\n");
	for(int i = 0; i < totalArr.size(); i++)
	{
		minPurity = std::min(minPurity, partPurity[i].purity);
		fprintf(fp, "%8lf\t%8d\t%8d\t\t%8d\t%s\n", partPurity[i].purity, partPurity[i].totalNumber, partPurity[i].dominateNumber, partPurity[i].dominateSpeciesId, partPurity[i].dominateOrganism.c_str());
	}
	fclose(fp);

	string outputFile2 = outputFile + ".accession.unpurity";
	ofstream ofs0(outputFile2);
	if(!reserveAllRepresentativeGenome){
		for(auto x : allAccessionTaxidArr){
			if(x.size() > 1){
				string repAccId = x[0][0].first + '\t' + to_string(x[0][0].second);
				ofs0 << repAccId << endl;
				for(int i = 1; i < x.size(); i++){
					for(auto y : x[i])
						ofs0 << '\t' << y.first << '\t' << y.second << endl;
				}
				ofs0 << endl;
			}
		}
	}
	else{
		for(auto x : allAccessionTaxidArr){
			if(x.size() > 1){
				for(auto y : x[0])
					ofs0 << y.first << '\t' << y.second << endl;
				for(int i = 1; i < x.size(); i++){
					for(auto y : x[i])
						ofs0 << '\t' << y.first << '\t' << y.second << endl;
				}
				ofs0 << endl;
			}
		}
	}
	ofs0.close();

	string outputFile3 = outputFile + ".accession.purity";
	ofstream ofs1(outputFile3);
	if(!reserveAllRepresentativeGenome){
		for(auto x : allAccessionTaxidArr){
			if(x.size() == 1){
				string repAccId = x[0][0].first + '\t' + to_string(x[0][0].second);
				ofs1 << repAccId << endl;
			}
		}
	}
	else{
		for(auto x : allAccessionTaxidArr){
			if(x.size() == 1){
				for(auto y : x[0])
					ofs1 << y.first << '\t' << y.second << endl;
				ofs1 << endl;
			}
		}
	}
	ofs0.close();


	//if(!reserveAllRepresentativeGenome){
	//assert(badAccessionTaxidArr.size() == representAccessionTaxidArr.size());
	//	for(int i = 0; i < representAccessionTaxidArr.size(); i++){
	//		string repAccId = representAccessionTaxidArr[i].first + '\t' + to_string(representAccessionTaxidArr[i].second);
	//		if(badAccessionTaxidArr[i].size() != 0){
	//			ofs << repAccId << endl;
	//			for(auto y : badAccessionTaxidArr[i]){
	//				ofs << '\t' << y.first << '\t' << y.second << endl;
	//			}
	//			ofs << endl;
	//		}
	//	}
	//}
	//else{

	//}

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
