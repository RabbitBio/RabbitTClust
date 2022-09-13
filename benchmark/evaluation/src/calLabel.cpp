/* Author: Xiaoming Xu
 * Data: 2022/6/9
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
#include "groundTruth.h"


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

void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions);

int updateLabel(vector< vector<LabNum> > &labNumArr, unordered_map<int, GlobalLabelInfo> &globalMap, int clustId, int &badLabel, vector<int> &resLabelArr);

void calLabelFile(string groundTruth, string clustFile, string labelFile);

void calLabelSequence(string groundTruth, string clustFile, string labelFile);

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
	string example = application + " bacteria.groundTruth -l bacteria.mst.clust bacteria.mst.label";
	args.push_back("groundTruth");
	args.push_back("sketchOption");
	args.push_back("clustFile");
	args.push_back("labelFile");
	descriptions.push_back("input file, groundTruth file, <accession, taxid, organismName> per line(first line is header)");
	descriptions.push_back("input option, sketch options, -l or -i, -l means sketchByFile, -i means sketchBySequence");
	descriptions.push_back("input file, cluster result file need to be labeled");
	descriptions.push_back("output file, label file according the groundTruth");

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
	string groundTruth = argv[1];
	string option = argv[2];
	string clustFile = argv[3];
	string labelFile = argv[4];

	bool sketchByFile;
	if(option == "-l")	sketchByFile = true;
	else if(option == "-i") sketchByFile = false;
	else{
		cerr << "error option: " << option << ", need -l or -i " << endl;
		return 1;
	}
	if(sketchByFile){
		calLabelFile(groundTruth, clustFile, labelFile);
	}
	else{
		calLabelSequence(groundTruth, clustFile, labelFile);
	}

  return 0;
}
void calLabelSequence(string groundTruth, string clustFile, string labelFile){
	//--------for groundTruth--------------
	unordered_map<string, int> seqName_taxid_map;
	unordered_map<int, string> taxid_organismName_map;
	getGroundTruthBySequence(groundTruth, seqName_taxid_map, taxid_organismName_map);

	//--------for cluster file--------------------------
	vector<int> ourClust;
	vector<int> standardClust;
	unordered_map<string, int> standardMap;
	unordered_map<int, int> curMap;
	vector<vector<LabNum>> labNumArr;
	vector<PosNum> posArr;
	int startPos = 0;
	string line;
	
	int numNotInGroundTruth = 0;
	ifstream ifs1(clustFile);
	if(!ifs1){
		cerr << "error open: " << clustFile << endl;
		exit(1);
	}
	while(getline(ifs1, line)){
		if(line[0] != '\t'){
			if(curMap.size() != 0){
				int clustSize = 0;
				vector<LabNum> curClustInfo;
				for(auto x : curMap){
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
			ss >> curId >> genomeId >> genomeSize >> genomeName;
			string key = genomeName;
			if(seqName_taxid_map.find(key) == seqName_taxid_map.end()){
				numNotInGroundTruth++;
				continue;
			}
			else{
				int curLabel = seqName_taxid_map[key];
				standardClust.push_back(curLabel);
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
			}
		}
	}
	if(curMap.size() != 0){
		int clustSize = 0;
		vector<LabNum> curClustInfo;
		for(auto x : curMap){
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

	//-------------for update labels------------------------------------
	unordered_map<int, GlobalLabelInfo> globalMap;
	int badLabel = -1;
	int clustNumber = labNumArr.size();
	vector<int> resLabelArr;
	resLabelArr.resize(clustNumber); 
	for(int i = 0; i < clustNumber; i++)
	{
		badLabel = updateLabel(labNumArr, globalMap, i, badLabel, resLabelArr); 
	}
	for(int i = 0; i < posArr.size(); i++)
	{
		int startPos = posArr[i].startPos;
		int clustSize = posArr[i].clustSize;
		for(int j = 0; j < clustSize; j++)
		{
			ourClust.push_back(resLabelArr[i]);
		}
	}
	cerr << "the number of which not in the groundTruth is: " << numNotInGroundTruth << endl;
	cerr << "the size of ourClust is: " << ourClust.size() << endl;
	cerr << "the size of standardClust is: " << standardClust.size() << endl;
	if(ourClust.size() != standardClust.size())
	{
		cerr << "the size of ourClust is not equal to the standardClust, exit()" << endl;
		return;
	}

	//--------------------for output labels-------------------------------------
	ofstream ofs(labelFile);
	for(int i = 0; i < ourClust.size(); i++)
		ofs << ourClust[i] << ' ';
	ofs << endl;
	for(int i = 0; i < standardClust.size(); i++)
		ofs << standardClust[i] << ' ';
	ofs << endl;
	ofs.close();

	ofstream ofs1(labelFile+".humanReadable");
	for(int i = 0; i < ourClust.size(); i++)
	{
		ofs1 << ourClust[i] << '\t' << standardClust[i] << endl;
	}
	ofs1.close();
}

void calLabelFile(string groundTruth, string clustFile, string labelFile){
	//--------for groundTruth--------------
	unordered_map<string, int> accession_taxid_map;
	unordered_map<int, string> taxid_organismName_map;
	getGroundTruthByFile(groundTruth, accession_taxid_map, taxid_organismName_map);

	//--------for cluster file--------------------------
	vector<int> ourClust;
	vector<int> standardClust;
	unordered_map<string, int> standardMap;
	unordered_map<int, int> curMap;
	vector<vector<LabNum>> labNumArr;
	vector<PosNum> posArr;
	int startPos = 0;
	
	int numNotInGroundTruth = 0;
	ifstream ifs1(clustFile);
	if(!ifs1){
		cerr << "error open: " << clustFile << endl;
		exit(1);
	}
	string line;
	while(getline(ifs1, line)){
		if(line[0] != '\t'){
			if(curMap.size() != 0){
				int clustSize = 0;
				vector<LabNum> curClustInfo;
				for(auto x : curMap){
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
			ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName;
			int startIndex = fileName.find_last_of('/');
			int endIndex = fileName.find_first_of('_', startIndex+5);
			if(endIndex == -1)	endIndex = fileName.find('.', startIndex+5);
			string key = fileName.substr(startIndex+1, endIndex-startIndex-1);
			if(accession_taxid_map.find(key) == accession_taxid_map.end()){
				numNotInGroundTruth++;
				continue;
			}
			else{
				int curLabel = accession_taxid_map[key];
				standardClust.push_back(curLabel);
				curMap.insert({curLabel, 0});
				curMap[curLabel]++;
			}
		}
	}
	if(curMap.size() != 0){
		int clustSize = 0;
		vector<LabNum> curClustInfo;
		for(auto x : curMap){
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

	//-------------for update labels------------------------------------
	unordered_map<int, GlobalLabelInfo> globalMap;
	int badLabel = -1;
	int clustNumber = labNumArr.size();
	vector<int> resLabelArr;
	resLabelArr.resize(clustNumber); 
	for(int i = 0; i < clustNumber; i++)
	{
		badLabel = updateLabel(labNumArr, globalMap, i, badLabel, resLabelArr); 
	}
	for(int i = 0; i < posArr.size(); i++)
	{
		int startPos = posArr[i].startPos;
		int clustSize = posArr[i].clustSize;
		for(int j = 0; j < clustSize; j++)
		{
			ourClust.push_back(resLabelArr[i]);
		}
	}
	cerr << "the number of which not in the groundTruth is: " << numNotInGroundTruth << endl;
	cerr << "the size of ourClust is: " << ourClust.size() << endl;
	cerr << "the size of standardClust is: " << standardClust.size() << endl;
	if(ourClust.size() != standardClust.size())
	{
		cerr << "the size of ourClust is not equal to the standardClust, exit()" << endl;
		return;
	}

	//--------------------for output labels-------------------------------------
	ofstream ofs(labelFile);
	for(int i = 0; i < ourClust.size(); i++)
		ofs << ourClust[i] << ' ';
	ofs << endl;
	for(int i = 0; i < standardClust.size(); i++)
		ofs << standardClust[i] << ' ';
	ofs << endl;
	ofs.close();

	//ofstream ofs1(labelFile+".humanReadable");
	//for(int i = 0; i < ourClust.size(); i++)
	//{
	//	ofs1 << ourClust[i] << '\t' << standardClust[i] << endl;
	//}
	//ofs1.close();

}

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
