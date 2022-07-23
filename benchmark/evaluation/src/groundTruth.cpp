#include "groundTruth.h"
#include <sstream>

void getGroundTruthBySequence(string groundTruth, unordered_map<string, int>& seqName_taxid_map, unordered_map<int, string>& taxid_organismName_map){

	ifstream ifs0(groundTruth);
	if(!ifs0){
		cerr << "error open: " << groundTruth << endl;
		exit(1);
	}
	string line;
	getline(ifs0, line);//for the header line
	while(getline(ifs0, line)){
		stringstream ss;
		string seqName, organismName(""), tmpStr;
		int taxid;
		ss << line;
		ss >> seqName >> taxid;
		while(ss >> tmpStr){
			organismName += tmpStr + ' ';
		}
		organismName.substr(0, organismName.length()-1);
		seqName_taxid_map.insert({seqName, taxid});
		taxid_organismName_map.insert({taxid, organismName});
	}
	ifs0.close();
}

void getGroundTruthByFile(string groundTruth, unordered_map<string, int>& accession_taxid_map, unordered_map<int, string>& taxid_organismName_map){

	ifstream ifs0(groundTruth);
	if(!ifs0){
		cerr << "error open: " << groundTruth << endl;
		exit(1);
	}
	string line;
	getline(ifs0, line);//for the header line
	while(getline(ifs0, line)){
		stringstream ss;
		string accession, organismName(""), tmpStr;
		int taxid;
		ss << line;
		ss >> accession >> taxid;
		while(ss >> tmpStr){
			organismName += tmpStr + ' ';
		}
		organismName.substr(0, organismName.length()-1);
		accession_taxid_map.insert({accession, taxid});
		taxid_organismName_map.insert({taxid, organismName});
	}
	ifs0.close();
}






