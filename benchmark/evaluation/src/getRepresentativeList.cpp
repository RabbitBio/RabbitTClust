/* Author: Xiaoming Xu
 * Data: 2022/7/16
 *
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


using namespace std;

void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions);

int main(int argc , char *argv[]){
	string application = argv[0];
	vector<string> args, descriptions;
	args.push_back(application);
	descriptions.push_back("the application name");

	//========= parameters need changing ========
	//The example is with parameters of specific numbers.
	//args is the tutorial names.
	string pwd = "RabbitTClust/benchmark/evaluation/src/getRepresentativeList.cpp";
	string dependency = "None";
	string example = application + " -l bacteria.greedy.clust bacteria_representative.list";
	args.push_back("-i/-l");
	args.push_back("clustFile");
	args.push_back("representative_list");
	descriptions.push_back("input parameter, sketch parameter for the cluster file, -i means sketchBySequence, -l means sketchByFile");
	descriptions.push_back("input file, the cluster result from RabbitTClust");
	descriptions.push_back("output file, the representative list of genomes file or sequences");

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
	string option = argv[1];
	if(option != "-l" && option != "-i"){
		cerr << "error option: " << option << ", need -l or -i option" << endl;
		return 1;
	}

	string clustFile = argv[2];
	string outputFile = argv[3];

	ifstream ifs(clustFile);
	string line;
	ofstream ofs(outputFile);
	bool isClust = false;
	while(getline(ifs, line)){
		if(line[0] != '\t'){
			isClust = true;
		}
		else if(isClust){
			isClust = false;
			stringstream ss;
			int curId, globalId;
			string length, fileName, seqName, comment, tmpComment;
			ss << line;
			if(option == "-l"){
				ss >> curId >> globalId >> length >> fileName >> seqName;
				ofs << fileName << endl;
			}
			else if(option == "-i"){
				ss >> curId >> globalId >> length >> seqName;
				ofs << seqName << endl;
			}
		}
		
	}
	ifs.close();

  return 0;
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
		fprintf(stderr, "\tThe %d parameter(%s) is %s\n", i, args[i].c_str(), descriptions[i].c_str());
	}
}
