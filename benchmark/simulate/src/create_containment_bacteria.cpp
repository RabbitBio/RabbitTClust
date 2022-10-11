/* Author: Xiaoming Xu
 * Data: 2022/5/31
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
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <zlib.h>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread);

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
	string pwd = "RabbitTClust/benchmark/simulate/src/create_containment_bacteria.cpp";
	string dependency = "kseq.h";
	string example = application + " simulateList num_of_clust num_genomes_each_clust simulatePath";
	args.push_back("simulateList");
	args.push_back("num_of_clust");
	args.push_back("num_genomes_each_clust");
	args.push_back("simulatePath");
	descriptions.push_back("input file, genome file list, one genome path per line");
	descriptions.push_back("input parameter, the number of clusters");
	descriptions.push_back("input parameter, the number of genomes in each cluster");
	descriptions.push_back("output path, the output file path");

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
	string inputList = argv[1];
	int numClusts = stoi(argv[2]);
	int numGenomePerClust = stoi(argv[3]);
	string outputPath = argv[4];
	string cmd0 = "mkdir -p " + outputPath;
	system(cmd0.c_str());

	//get randoms
	vector<double> randArr;
	for(int i = 0; i < numClusts* numGenomePerClust; i++)
	{
		int randNumber = rand() % 1000;
		randArr.push_back((double)randNumber/1000);
	}

	fstream fs0(inputList);
	unordered_set<string> seedSet;
	vector<string> seedArr;
	string line;
	
	while(getline(fs0, line))
	{
		int startIndex = line.find('/');
		int endIndex = line.find('.');
		string key = line.substr(startIndex + 1, endIndex-startIndex-1);
		if(seedSet.find(key) == seedSet.end())
		{
			seedArr.push_back(line);
			seedSet.insert(key);
		}
	}
	//cerr << "the size of seedArr is: " << seedArr.size() << endl;
	
	//string groundTruthFile = outputPath + "/groundTruth";
	//ofstream ofs(groundTruthFile);
	int seedArrSize = seedArr.size();
	int minNumber = std::min(numClusts, seedArrSize);
	for(int i = 0; i < minNumber; i++)
	{
		string cp_command = "cp " + seedArr[i] + " " + outputPath + '/';
		cerr << cp_command << endl;
		system(cp_command.c_str());
		gzFile fp1 = gzopen(seedArr[i].c_str(), "r");
		if(!fp1){
			cerr << "cannot open file: " << seedArr[i] << endl;
			continue;
		}
		kseq_t *ks1 = kseq_init(fp1);
		vector<string> bufArr;

		FILE *fpArr[numGenomePerClust];
		int startIndex = seedArr[i].find('/');
		string keyName = seedArr[i].substr(startIndex+1);
		int indexStr = keyName.find_last_of('.');
		keyName = keyName.substr(0, indexStr);

		for(int j = 0; j < numGenomePerClust; j++)
		{
			string outputName = outputPath + '/' + keyName + '.' + to_string(j) + ".fna";
			fpArr[j] = fopen(outputName.c_str(), "w");
			string writeBuffer("");
			bufArr.push_back(writeBuffer);
		}

		int index = 0;
		while(1)
		{
			int length = kseq_read(ks1);
			if(length < 0) break;
			string name = ks1->name.s;
			string comment = ks1->comment.s;
			string content = ks1->seq.s;
			string headLine = '>' + name + ' ' + comment + '\n';
			for(int j = 0; j < numGenomePerClust; j++)
			{
				bufArr[j] += headLine;

				int readLength = length * randArr[i*numGenomePerClust + j];
				//cerr << "the readLength is: " << readLength << endl;
				for(int k = 0; k < readLength; k +=80)
				{
					int actualLen = std::min(80, readLength - k);
					string tmpContent = content.substr(k, actualLen);
					bufArr[j] += tmpContent + '\n';
				}
				index++;
			}
		}
		for(int j = 0; j < numGenomePerClust; j++)
		{
			fwrite(bufArr[j].c_str(), sizeof(char), bufArr[j].length(), fpArr[j]);
			fclose(fpArr[j]);
		}

		gzclose(fp1);
		kseq_destroy(ks1);
	}
	//ofs.close();

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
