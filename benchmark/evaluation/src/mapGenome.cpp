/* There maybe not only one sequence in the genome file. 
 * mapGenome.cpp is used to check whether multi sequences within a genome file have the same nomenclature type in the comment description or not.
 * The result have show that all the sequences within the same genome file have the same nomenclature type.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include "kseq.h"
#include <zlib.h>
#include <sstream>
#include <unordered_map>
#include <vector>


KSEQ_INIT(gzFile, gzread);
using namespace std;
struct Type{
	string type;
	int number;
};


int main(int argc, char *argv[]){
	if(argc < 2) return 1;	
	string inputFile = argv[1];
	string outputFile = "mapType.out";
	fstream fs(inputFile);
	string line;
	vector<string> fileList;

	while(getline(fs, line))
	{
		fileList.push_back(line);
	}
	cout << "the size of fileList: " << fileList.size() << endl;

	vector<Type> genomeType[fileList.size()];
	FILE *fp = fopen(outputFile.c_str(), "w");
	int tmpIndex = 0;
	int subIndex = fileList.size() / 100;
	
	//#pragma omp parallel for num_threads(48)
	for(int i = 0; i < fileList.size(); i++)
	{
		if(i > tmpIndex * subIndex)
		{
			tmpIndex++;
			cerr << tmpIndex << endl;
		}
		gzFile fp1 = gzopen(fileList[i].c_str(), "r");
		kseq_t * ks1;
		ks1 = kseq_init(fp1);
		//cerr << "read the file " << fileList[i] << endl;
		unordered_map<string, int> genomeMap;
		while(1)
		{
			int length = kseq_read(ks1);
			if(length < 0) break;
			string name = ks1->name.s;
			string comment = ks1->comment.s;
			stringstream ss;	
			ss << comment;
			string type0, type1, type2;
			ss >> type0 >> type1 >> type2;
			if(type0.substr(0, 10) == "UNVERIFIED")
			{
				type0 = type1;
				type1 = type2;
			}
			if(type0.back() == ',') type0.pop_back();
			if(type1.back() == ',') type1.pop_back();

			string key = type0 + '\t' + type1;
			genomeMap.insert({key, 0});
			genomeMap[key]++;
		}
		gzclose(fp1);
		kseq_destroy(ks1);
		if(genomeMap.size() != 1)
		{
			cerr << "there are not only one class in the file: " << fileList[i] << endl;
			for(auto x : genomeMap)
			{
				cerr << "\t" << x.first << "\t" << x.second << endl;
			}
		}

		for(auto x : genomeMap)
		{
			//cout << x.first << "\t" << x.second << endl;
			fprintf(fp, "%s\t%d\n", x.first.c_str(), x.second);
			genomeType[i].push_back({x.first, x.second});
		}
		fprintf(fp, "\n");
		unordered_map<string, int>().swap(genomeMap);
	}
	cerr << "finished" << endl;

	return 0;
}
