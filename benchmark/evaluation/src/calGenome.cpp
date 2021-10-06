/* There are some genome files in the bacteriaList from refseq having the same nomenclature type0.
 * claGenome.cpp is used to calculate number of each nomenclature type in the bacteriaList.
 * The sum of all the number of type is the number of genome in bacteriaList.
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
	string outputFile = "genomeType2.out";
	fstream fs(inputFile);
	string line;
	vector<string> fileList;

	while(getline(fs, line))
	{
		fileList.push_back(line);
	}
	cout << "the size of fileList: " << fileList.size() << endl;

	//vector<Type> genomeType[fileList.size()];
	unordered_map<string, int> genomeType;
	FILE *fp = fopen(outputFile.c_str(), "w");
	int tmpIndex = 0;
	int subIndex = fileList.size() / 100;
	
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

		//only read the first sequence, since we have checked the species in one genome file is the same.
		int length = kseq_read(ks1);
		if(length < 0) continue;
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
		genomeType.insert({key, 0});
		genomeType[key]++;
		
		gzclose(fp1);
		kseq_destroy(ks1);

	}
	int totalNumber = 0;
	for(auto x : genomeType)
	{
		fprintf(fp, "%s\t%d\n", x.first.c_str(), x.second);
		totalNumber += x.second;
	}
	fprintf(fp, "\n");
	unordered_map<string, int>().swap(genomeType);
	cerr << totalNumber << endl;
	cerr << "finished" << endl;

	return 0;
}
