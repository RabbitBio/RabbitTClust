/* Input file is the inputFile list.
 * Output file(nonRedundantList) contains the fileName for each line.
 * All the fileNames in the output file(nonRedundantList) have very low similarity with each other.
 * The nonRedundantList is to reduce the distraction of similar genomes with generating simulated database(very different genome size).
 *
 */

#include <iostream>
#include "kseq.h"
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <sstream>
#include <unordered_map>


KSEQ_INIT(gzFile, gzread);
using namespace std;

//struct Info{
//	string name;
//	string comment;
//	int length;
//};


int main(int argc, char * argv[]){
	
	if(argc < 2) return 1;
	string fileName = argv[1];

	fstream fs(fileName);
	if(!fs) return 1;

	vector<string> fileList;
	string tmpName;
	while(getline(fs, tmpName))
		fileList.push_back(tmpName);


	//vector< vector<Info> > clusters;
	string outputFile = "nonRedundantList";
	FILE *fp = fopen(outputFile.c_str(), "w");

	unordered_map<string, int> mapFile;
	cerr << "the size of fileList is: " << fileList.size() << endl;
	int index = 0;

	for(int i = 0; i < fileList.size(); i++)
	{
		gzFile fp1;
		kseq_t * ks1;
		fp1 = gzopen(fileList[i].c_str(), "r");
		if(!fp1){
			cerr << "can not open: " << fileList[i] << endl;
			return 1;
		}
		ks1 = kseq_init(fp1);

		int length = kseq_read(ks1);
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
		auto got = mapFile.find(key);
		if(got == mapFile.end())
		{
			fprintf(fp, "%s\n", fileList[i].c_str());
			index++;
		}
		mapFile.insert({key, 0});
		mapFile[key]++;

		gzclose(fp1);
		kseq_destroy(ks1);

	}
	cerr << "after redundant removement, the size is: " << index << endl;

	return 0;
}

