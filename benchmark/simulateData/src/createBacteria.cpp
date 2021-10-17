/* Input file is the inputFile list.
 * Output file(subBacteria.fna) is the sub-set of bacteria database since the gclust cannot finish the whole bacteria genomes cluster in pratical time.
 * Each sequence in the subBacteria.fna represent one genome in bacteria database.
 */

#include <iostream>
#include "kseq.h"
#include <zlib.h>
#include <string>
#include <fstream>
#include <vector>


KSEQ_INIT(gzFile, gzread);

using namespace std;


int main(int argc, char * argv[]){
	if(argc < 2)
	{
		cerr << "need input fileList, exit!" << endl;
		return 1;
	}
	string inputFile = argv[1];
	fstream fs(inputFile);
	if(!fs)
	{	
		cerr << "cannot open: " << inputFile << ", exit!" << endl;
		return 1;
	}
	string line;
	vector<string> fileList;
	while(getline(fs, line))
	{
		fileList.push_back(line);
	}

	string outputFile = "subBacteria.fna";
	ofstream ofile;
	ofile.open(outputFile);
	uint64_t totalLength = 0;
	uint64_t lengthThreshold = (uint64_t)19848 * (1 << 20);
	cerr << "length threshold is: " << lengthThreshold << endl;
	int index = 0;

	for(int i = 0; i < fileList.size(); i++)
	{
		gzFile fp1 = gzopen(fileList[i].c_str(), "r");
		if(!fp1)
		{
			cerr << "cannot open: " << fileList[i] << ", exit! " << endl;
			return 1;
		}
		kseq_t * ks1 = kseq_init(fp1);
		int length = kseq_read(ks1);
		if(length < 0) break;
		string name = ks1->name.s;
		string comment = ks1->comment.s;
		string content = ks1->seq.s;
		string headLine = '>' + name + ' ' + comment;
		ofile << headLine << endl;
		for(int j = 0; j < length; j+=80)
		{
			string tmpContent = content.substr(j, 80);
			ofile << tmpContent << endl;
		}
		totalLength += length;
		index++;
		//if(totalLength >= 19848 * 1<<20 && index >= 112111) break;
		if(totalLength >= lengthThreshold) break;
		gzclose(fp1);
		kseq_destroy(ks1);

		
	}//end for
	cerr << "the totalNumber of sequence is: " << index << endl;

	return 0;
}


