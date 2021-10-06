/* simulate.cpp is used to generator subGenomes with very different size to the origin refGenome.
 * parameters:
 * 	numPerRef: number of subGenomes generated.
 * 	numPerRef should not too big(bigger than 3) since that will make the subGenomes without the very different size.
 * The size of subGenomes is generated with a random proportation of range 0.0 to 1.0. 
 * The  generated subGenomes are saved into the same path with the origin refGenome.
 * 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "kseq.h"
#include <zlib.h>
#include <random>

KSEQ_INIT(gzFile, gzread);
using namespace std;


int main(int argc, char* argv[]){
	if(argc < 2)
	{
		cerr <<"need input fileList, exit! " << endl;
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

	//random generator
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	int numPerRef = 3;
	int nrolls = fileList.size();//number of subGenomes generated.
	double randArr[nrolls * numPerRef];
	int globalIndex = 0;
	for(int i = 0; i < nrolls * numPerRef; i++)
	{
		double tmpDouble = distribution(generator);
		randArr[i] = tmpDouble;
		//fprintf(stdout, "%d\t%lf\n", i, tmpDouble);
	}


	for(int i = 0; i < fileList.size(); i++)
	{
		gzFile fp1 = gzopen(fileList[i].c_str(), "r");
		if(!fp1)
		{
			cerr << "cannot open: " << fileList[i] << ", exit!" << endl;
			return 1;
		}
		kseq_t * ks1;
		ks1 = kseq_init(fp1);

		ofstream ofile[numPerRef];
		for(int j = 0; j < numPerRef; j++)
		{
			ofile[j].open(fileList[i] + to_string(j));
		}
		while(1){
			int length = kseq_read(ks1);
			if(length < 0) break;
			string name = ks1->name.s;
			string comment = ks1->comment.s;
			string content = ks1->seq.s;
			for(int j = 0; j < numPerRef; j++)
			{
				string headLine = '>' + name + ' ' + comment;
				ofile[j] << headLine << endl;
				for(int k = 0; k < length*randArr[i*numPerRef + j]; k+=80)
				{
					string tmpContent = content.substr(k, 80);
					ofile[j] << tmpContent << endl;
				}
			}
		}//end while read the file
		for(int j = 0; j < numPerRef; j++)
		{
			ofile[j].close();
		}
		gzclose(fp1);
		kseq_destroy(ks1);
	}


	return 0;
}
