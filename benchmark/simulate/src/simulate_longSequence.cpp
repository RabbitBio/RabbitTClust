/* Author: Xiaoming Xu
 * Data: 2022/5/12
 * 
 * See the LICENSE.txt file included with this software for licence information.
 */
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <assert.h>

using namespace std;
void printInfo(string pwd, string example, vector<string> args, vector<string> descriptions){
	assert(args.size() == descriptions.size());
	cerr << endl;
	cerr << "example: " << example << endl;
	cerr << endl;
	cerr << "source file path: " << pwd << endl;
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


int main(int argc , char *argv[]){
	string pwd = "RabbitTClust/benchmark/simulate/src/simulate_longSequence.cpp";
	string application = argv[0];
	string example = application + " 10 20 300 1000000 simulate_10_20_300_1M";
	vector<string> args, descriptions;
	args.push_back(application);
	args.push_back("mutation_rate*1000(integer)");
	args.push_back("numSeedSeqs");
	args.push_back("numEachClusts");
	args.push_back("seqLength");
	args.push_back("output");
	descriptions.push_back("the application name");
	descriptions.push_back("the mutation rate");
	descriptions.push_back("the number of seed sequence (number of clusters)");
	descriptions.push_back("the number sequence in a cluster generate from each seed sequence");
	descriptions.push_back("the approximate length for each sequence");
	descriptions.push_back("the prefix name for groundTruth, seedSequences and totalSimulateSequence fasta files");

	assert(args.size() == descriptions.size());

  if(argc != args.size()) {
		printInfo(pwd, example, args, descriptions);
    return -1;
  }
	else if(argc == 2 && (argv[1] == "-h" || argv[1] == "--help"))
	{
		printInfo(pwd, example, args, descriptions);
		return 1;
	}

  const char nucs[4] = { 'A','T','G','C'};

  int erate = std::atoi(argv[1]);
	int numClusts = std::stoi(argv[2]);
	int numEachClusts = std::stoi(argv[3]);
	int seqLength = std::stoi(argv[4]);
	string outPrefix = argv[5];
	string outSeedFile = outPrefix + "_seed.fna";
	string outTotalFile = outPrefix + "_total.fna";
	string outGroundTruth = outPrefix + "_groundTruth";

	cerr << "the error rate is: " << double(erate)/1000 << endl;
	cerr << "the number of clusters is: " << numClusts << endl;
	cerr << "the number of sequences in each cluster is: " << numEachClusts << endl;
	cerr << "the approximate sequence length is: " << seqLength << endl;
	cerr << "the output seed sequences file is: " << outSeedFile << endl;
	cerr << "the output total sequences file is: " << outTotalFile << endl;
	cerr << "the groundTruth file is: " << outGroundTruth << endl;

	FILE * fp0 = fopen(outSeedFile.c_str(), "w");
	FILE * fp1 = fopen(outTotalFile.c_str(), "w");
	FILE * fp2 = fopen(outGroundTruth.c_str(), "w");

	string key1 = "seqName";
	string key2 = "taxid";
	fprintf(fp2, "%s\t%s\n",key1.c_str(), key2.c_str());
	
	for(int i = 0; i < numClusts; i++)
	{
		struct timeval tv;
		gettimeofday(&tv, NULL);
  	srand(tv.tv_usec);
		string seqName = ">seq_" + to_string(i);
		string groundTruthName = seqName.substr(1);
		fprintf(fp2, "%s\t%d\n", groundTruthName.c_str(), i);
		string seqComment = "Seed sequence " + to_string(i) + " to generate mutations";
		string infoLine = seqName + '\t' + seqComment + '\n';
		string seedSeq("");
		for(int i = 0; i < seqLength; i++)
		{
			char newC = nucs[random()%4];
			seedSeq += newC;
		}
		//output the seed sequence into outSeedFile
		int infoLineLen = infoLine.length();
		fwrite(infoLine.c_str(), sizeof(char), infoLineLen, fp0);
		fwrite(infoLine.c_str(), sizeof(char), infoLineLen, fp1);
		int seedLen = seedSeq.length();
		string outSeedSeq("");
		for(int k = 0; k < seedLen; k += 80)
		{
			int curLength = std::min(80, seedLen-k);
			string tmpLine = seedSeq.substr(k, curLength);
			outSeedSeq += tmpLine + '\n';
		}
		int outSeedSeqLen = outSeedSeq.length();
		fwrite(outSeedSeq.c_str(), sizeof(char), outSeedSeqLen, fp0);
		fwrite(outSeedSeq.c_str(), sizeof(char), outSeedSeqLen, fp1);


		//for generate mutation sequences
		for(int j = 0; j < numEachClusts; j++)
		{
			string mutationName = seqName + "_mutation_" + to_string(j);
			string groundTruthMuName = mutationName.substr(1);
			fprintf(fp2, "%s\t%d\n", groundTruthMuName.c_str(), i);
			string mutationComment = "mutation sequence " + to_string(j) + " from seedSequence " + to_string(i);
			string mutaInfoLine = mutationName + '\t' + mutationComment + '\n';
			string mutationSeq("");
			for(int t = 0; t < seedSeq.length(); t++)
			{
				if(random()%1000 < erate){
					int mut = random()%3;
					if(mut == 0)//sub
					{
						while(1){
							char newc = nucs[random()%4];
							if(newc != seedSeq[t]){
								mutationSeq += newc;
								break;
							}
						}
					}
					else if(mut == 1)// ins
					{
						mutationSeq += nucs[random()%4];
						t = t - 1;
					}
					else//del
						continue;
				}//end if mutation
				else// no mutation
					mutationSeq += seedSeq[t];
			}
			int mutaInfoLineLen = mutaInfoLine.length();
			fwrite(mutaInfoLine.c_str(), sizeof(char), mutaInfoLineLen, fp1);
			int mutationLen = mutationSeq.length();
			string outMutationSeq("");
			for(int k = 0; k < mutationLen; k += 80)
			{
				int curLength = std::min(80, mutationLen-k);
				string tmpLine = mutationSeq.substr(k, curLength);
				outMutationSeq += tmpLine + '\n';
			}
			int outMutationSeqLen = outMutationSeq.length();
			fwrite(outMutationSeq.c_str(), sizeof(char), outMutationSeqLen, fp1);
		}
	}
	fclose(fp0);
	fclose(fp1);
	fclose(fp2);

	cerr << "finish generate mutation files with multithread " << endl;

  return 0;
}

