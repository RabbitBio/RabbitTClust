/* The evaluation of cluster is based on the precision and recall.
 * CalResult1.cpp is used to calculate the evaluation of precision and recall.
 * The label of genome is as the first keyword of nomenclature of gene feature.
 * All the concepts refer to https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html 
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>

using namespace std;

double getPurity(vector< vector<uint64_t> > numClust)
{
	uint64_t totalElement = 0;
	uint64_t purityElement = 0;
	for(uint64_t i = 0; i < numClust.size(); i++)
	{
		uint64_t tmpPurity = 0;
		for(uint64_t j = 0; j < numClust[i].size(); j++)
		{
			totalElement += numClust[i][j];
			tmpPurity = std::max(tmpPurity, numClust[i][j]);
		}
		purityElement += tmpPurity;
	}

	cerr << "the purityElement is: " << purityElement << endl;
	cerr << "the totalElement is: " << totalElement << endl;
	double result = (double) purityElement / totalElement;
	return result;
}

inline uint64_t combin(uint64_t x)
{
	if(x < 2) return 0;
	return (x * (x-1))/2;
}

uint64_t getTP(vector< vector<uint64_t> > numClust)
{
	uint64_t numTP = 0;
	for(uint64_t i = 0; i < numClust.size(); i++)
	{
		for(uint64_t j = 0; j < numClust[i].size(); j++)
		{
			numTP += combin(numClust[i][j]);
		}
	}
	return numTP;
}

uint64_t getFP(vector< vector<uint64_t> > numClust)
{
	uint64_t numTP_FP = 0;
	for(uint64_t i = 0; i < numClust.size(); i++)
	{
		uint64_t clustElement = 0;
		for(uint64_t j = 0; j < numClust[i].size(); j++)
		{
			clustElement += numClust[i][j];
		}
		numTP_FP += combin(clustElement);
	}
	uint64_t numTP = getTP(numClust);
	uint64_t numFP = numTP_FP - numTP;
	return numFP;
}

uint64_t getFN(vector< vector<uint64_t> > numClust, unordered_map<string, uint64_t> mapClust)
{
	uint64_t numTP_FN = 0;
	for(auto x : mapClust)
	{
		numTP_FN += combin(x.second);
	}
	uint64_t numTP = getTP(numClust);
	uint64_t numFN = numTP_FN - numTP;
	return numFN;
}

uint64_t getTN(vector< vector<uint64_t> > numClust, unordered_map<string, uint64_t> mapClust)
{
	uint64_t totalElement = 0;
	for(auto x : mapClust)
	{
		totalElement += x.second;
	}
	uint64_t numTP_FP_FN_TN = combin(totalElement);
	uint64_t numTP = getTP(numClust);
	uint64_t numFP = getFP(numClust);
	uint64_t numFN = getFN(numClust, mapClust);
	uint64_t numTN = numTP_FP_FN_TN - numTP - numFP - numFN;
	return numTN;
}

double getRI(vector< vector<uint64_t> > numClust, unordered_map<string, uint64_t> mapClust)
{
	uint64_t tp = getTP(numClust);
	uint64_t fp = getFP(numClust);
	uint64_t fn = getFN(numClust, mapClust);
	uint64_t tn = getTN(numClust, mapClust);
	cerr << "tp: " << tp << endl;
	cerr << "fp: " << fp << endl;
	cerr << "fn: " << fn << endl;
	cerr << "tn: " << tn << endl;
	double ri = (double) (tp + tn) / (tp + fp + fn + tn);
	double precision = (double)tp / (tp + fp);
	double recall = (double)tp / (tp + fn);
	double ari = (double)2 * (tp*tn - fn*fp) / ((tp+fn)*(fn+tn) + (tp+fp)*(fp+tn));
	double beta = 1.0;
	double sqrBeta = beta * beta;
	double F1 = (1 + sqrBeta) * (precision * recall) / (sqrBeta *precision + recall);
	cout << "the precision is: " << precision << endl;
	cout << "the recall is: " << recall << endl;
	cout << "the ari is: " << ari << endl;
	cout << "the F1(F beta) is: " << F1 << endl;
	return ri;
}


int main(int argc, char* argv[]){
	if(argc < 3){
		cerr << "run with: ./calResult1 -l(-i) result.out" << endl;
		cerr << "where -l means cluster by genomes, -i means cluster by sequences" << endl;
		return 1;
	}
	string argument = argv[1];
	string inputFile = argv[2];

	fstream fs(inputFile);
	string line;
	unordered_map<string, uint64_t> mapClust;
	vector< vector<uint64_t> > numClust;
	unordered_map<string, uint64_t> tmpMap;
	while(getline(fs, line))
	{
		if(line.length() == 0) continue;
		if(line[0] != '\t')//end the last cluster
		{
			if(tmpMap.size() != 0)
			{
				vector<uint64_t> tmpClust;
				for(auto x : tmpMap)
				{
					tmpClust.push_back(x.second);
				}
				unordered_map<string, uint64_t>().swap(tmpMap);
				numClust.push_back(tmpClust);
				vector<uint64_t>().swap(tmpClust);
			}
		}
		else
		{
			stringstream ss;
			ss << line;
			int curId, genomeId;
			string genomeSize, fileName, genomeName;
			string type0, type1, type2;
			if(argument == "-l")
				ss >> curId >> genomeId >> genomeSize >> fileName >> genomeName >> type0 >> type1 >> type2;
			else if(argument == "-i")
				ss >> curId >> genomeId >> genomeSize >>  genomeName >> type0 >> type1 >> type2;
			else
			{
				cerr << "err input argument: " << argument <<endl;
				return 1;
			}
			if(type0.substr(0, 10) == "UNVERIFIED")
			{
				type0 = type1;
				type1 = type2;
			}
			if(type0.back() == ',') type0.pop_back();
			if(type1.back() == ',') type1.pop_back();

			string key = type0;
			tmpMap.insert({key, 0});
			tmpMap[key]++;
			mapClust.insert({key, 0});
			mapClust[key]++;
		}

	}//end while

	double result = getPurity(numClust);
	cout << "the purity is: " << result << endl;
	double ri = getRI(numClust, mapClust);
	cout << "the ri is: " << ri << endl;

	//for(uint64_t i = 0; i < numClust.size(); i++)
	//{
	//	for(uint64_t j = 0; j < numClust[i].size(); j++)
	//	{
	//		cout << numClust[i][j] << '\t';
	//	}
	//	cout << endl;
	//}


	return 0;
}
