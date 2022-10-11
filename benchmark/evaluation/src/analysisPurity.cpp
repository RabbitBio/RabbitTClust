/* Author: Xiaoming Xu
 * Data: 2022/8/2
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

using namespace std;

//#define LEVEL "family"
#define LEVEL "genus"
//#define LEVEL "species"

void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions);

int main(int argc , char *argv[]){
	string application = argv[0];
	vector<string> args, descriptions;
	args.push_back(application);
	descriptions.push_back("the application name");

	//========= parameters need changing ========
	//The example is with parameters of specific numbers.
	//args is the tutorial names.
	string pwd = "/RabbitTClust/benchmark/evaluation/src/analysisPurity.cpp";
	string dependency = "None";
	string example = application + " purity.accession nodes.dmp outputFile";
	args.push_back("nodes.dmp");
	args.push_back("purity.accession");
	args.push_back("outputFile");
	descriptions.push_back("input file, the taxonomy nodes.dmp file");
	descriptions.push_back("input file, the purity solved result file, from the calPurity file <accession, species-taxid> per line");
	descriptions.push_back("output file, the final output file");

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
	string nodesFile = argv[1];
	string inputFile = argv[2];
	string outputFile = argv[3];
	string line;
	ifstream ifs0(nodesFile);
	unordered_map<int, pair<int, string>> id_preIDRank_map;
	while(getline(ifs0, line)){
		vector<string> vstr;
		boost::split(vstr, line, boost::is_any_of("\t|"), boost::token_compress_on);
		int curId = stoi(vstr[0]);
		int preId = stoi(vstr[1]);
		string rank = vstr[2];
		pair<int, string> p(preId, rank);
		id_preIDRank_map.insert({curId, p});
	}
	ifs0.close();
	cerr << "the size of id_preIDRank_map is: " << id_preIDRank_map.size() << endl;

	ifstream ifs1(inputFile);
	string outputFile0 = outputFile + ".same";
	string outputFile1 = outputFile + ".diff";
	string outputFile2 = outputFile + ".same0";
	ofstream ofs0(outputFile0);
	ofstream ofs1(outputFile1);
	ofstream ofs2(outputFile2);
	ofs0 << "label\taccession\tspecies\tno_rank\tgenus\tfamily\torder" << endl;
	ofs1 << "label\taccession\tspecies\tno_rank\tgenus\tfamily\torder" << endl;
	ofs2 << "label\taccession\tspecies\tno_rank\tgenus\tfamily\torder" << endl;

	//string repAccession;
	vector<string> repAccessionArr;
	vector<string> badAccessionArr;
	unordered_map<string, int> repClass;
	vector<unordered_map<string, int>> badClassArr;

	int index = 0;
	string cmpLevel = LEVEL;

	while(getline(ifs1, line)){
		//cout << index++ << endl;
		bool getGenus = false;
		if(line.length() == 0){ //finish a cluster
			bool isEqual = true;
			for(auto x : badClassArr){
				if(x[cmpLevel] != repClass[cmpLevel]){
					isEqual = false;
					break;
				}
			}
			if(isEqual){
				if(repClass[cmpLevel] != 0){
					for(auto repAccession : repAccessionArr){
						ofs0 << '+' << '\t' << repAccession << '\t';
						ofs0 << repClass["species"] << '\t' << repClass["no_rank"] << '\t' << repClass["genus"] << '\t' << repClass["family"] << '\t' << repClass["order"] << endl;
					}
					for(int i = 0; i < badClassArr.size(); i++){
						unordered_map<string, int> x = badClassArr[i];
						ofs0 << '-' <<'\t' << badAccessionArr[i] << '\t' << x["species"] << '\t' << x["no_rank"] << '\t' << x["genus"] << '\t' << x["family"] << '\t' << x["order"] << endl;
					}
					ofs0 << endl;
				}
				else{
					for(auto repAccession : repAccessionArr){
						ofs2 << '+' << '\t' << repAccession << '\t';
						for(auto x : repClass){
							ofs2 << x.second << "(" << x.first << ")" << '\t';
						}
						ofs2 << endl;
					}
					for(int i = 0; i < badClassArr.size(); i++){
						ofs2 << '-' << '\t' << badAccessionArr[i] << '\t';
						for(auto x : badClassArr[i]){
							ofs2 << x.second << "(" << x.first << ")" << '\t';
						}
						ofs2 << endl;
					}
					ofs2<< endl;
				}
			}
			else{
				for(auto repAccession : repAccessionArr){
					ofs1 << '+' << '\t' << repAccession << '\t';
					ofs1 << repClass["species"] << '\t' << repClass["no_rank"] << '\t' << repClass["genus"] << '\t' << repClass["family"] << '\t' << repClass["order"] << endl;
				}
				bool hasEqual = false;
				for(int i = 0; i < badClassArr.size(); i++){
					unordered_map<string, int> x = badClassArr[i];
					if(x[cmpLevel] != repClass[cmpLevel]){
						ofs1 << '-' << '\t' << badAccessionArr[i] << '\t' << x["species"] << '\t' << x["no_rank"] << '\t' << x["genus"] << '\t' << x["family"] << '\t' << x["order"] << endl;
					}
					else{
						if(!hasEqual){
							if(repClass[cmpLevel] != 0){
								for(auto repAccession : repAccessionArr){
									ofs0 << '+' << '\t' << repAccession << '\t';
									ofs0 << repClass["species"] << '\t' << repClass["no_rank"] << '\t' << repClass["genus"] << '\t' << repClass["family"] << '\t' << repClass["order"] << endl;
								}
							}
							else{
								for(auto repAccession : repAccessionArr){
									ofs2 << '+' << '\t' << repAccession << '\t';
									ofs2 << repClass["species"] << '\t' << repClass["no_rank"] << '\t' << repClass["genus"] << '\t' << repClass["family"] << '\t' << repClass["order"] << endl;
								}
							}
							hasEqual = true;
						}
						if(repClass[cmpLevel] != 0){
							ofs0 << '-' << '\t' << badAccessionArr[i] << '\t' << x["species"] << '\t' << x["no_rank"] << '\t' << x["genus"] << '\t' << x["family"] << '\t' << x["order"] << endl;
						}
						else{
							ofs2 << '-' << '\t' << badAccessionArr[i] << '\t' << x["species"] << '\t' << x["no_rank"] << '\t' << x["genus"] << '\t' << x["family"] << '\t' << x["order"] << endl;
						}
					}
				}//end for loop of badClassArr
				if(hasEqual){
					if(repClass[cmpLevel] != 0)
						ofs0 << endl;
					else
						ofs2 << endl;
				}
				ofs1 << endl;

			}
			unordered_map<string, int>().swap(repClass);
			vector<unordered_map<string, int>>().swap(badClassArr);
			vector<string>().swap(badAccessionArr);
			vector<string>().swap(repAccessionArr);
			continue;
		}

		string accession;
		int curId;
		stringstream ss;
		ss << line;
		ss >> accession >> curId;
		if(id_preIDRank_map.count(curId) == 0){
			cerr << "the id: " << curId  << " is not in the taxonomy" << endl;
			continue;
		}
		if(line[0] != '\t'){
			repAccessionArr.push_back(accession);
			//repAccession = accession;
			if(id_preIDRank_map.count(curId) > 0){
				string rank = id_preIDRank_map[curId].second;
				repClass.insert({rank, curId});
				repClass[rank] = curId;
			}
			while(id_preIDRank_map.count(curId) > 0 && curId != 1){
				curId = id_preIDRank_map[curId].first;
				string rank = id_preIDRank_map[curId].second;
				if(rank == "no rank"){
					if(!getGenus){
						repClass.insert({rank, curId});
						repClass[rank] = curId;
					}
				}
				else if(rank == "genus"){
					getGenus = true;
				}
				repClass.insert({rank, curId});
				repClass[rank] = curId;
			}
		}//end if line[0] != \t
		else{
			badAccessionArr.push_back(accession);
			unordered_map<string, int> curBadClass;
			if(id_preIDRank_map.count(curId) > 0){
				string rank = id_preIDRank_map[curId].second;
				curBadClass.insert({rank, curId});
				curBadClass[rank] = curId;
			}
			while(id_preIDRank_map.count(curId) > 0 && curId != 1){
				curId = id_preIDRank_map[curId].first;
				string rank = id_preIDRank_map[curId].second;
				if(rank == "no rank"){
					if(!getGenus){
						curBadClass.insert({rank, curId});
						curBadClass[rank] = curId;
					}
				}
				else if(rank == "genus"){
					getGenus = true;
				}
				curBadClass.insert({rank, curId});
				curBadClass[rank] = curId;
			}
			badClassArr.push_back(curBadClass);
			unordered_map<string, int>().swap(curBadClass);
		}

	}
	ifs1.close();
	ofs0.close();
	ofs1.close();
	ofs2.close();
	

	cerr << "finished" << endl;
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
