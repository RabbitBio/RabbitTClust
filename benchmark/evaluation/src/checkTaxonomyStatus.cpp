/* Author: Xiaoming Xu
 * Data: 2022/8/3
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

void printInfo(string pwd, string dependency, string example, vector<string> args, vector<string> descriptions);

int main(int argc , char *argv[]){
	string application = argv[0];
	vector<string> args, descriptions;
	args.push_back(application);
	descriptions.push_back("the application name");

	//========= parameters need changing ========
	//The example is with parameters of specific numbers.
	//args is the tutorial names.
	string pwd = "RabbitTClust/benchmark/evaluation/src/checkTaxonomyStatus.cpp";
	string dependency = "None";
	string example = application + " ANI_report_prokaryotes.txt greedy.ans output";
	args.push_back("ANI_report_prokaryotes.txt");
	args.push_back("greedy.ans");
	args.push_back("output");
	descriptions.push_back("input file, the ANI_report_prokaryotes file, <genbank-accession, species-taxid, best-match-species-taxid, best-match-status, excluded-from-refseq> per line");
	descriptions.push_back("input file, the analysis from the analysisPurity, <index(+/-), genbank-accession, species-taxid, no-rank, genus-taxid, family-taxid, order-taxid> per line");
	descriptions.push_back("output file, the output file of result");

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
	string aniFile = argv[1];
	string anaFile = argv[2];
	string outputFile = argv[3];

	string line;
	unordered_set<string> accessionSet;
	unordered_map<string, int> accSpeciesTaxidMap;
	unordered_map<string, string> accExcludedFromRefseqMap;
	unordered_map<string, int> accBestMatchSpeciesTaxidMap;
	unordered_map<string, string> accBestMatchStatusMap;
	unordered_map<string, double> accQcoverageMap;
	unordered_map<string, double> accScoverageMap;
	ifstream ifs0(aniFile);
	getline(ifs0, line);//header line
	int index = 0;
	int calNumBestMatchStatus = 0;
	int calNumExcludedFromRefseq = 0;
	while(getline(ifs0, line)){
		string accession, excluded_from_refseq, best_match_status;
		int species_taxid,  best_match_species_taxid;
		double qcoverage, scoverage;
		vector<string> vstr;
		boost::split(vstr, line, boost::is_any_of("\t"), boost::token_compress_on);
		accession = vstr[0];
		species_taxid = vstr[1] != "na" ? stoi(vstr[1]) : 0;
		best_match_species_taxid = vstr[2] != "na" ? stoi(vstr[2]) : 0;
		best_match_status = vstr[3];
		excluded_from_refseq = vstr[4];
		qcoverage = vstr[5] != "na" ? stod(vstr[5]) : 0.0;
		scoverage = vstr[6] != "na" ? stod(vstr[6]) : 0.0;
		if(best_match_status == "species-match")	calNumBestMatchStatus++;
		if(excluded_from_refseq == "na"){
			calNumExcludedFromRefseq++;
			//cout << accession << endl;
		}

		accSpeciesTaxidMap.insert({accession, species_taxid});
		accExcludedFromRefseqMap.insert({accession, excluded_from_refseq});
		accBestMatchSpeciesTaxidMap.insert({accession, best_match_species_taxid});
		accBestMatchStatusMap.insert({accession, best_match_status});
		accQcoverageMap.insert({accession, qcoverage});
		accScoverageMap.insert({accession, scoverage});
		accessionSet.insert({accession});
	}
	ifs0.close();
	cerr << "the size of accSpeciesTaxidMap is: " << accSpeciesTaxidMap.size() << endl;
	cerr << "the best_match_status of species_match is: " << calNumBestMatchStatus << ", the percent is: " << (double)calNumBestMatchStatus / accSpeciesTaxidMap.size() << endl;
	cerr << "the excluded_from_refseq of na is: " << calNumExcludedFromRefseq << ", the percent is: " << (double)calNumExcludedFromRefseq / accSpeciesTaxidMap.size() << endl;
	//exit(0);

	ifstream ifs1(anaFile);
	string outputFile0 = outputFile + ".species_taxid.check";
	string outputFile1 = outputFile + ".best_match_species_taxid.check";
	string outputFile2 = outputFile + ".exclude_from_refseq.check";
	string outputFile3 = outputFile + ".best_match_status.check";
	string outputFile4 = outputFile + ".perfect.check";
	string outputFile5 = outputFile + ".coverage.check";
	ofstream ofs0(outputFile0);
	ofstream ofs1(outputFile1);
	ofstream ofs2(outputFile2);
	ofstream ofs3(outputFile3);
	ofstream ofs4(outputFile4);
	ofstream ofs5(outputFile5);
	ofs0 << "label\taccession\tassembly_taxid\ttaxonomy_taxid" << endl;
	ofs1 << "label\taccession\tassembly_taxid\tbest_match_species_taxid" << endl;
	ofs2 << "label\taccession\texclude_from_refseq" << endl;
	ofs3 << "label\taccession\tbest_match_status" << endl;
	ofs4 << "label\taccession\tassembly_taxid" << endl;
	ofs5 << "label\taccession\tqcoverage\tscoverage" << endl;
	getline(ifs1, line);//header line
	int numNotInTaxonomy = 0;
	int totalNumber_rep = 0;
	int totalNumber_bad = 0;
	int perfectNum_rep = 0;
	int perfectNum_bad = 0;
	int numNotEqualTaxid_rep = 0;
	int numNotEqualTaxid_bad = 0;
	int numNotEqualBestMatch_rep = 0;
	int numNotEqualBestMatch_bad = 0;
	int numNotBestMatchSpeciesLevel_rep = 0;
	int numNotBestMatchSpeciesLevel_bad = 0;
	int numExclude_from_refSeq_rep = 0;
	int numExclude_from_refSeq_bad = 0;

	vector<string> matchStatusArr;
	matchStatusArr.push_back("species-match");
	matchStatusArr.push_back("subspecies-match");
	matchStatusArr.push_back("synonym-match");
	matchStatusArr.push_back("derived-species-match");
	matchStatusArr.push_back("genus-match");
	matchStatusArr.push_back("approved-mismatch");
	matchStatusArr.push_back("mismatch");
	matchStatusArr.push_back("below-threshold-match");
	matchStatusArr.push_back("below-threshold-mismatch");
	matchStatusArr.push_back("low-coverage");
	
	vector<int> numMatchStatusRepArr;
	vector<int> numMatchStatusBadArr;
	for(int i = 0; i < 10; i++){
		numMatchStatusRepArr.push_back(0);
		numMatchStatusBadArr.push_back(0);
	}

	//int rep_numSpeciesMatch = 0;
	//int rep_numSubSpeciesMatch = 0;
	//int rep_numSynonymMatch = 0;
	//int rep_numDerivedSpeciesMatch = 0;
	//int rep_numGenusMatch = 0;
	//int rep_numApprovedMismatch = 0;
	//int rep_numMismatch_rep = 0;
	//int rep_numBelowThresholdMatch = 0;
	//int rep_numBelowThresholdMismatch = 0;
	//int rep_numLowCoverage = 0;

	//int bad_numSpeciesMatch = 0;
	//int bad_numSubSpeciesMatch = 0;
	//int bad_numSynonymMatch = 0;
	//int bad_numDerivedSpeciesMatch = 0;
	//int bad_numGenusMatch = 0;
	//int bad_numApprovedMismatch = 0;
	//int bad_numMismatch = 0;
	//int bad_numBelowThresholdMatch = 0;
	//int bad_numBelowThresholdMismatch = 0;
	//int bad_numLowCoverage = 0;

	while(getline(ifs1, line)){
		if(line.length() == 0){
			ofs0 << endl;
			ofs1 << endl;
			ofs2 << endl;
			ofs3 << endl;
			//ofs4 << endl;
			ofs5 << endl;
			//cout << endl;
			continue;
		}
		string index, accession;
		int species, no_rank, genus, family, order;
		stringstream ss;
		ss << line;
		ss >> index >> accession >> species >> no_rank >> genus >> family >> order;
		if(accessionSet.count(accession) == 0){
			//cerr << "the accession is not in the taxonomy file" << endl;
			//cout << line << endl;
			numNotInTaxonomy++;
			continue;
		}

		int taxSID = accSpeciesTaxidMap[accession];
		int taxBMID = accBestMatchSpeciesTaxidMap[accession];
		string taxEFR = accExcludedFromRefseqMap[accession];
		string taxBMS = accBestMatchStatusMap[accession];
		double qcoverage = accQcoverageMap[accession];
		double scoverage = accScoverageMap[accession];
		if(index == "+"){
			totalNumber_rep++;
			if(species != taxSID) numNotEqualTaxid_rep++;
			if(taxSID != taxBMID) numNotEqualBestMatch_rep++;
			if(taxEFR != "na") numExclude_from_refSeq_rep++;
			if(taxBMS != "species-match") numNotBestMatchSpeciesLevel_rep++;

			for(int i = 0; i < 10; i++){
				if(taxBMS == matchStatusArr[i])	numMatchStatusRepArr[i]++;
			}

			if(species == taxSID && taxSID == taxBMID && taxEFR == "na" && taxBMS == "species-match"){
				perfectNum_rep++;
				ofs4 << line << endl;
			}
		}
		else{
			totalNumber_bad++;
			if(species != taxSID) numNotEqualTaxid_bad++;
			if(taxSID != taxBMID) numNotEqualBestMatch_bad++;
			if(taxEFR != "na") numExclude_from_refSeq_bad++;
			if(taxBMS != "species-match") numNotBestMatchSpeciesLevel_bad++;
			
			for(int i = 0; i < 10; i++){
				if(taxBMS == matchStatusArr[i])	numMatchStatusBadArr[i]++;
			}

			//if(taxBMS == "mismatch"|| taxBMS == "below-threshold-match" || taxBMS == "below-threshold-mismatch" || taxBMS == "low-coverage") numNotBestMatchSpeciesLevel_bad++;
			//if(species == taxSID && taxSID == taxBMID && taxEFR == "na" && taxBMS == "species-match"){
			//if(species == taxSID &&  taxEFR == "na" && taxBMS == "species-match"){
			if(species == taxBMID && taxEFR == "na"){
				perfectNum_bad++;
				ofs4 << line << endl;
			}
		}

		ofs0 << index << '\t' << accession << '\t' << species << '\t' << taxSID << endl;
		ofs1 << index << '\t' << accession << '\t' << species << '\t' << taxBMID << endl;
		ofs2 << index << '\t' << accession << '\t' << taxEFR << endl;
		ofs3 << index << '\t' << accession << '\t' << taxBMS << endl;
		ofs5 << index << '\t' << accession << '\t' << qcoverage << '\t' << scoverage << endl;
	}
	ofs0.close();
	ofs1.close();
	ofs2.close();
	ofs3.close();
	ofs5.close();

	cerr << "finished" << endl;
	cerr << "the number not in the taxonomy is: " << numNotInTaxonomy << endl;
	cerr << "for representative genomes, the total number is: " << totalNumber_rep << endl;
	cerr << "\tthe numNotEqualtaxid of assembly-summary and taxonomy is: " << numNotEqualTaxid_rep << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotEqualTaxid_rep / totalNumber_rep << endl;
	cerr << "\tthe numNotEqualBestMatch of species-taxid and best-species-taxid is: " << numNotEqualBestMatch_rep << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotEqualBestMatch_rep / totalNumber_rep << endl;
	cerr << "\tthe numExclude_from_refSeq_rep is: " << numExclude_from_refSeq_rep << endl;
	cerr << "\t\tthe percentange is: " << (double)numExclude_from_refSeq_rep / totalNumber_rep << endl;
	cerr << "\tthe numNotBestMatchSpeciesLevel_rep is: " << numNotBestMatchSpeciesLevel_rep << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotBestMatchSpeciesLevel_rep / totalNumber_rep << endl;
	cerr << "\tthe perfectNum_rep is: " << perfectNum_rep << endl;
	cerr << "\t\tthe percentange is: " << (double)perfectNum_rep / totalNumber_rep << endl;

	cerr << "for bad genomes, the total number is: " << totalNumber_bad << endl;
	cerr << "\tthe numNotEqualtaxid of assembly-summary and taxonomy is: " << numNotEqualTaxid_bad << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotEqualTaxid_bad / totalNumber_bad << endl;
	cerr << "\tthe numNotEqualBestMatch of species-taxid and best-species-taxid is: " << numNotEqualBestMatch_bad << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotEqualBestMatch_bad / totalNumber_bad << endl;
	cerr << "\tthe numExclude_from_refSeq_bad is: " << numExclude_from_refSeq_bad << endl;
	cerr << "\t\tthe percentange is: " << (double)numExclude_from_refSeq_bad / totalNumber_bad << endl;
	cerr << "\tthe numNotBestMatchSpeciesLevel_bad is: " << numNotBestMatchSpeciesLevel_bad << endl;
	cerr << "\t\tthe percentange is: " << (double)numNotBestMatchSpeciesLevel_bad / totalNumber_bad << endl;
	cerr << "\tthe perfectNum_bad is: " << perfectNum_bad << endl;
	cerr << "\t\tthe percentange is: " << (double)perfectNum_bad / totalNumber_bad << endl;


	cerr << "====================================================================================" << endl;
	int tmpGoodTotalNumber = 0;
	for(int i = 0; i < 10; i++){
		cerr << "the number of rep " << matchStatusArr[i] << " is: " << numMatchStatusRepArr[i] << ", and percent is: " << (double)numMatchStatusRepArr[i] / totalNumber_rep << endl;
		tmpGoodTotalNumber += numMatchStatusRepArr[i];
	}
	cerr << "the total good number is: " << tmpGoodTotalNumber << endl;

	cerr << "====================================================================================" << endl;
	int tmpBadTotalNumber = 0;	
	for(int i = 0; i < 10; i++){
		cerr << "the number of bad " << matchStatusArr[i] << " is: " << numMatchStatusBadArr[i] << ", and percent is: " << (double)numMatchStatusBadArr[i] / totalNumber_bad << endl;
		tmpBadTotalNumber += numMatchStatusBadArr[i];
	}
	cerr << "the total bad number is: " << tmpBadTotalNumber << endl;


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
