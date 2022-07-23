#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>

using namespace std;

void getGroundTruthBySequence(string groundTruth, unordered_map<string, int>& seqName_taxid_map, unordered_map<int, string>& taxid_organismName_map);


void getGroundTruthByFile(string groundTruth, unordered_map<string, int>& accession_taxid_map, unordered_map<int, string>& taxid_organismName_map);







