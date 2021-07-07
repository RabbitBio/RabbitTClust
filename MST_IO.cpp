#include "MST_IO.h"
//#include "parameter.h"
using namespace std;

void MST2Cluster(string inputFile, string inputFile1, double threshold)
{
	cout << "input files is MST information and genome informations" << endl;

	fstream fs(inputFile);
	fstream fs1(inputFile1);
	if(!fs){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
		printUsage();
		exit(1);
	}
	if(!fs1){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile1.c_str());
		printUsage();
		exit(1);
	}

	vector<EdgeInfo> mst;
	string line;
	//read the GenomeInfo and MST
	while(getline(fs, line)){
		stringstream ss;
		ss << line;
		int preNode, sufNode;
		double distance;
		ss >> preNode >> sufNode >> distance;

		EdgeInfo tmpEdge;
		tmpEdge.preNode = preNode;
		tmpEdge.sufNode = sufNode;
		tmpEdge.dist = distance;
		//printf("<%d, %d, %lf>\n", tmpEdge.preNode, tmpEdge.sufNode, tmpEdge.dist); 
		mst.push_back(tmpEdge);
	}//end while

	vector<SimilarityInfo> similarityInfos;
	//TODO: transfer the paremeter from the GenomeInfo file.
	if(!getline(fs1, line)){
		fprintf(stderr, "the inputFile: %s is NULL\n", inputFile1.c_str());
		exit(1);
	}
	cerr << "the sketchByFile is: " << line << endl;
	bool sketchByFile = stoi(line);
	cerr << "the sketchByFile is: " << sketchByFile << endl;


	while(1){
		if(!getline(fs1, line)) break;
		stringstream ss;
		ss << line;
		int id, length;
		string name, comment;
		ss >> id >> name >> length;
		
		getline(fs1, comment);
		

		SimilarityInfo tmpSimilarityInfo;
		tmpSimilarityInfo.id = id;
		tmpSimilarityInfo.name = name;
		tmpSimilarityInfo.comment = comment;
		tmpSimilarityInfo.length = length;
		similarityInfos.push_back(tmpSimilarityInfo);
	}

	vector<EdgeInfo> forest = generateForest(mst, threshold);

	vector<vector<int> > cluster = generateCluster(forest, mst.size()+1);

	printResult(cluster, similarityInfos, sketchByFile);

}

void printResult(vector< vector<int> > cluster, vector<SimilarityInfo> similarityInfos, bool sketchByFile)
{
	cerr << "start the output " << endl;
	
	if(sketchByFile)
	{
		for(int i = 0; i < cluster.size(); i++){
			printf("the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++){
				//if(cluster[i][j] > similarityInfos.size()) continue;
				cout << j << '\t' << cluster[i][j] << '\t' << similarityInfos[cluster[i][j]].length << "nt\t" << similarityInfos[cluster[i][j]].name << endl;
				string tmpComment("");
				for(int k = 1; k < similarityInfos[cluster[i][j]].comment.length(); k++){
					if(similarityInfos[cluster[i][j]].comment[k] != '>'){
						tmpComment += similarityInfos[cluster[i][j]].comment[k];	
					}
					else{
						cout << "\t\t\t\t" << tmpComment << endl;
						tmpComment = "";
					}
				}
				if(tmpComment.length() != 0){
					cout << "\t\t\t\t" << tmpComment << endl;
					tmpComment = "";
				}

			}
			cout << endl;
		}
	}//end sketchByFile
	else//sketch by sequence
	{
		for(int i = 0; i < cluster.size(); i++){
			printf("the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++){
				cout << j << '\t' << cluster[i][j] << '\t' << similarityInfos[cluster[i][j]].length << "nt\t" << similarityInfos[cluster[i][j]].name << ' ' << similarityInfos[cluster[i][j]].comment << endl;

			}
		}
	}//end sketchBySequence

}

void saveMST(string inputFile, string sketchFunc, vector<SimilarityInfo> similarityInfos, vector<EdgeInfo> mst, bool sketchByFile)
{
	//save the matching of graph id and genomeInfo 
	cerr << "save the genomeInfo into: " << inputFile+sketchFunc+"GenomeInfo" << endl;
	ofstream ofile;
	ofile.open(inputFile+sketchFunc+"GenomeInfo");
	if(sketchByFile)
		ofile << '1' << endl;
	else
		ofile << '0' << endl;
	for(int i = 0; i < similarityInfos.size(); i++){
		ofile << similarityInfos[i].id << ' ' << similarityInfos[i].name << ' ' << ' ' << similarityInfos[i].length << endl;
		ofile << similarityInfos[i].comment << endl;
	}
	ofile.close();

	//save the mst
	cerr << "save the MSTInfo into: " << inputFile+sketchFunc+"MSTInfo" << endl;
	ofstream ofile1;
	ofile1.open(inputFile+sketchFunc+"MSTInfo");
	for(int i = 0; i < mst.size(); i++){
		//printf("<%d, %d, %lf>\t%s\t%s\n", mst[i].preNode, mst[i].sufNode, mst[i].dist, similarityInfos[mst[i].preNode].name.c_str(), similarityInfos[mst[i].sufNode].name.c_str());
		ofile1 << mst[i].preNode << ' ' << mst[i].sufNode << ' ' << mst[i].dist << endl;
	}
	cout << endl;
	ofile1.close();

}
