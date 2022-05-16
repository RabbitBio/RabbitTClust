#include "MST_IO.h"
#include <ctime>
using namespace std;


inline bool cmpSketchLength(ClusterInfo c1, ClusterInfo c2){
	return c1.length > c2.length;
}

void MST2Cluster(string inputFile, string inputFile1, string outputFile, double threshold)
{
	cerr << "input files is MST information and genome informations" << endl;

	fstream fs0(inputFile);//Genome Info
	fstream fs1(inputFile1);//MST Info
	if(!fs0){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
		printUsage();
		exit(1);
	}
	if(!fs1){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile1.c_str());
		printUsage();
		exit(1);
	}


	string line;

	//read the fileSeqs.info or seq.info
	//TODO: transfer the paremeter from the GenomeInfo file.
	if(!getline(fs0, line)){
		fprintf(stderr, "the inputFile: %s is NULL\n", inputFile1.c_str());
		exit(1);
	}

	//sketchByfile is 1 else is 0
	bool sketchByFile = stoi(line);
	cerr << "the sketchByFile is: " << sketchByFile << endl;

	vector<SketchInfo> sketches;

	if(sketchByFile)
	{
		while(1)
		{
			if(!getline(fs0, line)) break;
			string fileName, name, comment, tmp;
			int strand, totalLength;
			stringstream ss;
			ss << line;
			ss >> fileName >> name >> strand >> totalLength;
			while(ss >> tmp) 
				comment += tmp + ' ';
			comment = comment.substr(0, comment.length()-1);
			Vec_SeqInfo curFileSeqs;
			SequenceInfo curSeq{name, comment, strand, totalLength};
			curFileSeqs.push_back(curSeq);
			SketchInfo curSketch;
			curSketch.fileName = fileName;
			curSketch.totalSeqLength = totalLength;
			curSketch.fileSeqs = curFileSeqs;
			sketches.push_back(curSketch);
		}
	}
	else//sketch sequences
	{
		while(1)
		{
			if(!getline(fs0, line)) break;
			string name, comment, tmp;
			int strand, length;
			stringstream ss;
			ss << line;
			ss >> name >> strand >> length;
			while (ss >> tmp)
				comment += tmp + ' ';
			comment = comment.substr(0, comment.length()-1);
			SequenceInfo curSeq{name, comment, strand, length};
			SketchInfo curSketch;
			curSketch.seqInfo = curSeq;
			sketches.push_back(curSketch);
		}
	}

	vector<EdgeInfo> mst;
	//read the MSTInfo and GenomeInfo
	getline(fs1, line);//sketchByFile or not
	cerr << line << endl;
	getline(fs1, line);//Sketch Function
	cerr << line << endl;
	getline(fs1, line);//SketchSize
	cerr << line << endl;
	getline(fs1, line);//kmerSize
	cerr << line << endl;

	while(getline(fs1, line)){
		stringstream ss;
		ss << line;
		int preNode, sufNode;
		double distance;
		ss >> preNode >> sufNode >> distance;

		EdgeInfo tmpEdge{preNode, sufNode, distance};
		mst.push_back(tmpEdge);
	}//end while


	vector<EdgeInfo> forest = generateForest(mst, threshold);

	vector<vector<int> > cluster = generateCluster(forest, mst.size()+1);

	printResult(cluster, sketches, sketchByFile, outputFile);

}

void printResult(vector< vector<int> > cluster, vector<SketchInfo> sketches, bool sketchByFile, string outputFile)
{
	cerr << "output the result into: " << outputFile << endl;
	FILE *fp = fopen(outputFile.c_str(), "w");
	
	if(sketchByFile)
	{
		for(int i = 0; i < cluster.size(); i++){
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j];
				fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileName.c_str(),  sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());

			}
			fprintf(fp, "\n");
		}
	}//end sketchByFile

	else//sketch by sequence
	{
		for(int i = 0; i < cluster.size(); i++){
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j];		
				fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].seqInfo.length, sketches[curId].seqInfo.name.c_str(), sketches[curId].seqInfo.comment.c_str());
			}
			fprintf(fp, "\n");
		}
	}//end sketchBySequence
	fclose(fp);

}

void saveMST(string inputFile, string sketchFunc, bool isContainment, int containCompress, vector<SketchInfo> sketches, vector<EdgeInfo> mst, bool sketchByFile, int sketchSize, int kmerSize)
{
	//save the matching of graph id and genomeInfo 
	string folderPath = currentDataTime();
	string command = "mkdir -p " + folderPath;
	system(command.c_str());
	string prefixName = inputFile;
	std::size_t found = prefixName.find_last_of('/');//if not found, return std::string::npos(-1);
	prefixName = prefixName.substr(found+1);

	cerr << "save the genomeInfo into: " << folderPath << '/' << prefixName+sketchFunc+"GenomeInfo" << endl;
	ofstream ofile;
	ofile.open(folderPath + '/' + prefixName+sketchFunc+"GenomeInfo");
	if(sketchByFile)
		ofile << '1' << endl;
	else
		ofile << '0' << endl;
	
	if(sketchByFile)
	{
		for(int i = 0; i < sketches.size(); i++)
		{
			Vec_SeqInfo curFileSeqs = sketches[i].fileSeqs;
			ofile << sketches[i].fileName << ' ' << curFileSeqs[0].name << ' ' << curFileSeqs[0].strand << ' ' << sketches[i].totalSeqLength << ' ' << curFileSeqs[0].comment << '\n';
		}
	}
	else//sketchBySequence
	{
		for(int i = 0; i < sketches.size(); i++)
		{
			SequenceInfo curSeq = sketches[i].seqInfo;
			ofile << curSeq.name << ' ' << ' ' << curSeq.strand << ' ' << curSeq.length << curSeq.comment << endl;
		}

	}
	ofile.close();

	//save the mst
	cerr << "save the MSTInfo into: " << folderPath << '/' << prefixName+sketchFunc+"MSTInfo" << endl;
	ofstream ofile1;
	ofile1.open(folderPath + '/' + prefixName+sketchFunc+"MSTInfo");

	if(sketchByFile)
		ofile1 << "sketch by File! " << endl;
	else 
		ofile1 << "sketch by Sequence!" << endl;

	string tmpFunc = sketchFunc;
	if(sketchFunc == "MinHash")
	{
	if(isContainment)
		tmpFunc += " containment";
	else
		tmpFunc += " mashDistance";
	}
	ofile1 << "the sketch function is: " << tmpFunc << endl;
	if(isContainment)
		ofile1 << "The sketchSize is in proportion with 1/" << containCompress << endl;
	else
		ofile1 << "The sketchSize is: " << sketchSize << endl;
	ofile1 << "The kmerSize is: " << kmerSize << endl;
	for(int i = 0; i < mst.size(); i++){
		ofile1 << mst[i].preNode << ' ' << mst[i].sufNode << ' ' << mst[i].dist << endl;
	}
	cout << endl;
	ofile1.close();

}
