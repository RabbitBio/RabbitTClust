#include "MST_IO.h"
#include <ctime>
//#include "parameter.h"
using namespace std;

const string currentDataTime(){
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y_%m_%d_%H-%M-%S", &tstruct);

	return buf;
}

inline bool cmpSketchLength(ClusterInfo c1, ClusterInfo c2){
	return c1.length > c2.length;
}

void MST2Cluster(string inputFile, string inputFile1, string outputFile, double threshold)
{
	cout << "input files is MST information and genome informations" << endl;

	fstream fs(inputFile);//MST Info
	fstream fs1(inputFile1);//Genome Info
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
	//read the MSTInfo and GenomeInfo
	getline(fs, line);//sketchByFile or not
	cerr << line << endl;
	getline(fs, line);//Sketch Function
	cerr << line << endl;
	getline(fs, line);//SketchSize
	cerr << line << endl;
	getline(fs, line);//kmerSize
	cerr << line << endl;

	while(getline(fs, line)){
		stringstream ss;
		ss << line;
		int preNode, sufNode;
		double distance;
		ss >> preNode >> sufNode >> distance;

		EdgeInfo tmpEdge{preNode, sufNode, distance};
		mst.push_back(tmpEdge);
	}//end while

	//read the fileSeqs.info or seq.info
	//TODO: transfer the paremeter from the GenomeInfo file.
	if(!getline(fs1, line)){
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
			if(!getline(fs1, line)) break;
			int sequenceNumber = stoi(line);
			string fileName;
			getline(fs1, fileName);
			Vec_SeqInfo curFileSeqs;
			uint64_t totalLength = 0;
			//this is the sequences info in one file
			for(int i = 0; i < sequenceNumber; i++)
			{
				getline(fs1, line);
				string name, comment;
				int strand, length;
				stringstream ss;
				ss << line;
				ss >> name >> strand >> length;
				getline(fs1, comment);
				totalLength += length;
				SequenceInfo curSeq{name, comment, strand, length};
				curFileSeqs.push_back(curSeq);
			}
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
			if(!getline(fs1, line)) break;
			string name, comment;
			int strand, length;
			stringstream ss;
			ss << line;
			ss >> name >> strand >> length;
			getline(fs1, comment);
			SequenceInfo curSeq{name, comment, strand, length};
			SketchInfo curSketch;
			curSketch.seqInfo = curSeq;
			sketches.push_back(curSketch);
		}
	}

	vector<EdgeInfo> forest = generateForest(mst, threshold);

	vector<vector<int> > cluster = generateCluster(forest, mst.size()+1);

	printResult(cluster, sketches, sketchByFile, outputFile);

}

void printResult(vector< vector<int> > clusterOrigin, vector<SketchInfo> sketches, bool sketchByFile, string outputFile)
{
	cerr << "output the result into: " << outputFile << endl;
	FILE *fp = fopen(outputFile.c_str(), "w");
	
	if(sketchByFile)
	{
		//sort the genomes within a cluster by degression of genome size.
		vector< vector<ClusterInfo> > cluster;
		for(int i = 0; i < clusterOrigin.size(); i++)
		{
			vector<ClusterInfo> tmpInfo;
			for(int j = 0; j < clusterOrigin[i].size(); j++)
			{
				int id = clusterOrigin[i][j];
				tmpInfo.push_back({id, sketches[id].totalSeqLength});
			}
			std::sort(tmpInfo.begin(), tmpInfo.end(), cmpSketchLength);
			cluster.push_back(tmpInfo);
			vector<ClusterInfo>().swap(tmpInfo);
		}

		for(int i = 0; i < cluster.size(); i++){
			//int tmpId = cluster[i][0].id;
			//fprintf(fp, "%s\n", sketches[tmpId].fileName.c_str());
			//continue;
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j].id;
				//printf("\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());
				//fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());
				fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileName.c_str(),  sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());
				//fprintf(fp, "\t%s\n", key.c_str());
				//fprintf(fp, "\t%20s\n", sketches[curId].fileName.c_str());
				//fprintf(fp, "\t%6d\t%6d\t%12dnt\t%30s%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileName.c_str(), sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());

			}
			//printf("\n");
			fprintf(fp, "\n");
		}
	}//end sketchByFile

	else//sketch by sequence
	{
		//sort the sequences within a cluster by degression of sequence length.
		vector< vector<ClusterInfo> > cluster;
		for(int i = 0; i < clusterOrigin.size(); i++)
		{
			vector<ClusterInfo> tmpInfo;
			for(int j = 0; j < clusterOrigin[i].size(); j++)
			{
				int id = clusterOrigin[i][j];
				tmpInfo.push_back({id, (uint64_t)sketches[id].seqInfo.length});
			}
			std::sort(tmpInfo.begin(), tmpInfo.end(), cmpSketchLength);
			cluster.push_back(tmpInfo);
			vector<ClusterInfo>().swap(tmpInfo);
		}

		for(int i = 0; i < cluster.size(); i++){
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j].id;		
				//printf("\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].seqInfo.length, sketches[curId].seqInfo.name.c_str(), sketches[curId].seqInfo.comment.c_str());
				fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].seqInfo.length, sketches[curId].seqInfo.name.c_str(), sketches[curId].seqInfo.comment.c_str());
			}
			//printf("\n");
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
			ofile << curFileSeqs.size() << endl;//get the current file sequence number
			ofile << sketches[i].fileName << endl;
			for(int j = 0; j < curFileSeqs.size(); j++)
			{
				/*	each line is name, strand, length, comment of the sequence.
				 *	place comment in a new line because there are space in the commment(the comment is not a single word.
				 *	When use the stringstream, it is sensitive to the space.
				 */												
				ofile << curFileSeqs[j].name << ' ' << ' ' << curFileSeqs[j].strand << ' ' << curFileSeqs[j].length << endl;
				ofile << curFileSeqs[j].comment << endl;
			}
		}
	}
	else//sketchBySequence
	{
		for(int i = 0; i < sketches.size(); i++)
		{
			SequenceInfo curSeq = sketches[i].seqInfo;
			ofile << curSeq.name << ' ' << ' ' << curSeq.strand << ' ' << curSeq.length << endl;
			ofile << curSeq.comment << endl;
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
	//	if(sketchByFile){
	//		//printf("{%6d,\t%6d,\t\t%lf }\t%60s\t%60s\n", mst[i].preNode, mst[i].sufNode, mst[i].dist, sketches[mst[i].preNode].fileName.c_str(), sketches[mst[i].sufNode].fileName.c_str());
	//		printf("{%6d,\t%6d,\t\t%lf }\t%60s\t%60s\n", sketches[mst[i].preNode].totalSeqLength, sketches[mst[i].sufNode].totalSeqLength, mst[i].dist, sketches[mst[i].preNode].fileName.c_str(), sketches[mst[i].sufNode].fileName.c_str());
	//	}
	//	else{
	//		printf("{%6d,\t%6d,\t\t%lf }\t%20s\t%20s\n", mst[i].preNode, mst[i].sufNode, mst[i].dist, sketches[mst[i].preNode].seqInfo.name.c_str(), sketches[mst[i].sufNode].seqInfo.name.c_str());
	//	}
		ofile1 << mst[i].preNode << ' ' << mst[i].sufNode << ' ' << mst[i].dist << endl;
	}
	cout << endl;
	ofile1.close();

}
