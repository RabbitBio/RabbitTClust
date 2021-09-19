#include "MST_IO.h"
//#include "parameter.h"
using namespace std;

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

void printResult(vector< vector<int> > cluster, vector<SketchInfo> sketches, bool sketchByFile, string outputFile)
{
	cerr << "output the result into: " << outputFile << endl;
	FILE *fp = fopen(outputFile.c_str(), "w");
	
	if(sketchByFile)
	{
		for(int i = 0; i < cluster.size(); i++){
			//printf("the cluster %d is: \n", i);
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j];
				//printf("\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());
				fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].totalSeqLength, sketches[curId].fileSeqs[0].name.c_str(), sketches[curId].fileSeqs[0].comment.c_str());

			}
			//printf("\n");
			fprintf(fp, "\n");
		}
	}//end sketchByFile

	else//sketch by sequence
	{
		for(int i = 0; i < cluster.size(); i++){
			//printf("the cluster %d is: \n", i);
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j];		
				//printf("\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].seqInfo.length, sketches[curId].seqInfo.name.c_str(), sketches[curId].seqInfo.comment.c_str());
				fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n", j, curId, sketches[curId].seqInfo.length, sketches[curId].seqInfo.name.c_str(), sketches[curId].seqInfo.comment.c_str());
			}
			//printf("\n");
			fprintf(fp, "\n");
		}
	}//end sketchBySequence
	fclose(fp);

}

void saveMST(string inputFile, string sketchFunc, vector<SketchInfo> sketches, vector<EdgeInfo> mst, bool sketchByFile, int sketchSize, int kmerSize)
{
	//save the matching of graph id and genomeInfo 
	cerr << "save the genomeInfo into: " << inputFile+sketchFunc+"GenomeInfo" << endl;
	ofstream ofile;
	ofile.open(inputFile+sketchFunc+"GenomeInfo");
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
	cerr << "save the MSTInfo into: " << inputFile+sketchFunc+"MSTInfo" << endl;
	ofstream ofile1;
	ofile1.open(inputFile+sketchFunc+"MSTInfo");

	if(sketchByFile)
		ofile1 << "sketch by File! " << endl;
	else 
		ofile1 << "sketch by Sequence!" << endl;

	ofile1 << "the sketch function is: " << sketchFunc << endl;
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
