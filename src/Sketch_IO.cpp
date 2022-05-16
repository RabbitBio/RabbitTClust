#ifdef GREEDY_CLUST
#include "Sketch_IO.h"
#include <fstream>
#include <sstream>
#include "greedy.h"
#include "MST_IO.h"


/* @brief											Saving the sketches including hash values, genome informations into output files.
 * @details										To save the genome informations into outputGenomeInfo. To save the sketches informations 
 * 														and hash values into outputSketchInfo.
 * 														The first line in outputGenomeInfo is to determine genome input as a list file(sketchByFile)
 * 														or a single file. The remain lines are specific genome infos including fileName, genomeName, 
 * 														genomeComment, and totalLength, each line a genome.
 * 														The first four lines in outputSketchInfo is to determine iscontainment, containCompress(sketchSize),
 * 														kmerSize and sketchFunc. The remain lines are specific sketches hash values, each line a genome sketch.
 * 														
 * @param[in]	sketches				sketches array need to save
 * @param[in] inputFile				input file name for the prefix name of output file
 * @param[in] sketchFunc			Sketch function, including MinHash and KSSD, used for output file name
 * @param[in] inContainment		whether is used for duplication detection
 * @param[in]	containCompress	the dimension reduction for containment sketches
 * @param[in] sketchByFile		whether the input genomes are list of files or a single file
 * @param[in] sketchSize			the sketch size in the sketches
 * @param[in]	kmerSize				the kmer size in the sketches
 */
void saveSketches(vector<SketchInfo> sketches, string inputFile, string sketchFunc, bool isContainment, int containCompress, bool sketchByFile, int sketchSize, int kmerSize)
{
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
			//ofile << curSeq.comment << endl;
		}
	}
	ofile.close();

	cerr << "save the SketchInfo into: " << folderPath << '/' << prefixName+sketchFunc+"SketchInfo" << endl;
	ofstream ofile1;
	ofile1.open(folderPath + '/' + prefixName + sketchFunc+"SketchInfo");
	if(isContainment){
		ofile1 << '1' << endl;
		ofile1 << containCompress << endl;
	}
	else{
		ofile1 << '0' << endl;
		ofile1 << sketchSize << endl;
	}
	
	ofile1 << kmerSize << endl;
	ofile1 << sketchFunc << endl;
	
	if(sketchFunc == "MinHash"){
		for(int i = 0; i < sketches.size(); i++){
			vector<uint64_t> hashArr = sketches[i].minHash->storeMinHashes();
			for(int j = 0; j < hashArr.size(); j++)
				ofile1 << hashArr[j] << '\t';
			ofile1 << endl;
		}
	}

	ofile1.close();

}


/* @brief										Get the sketches informations and hash values from the input files, corresponding to the saveSketches function.
 * @details									Thsi function is in corresponding to the function saveSketches;
 * @param[in] inputFile0		Genome information file including fileName, seqName, seqComment, seqStrand and totalSeqLength. 
 * 													The first line in inputFile0 is to determine sketchByFile or sketchBySequences(input as a fileList or a single genome file).
 * @param[in] inputFile1 		Sketch information file including isContainment, containCompress(sketchSize), kmerSize, sketchFunc and hashValues. 
 * 													The first four lines determine the isContainment, containCompress(sketchSize), kmerSize and sketch function. 
 * 													The remaining lines are sorted hash values one genome per line.
 * @param[in] threshold			Distance threshold for clustering.
 * @param[in] threads				Thread number of clustering.
 */
void Sketch2Clust(string inputFile0, string inputFile1, string outputFile, double threshold, int threads)
{
#ifdef Timer
	double t0 = get_sec();
#endif
	fstream fs0(inputFile0);//Genome Info
	fstream fs1(inputFile1);//Sketch Info
	if(!fs0){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile0.c_str());
		printUsage();
		exit(1);
	}
	if(!fs1){
		fprintf(stderr, "error open the inputFile: %s\n", inputFile1.c_str());
		printUsage();
		exit(1);
	}

	string line;
	getline(fs0, line);
	bool sketchByFile = stoi(line);//1 or 0
	getline(fs1, line);
	bool isContainment = stoi(line);//1 or 0
	getline(fs1, line);
	int containCompress(10000), sketchSize(1000);
	if(isContainment)
		containCompress = stoi(line);
	else
		sketchSize = stoi(line);
	getline(fs1, line);
	int kmerSize = stoi(line);
	getline(fs1, line);
	string sketchFunc = line;

	vector<SketchInfo> sketches;
	Sketch::MinHash * mh1;
	int sketchId = 0;

	while(getline(fs0, line)){
		SketchInfo tmpSketchInfo;

		stringstream ss;
		ss << line;
		string fileName, seqName, seqComment, tmpComment;
		int seqStrand, seqLength(0);
		uint64_t totalLength;
		if(sketchByFile){
			ss >> fileName >> seqName >> seqStrand >> totalLength; 
			while(ss >> tmpComment){
				seqComment += tmpComment + ' ';
			}
			seqComment = seqComment.substr(0, seqComment.length()-1);
			SequenceInfo tmpSeq{seqName, seqComment, seqStrand, seqLength};
			Vec_SeqInfo curFileSeqs;
			curFileSeqs.push_back(tmpSeq);
			tmpSketchInfo.fileName = fileName;
			tmpSketchInfo.totalSeqLength = totalLength;
			tmpSketchInfo.fileSeqs = curFileSeqs;
		}
		else{
			ss >> seqName >> seqStrand >> seqLength;
			while(ss >> tmpComment){
				seqComment += tmpComment + ' ';
			}
			seqComment = seqComment.substr(0, seqComment.length()-1);
			//getline(fs0, seqComment);//comment is a single line;
			SequenceInfo curSeq{seqName, seqComment, seqStrand, seqLength};
			tmpSketchInfo.seqInfo = curSeq;
		}

		getline(fs1, line);
		stringstream ss1;
		ss1 << line;
		vector<uint64_t> hashArr;
		uint64_t hashValue;
		while(ss1 >> hashValue){
			hashArr.push_back(hashValue);
		}
		
		if(sketchFunc == "MinHash"){
			if(isContainment){
				mh1 = new Sketch::MinHash(kmerSize, containCompress);
				tmpSketchInfo.isContainment = true;
			}
			else{
				mh1 = new Sketch::MinHash(kmerSize, sketchSize);
			}
			mh1->loadMinHashes(hashArr);
			tmpSketchInfo.minHash = mh1;
		}
		tmpSketchInfo.id = sketchId++;
		sketches.push_back(tmpSketchInfo);
	}//end while
#ifdef Timer
	double t1 = get_sec();
	cerr << "========time of load genome Infos and sketch Infos is: " << t1 - t0 << endl;
#endif

	vector<vector<int> > cluster = greedyCluster(sketches, sketchFunc, threshold, threads);
	printResult(cluster, sketches, sketchByFile, outputFile);
#ifdef Timer
	double t2 = get_sec();
	cerr << "========time of greedy incremental cluster is: " << t2 - t1 << endl;
#endif
	
}



#endif
