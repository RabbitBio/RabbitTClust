#include "MST_IO.h"
#include "Sketch_IO.h"
#include <ctime>
using namespace std;


inline bool cmpSketchLength(ClusterInfo c1, ClusterInfo c2){
	return c1.length > c2.length;
}

void loadDense(int** &denseArr, string folderPath, int& denseSpan, int& genome_number){
	string file_dense = folderPath + '/' + "mst.dense";
	FILE* fp_dense = fopen(file_dense.c_str(), "r");
	if(!fp_dense){
		cerr << "ERROR: saveDense(), cannot open the file: " << file_dense;
		exit(1);
	}
	fread(&genome_number, sizeof(int), 1, fp_dense);
	fread(&denseSpan, sizeof(int), 1, fp_dense);
	denseArr = new int*[denseSpan];
	for(int i = 0; i < denseSpan; i++){
		denseArr[i] = new int[genome_number];
		fread(denseArr[i], sizeof(int), genome_number, fp_dense);
	}
	fclose(fp_dense);
	cerr << "-----read the dense file from: " << file_dense << endl;
}

void loadANI(string folderPath, uint64_t* &aniArr, int sketch_func_id){
	if(sketch_func_id != 0 && sketch_func_id != 1){
		cerr << "ERROR: saveANI(), save ANI can only support MinHash and KSSD functions" << endl;
		return;
	}
	string file_ani = folderPath + '/' + "mst.ani";
	FILE* fp_ani = fopen(file_ani.c_str(), "r");
	if(!fp_ani){
		cerr << "ERROR: saveANI(), cannot open file: " << file_ani << endl;
		exit(1);
	}
	aniArr = new uint64_t[101];
	fread(aniArr, sizeof(uint64_t), 101, fp_ani);
	fclose(fp_ani);
	cerr << "-----read the ani file from: " << file_ani << endl;
}

void loadMST(string folderPath, vector<EdgeInfo>& mst)
{
	//load the mst edge 
	string file_mst = folderPath + '/' + "edge.mst";
	FILE* fp_mst = fopen(file_mst.c_str(), "r");
	if(!fp_mst){
		cerr << "ERROR: loadMST(), cannot open the file: " <<  file_mst << endl;
		exit(1);
	}
	size_t mst_size;
	fread(&mst_size, sizeof(size_t), 1, fp_mst);
	int preNode, sufNode;
	double dist;
	for(size_t i = 0; i < mst_size; i++){
		fread(&preNode, sizeof(int), 1, fp_mst);
		fread(&sufNode, sizeof(int), 1, fp_mst);
		fread(&dist, sizeof(double), 1, fp_mst);
		EdgeInfo tmpEdge{preNode, sufNode, dist};
		mst.push_back(tmpEdge);
		//cout << preNode << '\t' << sufNode << '\t' << dist << endl;
	}
	fclose(fp_mst);
	cerr << "-----read the mst file from " << file_mst << endl;
}

void printResult(vector<vector<int>>& cluster, vector<SketchInfo>& sketches, bool sketchByFile, string outputFile)
{
	//cerr << "output the result into: " << outputFile << endl;
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

void saveMST(vector<SketchInfo>& sketches, vector<EdgeInfo>& mst, string folderPath, bool sketchByFile){
	save_genome_info(sketches, folderPath, "mst", sketchByFile);
	string file_mst = folderPath + '/' + "edge.mst";
	FILE* fp_mst = fopen(file_mst.c_str(), "w+");
	if(!fp_mst){
		cerr << "ERROR: saveMST(), cannot open the file: " <<  file_mst << endl;
		exit(1);
	}
	size_t mst_size = mst.size();
	fwrite(&mst_size, sizeof(size_t), 1, fp_mst);
	for(size_t i = 0; i < mst.size(); i++){
		fwrite(&mst[i].preNode, sizeof(int), 1, fp_mst);
		fwrite(&mst[i].sufNode, sizeof(int), 1, fp_mst);
		fwrite(&mst[i].dist, sizeof(double), 1, fp_mst);
	}
	fclose(fp_mst);
	cerr << "-----save the mst into: " << folderPath << endl;
}

void saveDense(string folderPath, int** denseArr, int denseSpan, int genome_number){
	string file_dense = folderPath + '/' + "mst.dense";
	FILE* fp_dense = fopen(file_dense.c_str(), "w+");
	if(!fp_dense){
		cerr << "ERROR: saveDense(), cannot open the file: " << file_dense;
		exit(1);
	}
	fwrite(&genome_number, sizeof(int), 1, fp_dense);
	fwrite(&denseSpan, sizeof(int), 1, fp_dense);
	for(int i = 0; i < denseSpan; i++){
		fwrite(denseArr[i], sizeof(int), genome_number, fp_dense);
	}
	fclose(fp_dense);
	cerr << "-----save the dense file into: " << folderPath << endl;
}

void saveANI(string folderPath, uint64_t* aniArr, int sketch_func_id){
	
	if(sketch_func_id != 0 && sketch_func_id != 1){
		cerr << "ERROR: saveANI(), save ANI can only support MinHash and KSSD functions" << endl;
		return;
	}
	string file_ani = folderPath + '/' + "mst.ani";
	FILE* fp_ani = fopen(file_ani.c_str(), "w+");
	if(!fp_ani){
		cerr << "ERROR: saveANI(), cannot open file: " << file_ani << endl;
		exit(1);
	}
	fwrite(aniArr, sizeof(uint64_t), 101, fp_ani);
	fclose(fp_ani);
	cerr << "-----save the ani file into: " << file_ani << endl;
}

void print_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	string res_newick_tree = get_newick_tree(sketches, mst, sketch_by_file);
	FILE* fp_tree = fopen(output.c_str(), "w");
	if(!fp_tree){
		cerr << "ERROR: print_newick_tree(), cannot write file: " << output << endl;
		exit(1);
	}
	fprintf(fp_tree, "%s\n", res_newick_tree.c_str());
	fclose(fp_tree);
}



