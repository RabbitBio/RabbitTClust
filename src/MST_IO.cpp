#include "MST_IO.h"
#include "Sketch_IO.h"
#include "TreeFormat.h"
#include <ctime>
#include <cstdlib>
#include <exception>
#include <unordered_map>
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

void printKssdResult(vector<vector<int>>& cluster, vector<KssdSketchInfo>& sketches, bool sketchByFile, string outputFile, double threshold)
{
	//cerr << "output the result into: " << outputFile << endl;
	FILE *fp = fopen(outputFile.c_str(), "w");
	if(!fp){
		cerr << "Error in printKssdResult(), cannot open file: " << outputFile << endl;
		exit(1);
	}
	
	// Write threshold information at the beginning of the file
	if(threshold >= 0.0){
		fprintf(fp, "# Clustering threshold: %.6f\n", threshold);
		fprintf(fp, "# Total clusters: %zu\n", cluster.size());
		fprintf(fp, "#\n");
	}
	
	if(sketchByFile)
	{
		for(int i = 0; i < cluster.size(); i++){
			fprintf(fp, "the cluster %d is: \n", i);
			for(int j = 0; j < cluster[i].size(); j++)
			{
				int curId = cluster[i][j];
				if(curId < 0 || curId >= (int)sketches.size()){
					cerr << "ERROR: printKssdResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
					continue;
				}
				const char* seqName = "N/A";
				const char* seqComment = "N/A";
				if(!sketches[curId].fileSeqs.empty()){
					seqName = sketches[curId].fileSeqs[0].name.c_str();
					seqComment = sketches[curId].fileSeqs[0].comment.c_str();
				}
				fprintf(fp, "\t%5d\t%6d\t%12dnt\t%20s\t%20s\t%s\n",
					j, curId, sketches[curId].totalSeqLength,
					sketches[curId].fileName.c_str(), seqName, seqComment);
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
				if(curId < 0 || curId >= (int)sketches.size()){
					cerr << "ERROR: printKssdResult(), invalid index " << curId << " (sketches.size=" << sketches.size() << "), skipping." << endl;
					continue;
				}
				fprintf(fp, "\t%6d\t%6d\t%12dnt\t%20s\t%s\n",
					j, curId, sketches[curId].seqInfo.length,
					sketches[curId].seqInfo.name.c_str(),
					sketches[curId].seqInfo.comment.c_str());
			}
			fprintf(fp, "\n");
		}
	}//end sketchBySequence
	fclose(fp);

}

void printResult(vector<vector<int>>& cluster, vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, double threshold)
{
	//cerr << "output the result into: " << outputFile << endl;
	FILE *fp = fopen(outputFile.c_str(), "w");
	if(!fp){
		cerr << "Error in printResult(), cannot open file: " << outputFile << endl;
		exit(1);
	}
	
	// Write threshold information at the beginning of the file
	if(threshold >= 0.0){
		fprintf(fp, "# Clustering threshold: %.6f\n", threshold);
		fprintf(fp, "# Total clusters: %zu\n", cluster.size());
		fprintf(fp, "#\n");
	}
	
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

void saveKssdMST(vector<KssdSketchInfo>& sketches, vector<EdgeInfo>& mst, string folderPath, bool sketchByFile){
	save_kssd_genome_info(sketches, folderPath, "mst", sketchByFile);
	string file_mst = folderPath + '/' + "edge.mst";
	FILE* fp_mst = fopen(file_mst.c_str(), "w+");
	if(!fp_mst){
		cerr << "ERROR: saveKsdMST(), cannot open the file: " <<  file_mst << endl;
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
	cerr << "-----save the kssd mst into: " << folderPath << endl;
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

[[noreturn]] static void tree_output_error(
		const char* function_name,
		const string& output,
		const string& detail){
	cerr << "ERROR: " << function_name << "(), cannot write a valid tree to: "
	     << output << endl;
	cerr << "       " << detail << endl;
	exit(EXIT_FAILURE);
}

static void write_newick_file(
		const string& newick_tree,
		const string& output,
		const char* function_name){
	FILE* fp_tree = fopen(output.c_str(), "w");
	if(!fp_tree){
		tree_output_error(function_name, output, "cannot open the output file");
	}
	fprintf(fp_tree, "%s\n", newick_tree.c_str());
	fclose(fp_tree);
}

void print_kssd_newick_tree(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	try {
		write_newick_file(
			get_kssd_newick_tree(sketches, mst, sketch_by_file),
			output,
			"print_kssd_newick_tree");
	} catch(const exception& error) {
		tree_output_error("print_kssd_newick_tree", output, error.what());
	}
}

void print_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	try {
		write_newick_file(
			get_newick_tree(sketches, mst, sketch_by_file),
			output,
			"print_newick_tree");
	} catch(const exception& error) {
		tree_output_error("print_newick_tree", output, error.what());
	}
}

template <typename SketchType>
static vector<string> get_tree_labels(
		const vector<SketchType>& sketches,
		bool sketch_by_file){
	vector<string> labels;
	labels.reserve(sketches.size());
	for(const SketchType& sketch : sketches){
		labels.push_back(sketch_by_file ? sketch.fileName : sketch.seqInfo.name);
	}
	return labels;
}

static string escape_tsv_field(const string& value){
	string escaped;
	escaped.reserve(value.size());
	for(char c : value){
		switch(c){
			case '\\': escaped += "\\\\"; break;
			case '\t': escaped += "\\t"; break;
			case '\n': escaped += "\\n"; break;
			case '\r': escaped += "\\r"; break;
			default: escaped += c; break;
		}
	}
	return escaped;
}

template <typename SketchType>
static void print_phylip_tree_impl(
		const vector<SketchType>& sketches,
		const vector<EdgeInfo>& mst,
		bool sketch_by_file,
		const string& output,
		const char* function_name){
	if(sketches.size() < 2){
		tree_output_error(
			function_name,
			output,
			"PHYLIP tree output requires at least two taxa");
	}

	const vector<string> original_labels =
		get_tree_labels(sketches, sketch_by_file);
	tree_format::PhylipLabelSet phylip_labels;
	string newick_tree;
	try {
		phylip_labels = tree_format::make_phylip_labels(original_labels);
		newick_tree = tree_format::build_mst_newick(phylip_labels.encoded, mst);
	} catch(const exception& error) {
		tree_output_error(function_name, output, error.what());
	}

	FILE* fp_tree = fopen(output.c_str(), "w");
	if(!fp_tree){
		tree_output_error(function_name, output, "cannot open the output file");
	}

	const string label_map_path = output + ".labels.tsv";
	FILE* fp_labels = fopen(label_map_path.c_str(), "w");
	if(!fp_labels){
		fclose(fp_tree);
		tree_output_error(
			function_name, label_map_path, "cannot open the taxon-label mapping file");
	}

	// PHYLIP accepts the historical tree-count prefix.  There is one connected
	// MST, followed by the same topology and branch lengths in Newick notation.
	fprintf(fp_tree, "1\n%s\n", newick_tree.c_str());
	fclose(fp_tree);

	fprintf(fp_labels, "phylip_label\toriginal_label\n");
	for(size_t index = 0; index < original_labels.size(); ++index){
		const string escaped_original = escape_tsv_field(original_labels[index]);
		fprintf(fp_labels, "%s\t%s\n",
			phylip_labels.encoded[index].c_str(), escaped_original.c_str());
	}
	fclose(fp_labels);
	cerr << "-----write the PHYLIP taxon-label mapping into: "
	     << label_map_path << endl;
}

// PHYLIP-compatible tree file: first line = number of trees (1),
// second line = the Newick-formatted MST.
void print_phylip_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	print_phylip_tree_impl(
		sketches, mst, sketch_by_file, output, "print_phylip_tree");
}

void print_kssd_phylip_tree(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	print_phylip_tree_impl(
		sketches, mst, sketch_by_file, output, "print_kssd_phylip_tree");
}

// NEXUS format: TAXA block with TAXLABELS + TREES block with one Newick tree
static void write_nexus_tree(const string& newick_tree, const vector<string>& taxlabels, const string& output){
	FILE* fp = fopen(output.c_str(), "w");
	if(!fp){
		cerr << "ERROR: print_nexus_tree(), cannot write file: " << output << endl;
		exit(1);
	}
	fprintf(fp, "#NEXUS\n");
	fprintf(fp, "BEGIN TAXA;\n");
	fprintf(fp, "  DIMENSIONS NTAX=%zu;\n", taxlabels.size());
	fprintf(fp, "  TAXLABELS");
	for(size_t i = 0; i < taxlabels.size(); ++i){
		// NEXUS TAXLABELS: escape single quotes in names by doubling them
		string lab = taxlabels[i];
		size_t pos = 0;
		while((pos = lab.find('\'', pos)) != string::npos){
			lab.insert(pos, "'");
			pos += 2;
		}
		fprintf(fp, " '%s'", lab.c_str());
	}
	fprintf(fp, ";\n");
	fprintf(fp, "END;\n");
	fprintf(fp, "BEGIN TREES;\n");
	fprintf(fp, "  TREE tree_1 = [&R] %s\n", newick_tree.c_str());
	fprintf(fp, "END;\n");
	fclose(fp);
}

void print_nexus_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	try {
		string res_newick_tree = get_newick_tree(sketches, mst, sketch_by_file);
		vector<string> taxlabels;
		taxlabels.reserve(sketches.size());
		for(size_t i = 0; i < sketches.size(); ++i){
			taxlabels.push_back(sketch_by_file ? sketches[i].fileName : sketches[i].seqInfo.name);
		}
		write_nexus_tree(res_newick_tree, taxlabels, output);
	} catch(const exception& error) {
		tree_output_error("print_nexus_tree", output, error.what());
	}
}

void print_kssd_nexus_tree(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file, string output){
	try {
		string res_newick_tree = get_kssd_newick_tree(sketches, mst, sketch_by_file);
		vector<string> taxlabels;
		taxlabels.reserve(sketches.size());
		for(size_t i = 0; i < sketches.size(); ++i){
			taxlabels.push_back(sketch_by_file ? sketches[i].fileName : sketches[i].seqInfo.name);
		}
		write_nexus_tree(res_newick_tree, taxlabels, output);
	} catch(const exception& error) {
		tree_output_error("print_kssd_nexus_tree", output, error.what());
	}
}

void print_kssd_linkage_matrix(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, string output){
	int N = sketches.size();
	vector<LinkageRow> Z = get_linkage_from_mst(N, mst);

	FILE* fp = fopen(output.c_str(), "w");
	if(!fp){
		cerr << "ERROR: print_kssd_linkage_matrix(), cannot write file: " << output << endl;
		exit(1);
	}
	// 每行：cluster1  cluster2  dist  size
	for(const auto& row : Z){
		fprintf(fp, "%d\t%d\t%.6f\t%d\n", row.c1, row.c2, row.dist, row.size);
	}
	fclose(fp);
}

void print_linkage_matrix(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, string output){
	int N = sketches.size();
	vector<LinkageRow> Z = get_linkage_from_mst(N, mst);

	FILE* fp = fopen(output.c_str(), "w");
	if(!fp){
		cerr << "ERROR: print_linkage_matrix(), cannot write file: " << output << endl;
		exit(1);
	}
	for(const auto& row : Z){
		fprintf(fp, "%d\t%d\t%.6f\t%d\n", row.c1, row.c2, row.dist, row.size);
	}
	fclose(fp);
}
