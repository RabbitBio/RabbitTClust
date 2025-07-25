#include "sub_command.h"
#include <assert.h>
#include <mpi.h>
#include <sys/stat.h>
using namespace std;


size_t get_file_size(string file){
  struct stat file_stat;
  if(stat(file.c_str(), &file_stat) == -1){
    cerr << "ERROR: get_file_size(), failed to get file status of: " << file << endl;
    exit(1);
  }
  size_t file_size = file_stat.st_size;
  return file_size;
}

void build_message(char* &buffer, size_t& file_size, string file_name){
  file_size = get_file_size(file_name);
  buffer = new char[file_size];
  FILE * fp = fopen(file_name.c_str(), "r");
  assert(fp != NULL);
  size_t read_length = fread(buffer, sizeof(char), file_size, fp);
  assert(read_length == file_size);
}

void format_sketches(char* info_buffer, size_t info_size, char* hash_buffer, size_t hash_size, string folder_path, vector<KssdSketchInfo>& sketches, bool sketch_by_file, int threads){
  string info_file = folder_path + '/' + "kssd.info.sketch";
  string hash_file = folder_path + '/' + "kssd.hash.sketch";
  FILE* fp_info = fopen(info_file.c_str(), "w");
  FILE* fp_hash = fopen(hash_file.c_str(), "w");
  fwrite(info_buffer, sizeof(char), info_size, fp_info);
  fwrite(hash_buffer, sizeof(char), hash_size, fp_hash);
  fclose(fp_info);
  fclose(fp_hash);
  int sketch_func_id;
  KssdParameters info;
  bool cur_sketch_by_file = loadKssdSketches(folder_path, threads, sketches, info);
  assert(cur_sketch_by_file == sketch_by_file);
}


void format_sketches_index(char* info_buffer, size_t info_size, char* hash_buffer, size_t hash_size, char* index_buffer, size_t index_size, char* dict_buffer, size_t dict_size, string folder_path, vector<KssdSketchInfo>& sketches, bool sketch_by_file, int threads){
  string info_file = folder_path + '/' + "kssd.info.sketch";
  string hash_file = folder_path + '/' + "kssd.hash.sketch";
  string index_file = folder_path + '/' + "kssd.sketch.index";
  string dict_file = folder_path + '/' + "kssd.sketch.dict";
  FILE* fp_info = fopen(info_file.c_str(), "w");
  FILE* fp_hash = fopen(hash_file.c_str(), "w");
  FILE* fp_index = fopen(index_file.c_str(), "w");
  FILE* fp_dict = fopen(dict_file.c_str(), "w");
  fwrite(info_buffer, sizeof(char), info_size, fp_info);
  fwrite(hash_buffer, sizeof(char), hash_size, fp_hash);
  fwrite(index_buffer, sizeof(char), index_size, fp_index);
  fwrite(dict_buffer, sizeof(char), dict_size, fp_dict);
  fclose(fp_info);
  fclose(fp_hash);
  fclose(fp_index);
  fclose(fp_dict);
  int sketch_func_id;
  KssdParameters info;
  bool cur_sketch_by_file = loadKssdSketches(folder_path, threads, sketches, info);
  assert(cur_sketch_by_file == sketch_by_file);
}
void load_MST(string folderPath, vector<EdgeInfo>& mst)
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

void format_MST(int my_rank, char* edge_buffer, size_t edge_size, string folder_path, vector<EdgeInfo>& mst){
  string edge_file = folder_path + '/' + "edge.mst";
  FILE *fp_edge = fopen(edge_file.c_str(), "w");
  fwrite(edge_buffer, sizeof(char), edge_size, fp_edge);
  fclose(fp_edge);
  load_MST(folder_path, mst);
}

#ifdef GREEDY_CLUST








void append_clust_greedy(string folder_path, string input_file, string output_file, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
  int sketch_func_id_0; 
  vector<SketchInfo> pre_sketches; 
  bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0); 
  if(pre_sketch_by_file != sketch_by_file){
    cerr << "ERROR: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
    //cerr << "the output cluster file may not have the genome file name" << endl;
    exit(1);
  }
  int sketch_func_id_1, kmer_size, contain_compress, sketch_size, half_k, half_subk, drlevel;
  bool is_containment;
  read_sketch_parameters(folder_path, sketch_func_id_1, kmer_size, is_containment, contain_compress, sketch_size, half_k, half_subk, drlevel);
  assert(sketch_func_id_0 == sketch_func_id_1);
  cerr << "-----use the same sketch parameters with pre-generated sketches" << endl;
  if(sketch_func_id_0 == 0){
    cerr << "---the kmer size is: " << kmer_size << endl;
    if(is_containment)
      cerr << "---use the AAF distance (variable-sketch-size), the sketch size is in proportion with 1/" << contain_compress << endl;
    else 
      cerr << "---use the Mash distance (fixed-sketch-size), the sketch size is: " << sketch_size << endl;
  }
  else if(sketch_func_id_0 == 1){
    cerr << "---use the KSSD sketches" << endl;
    cerr << "---the half_k is: " << half_k << endl;
    cerr << "---the half_subk is: " << half_subk << endl;
    cerr << "---the drlevel is: " << drlevel << endl;
  }
  cerr << "---the thread number is: " << threads << endl;
  cerr << "---the threshold is: " << threshold << endl;
  string sketch_func;
  if(sketch_func_id_0 == 0) 
    sketch_func = "MinHash";
  else if(sketch_func_id_0 == 1)
    sketch_func = "KSSD";
  vector<SketchInfo> append_sketches;
  string append_folder_path;
  compute_sketches(append_sketches, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, false, threads);
  vector<SketchInfo> final_sketches;
  final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
  final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
  if(sketch_by_file)
    sort(final_sketches.begin(), final_sketches.end(), cmpGenomeSize);
  else
    sort(final_sketches.begin(), final_sketches.end(), cmpSeqSize);
  vector<SketchInfo>().swap(pre_sketches);
  vector<SketchInfo>().swap(append_sketches);
  string new_folder_path = currentDataTime();
  if(!no_save){
    string command = "mkdir -p " + new_folder_path;
    system(command.c_str());
    saveSketches(final_sketches, new_folder_path, sketch_by_file, sketch_func, is_containment, contain_compress, sketch_size, kmer_size);
  }
  vector<vector<int>> cluster;
  cluster = greedyCluster(final_sketches, sketch_func_id_0, threshold, threads);
  printResult(cluster, final_sketches, sketch_by_file, output_file);
  cerr << "-----write the cluster result into: " << output_file << endl;
  cerr << "-----the cluster number of " << output_file << " is: " << cluster.size() << endl;
}
#endif

#ifndef GREEDY_CLUST
void append_clust_mst_fast(string folder_path, string input_file, string output_file, bool is_newick_tree, bool no_dense, bool sketch_by_file, bool isContainment, int min_len, bool no_save, double threshold, int threads){
  bool isSave = !no_save;
  int sketch_func_id_0; 
  vector<KssdSketchInfo> pre_sketches; 
  KssdParameters pre_info;
  bool pre_sketch_by_file = loadKssdSketches(folder_path, threads, pre_sketches, pre_info); 
  if(pre_sketch_by_file != sketch_by_file){
    cerr << "Warning: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
    cerr << "the output cluster file may not have the genome file name" << endl;
  }
  vector<EdgeInfo> pre_mst;
  load_MST(folder_path, pre_mst);
  int kmer_size = pre_info.half_k * 2;
  int drlevel = pre_info.drlevel;

  cerr << "---the thread number is: " << threads << endl;
  cerr << "---the threshold is: " << threshold << endl;

  vector<KssdSketchInfo> append_sketches;
  KssdParameters append_info;
  string append_folder_path;
  compute_kssd_sketches(append_sketches, append_info, isSave, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, drlevel, threads);

  vector<KssdSketchInfo> final_sketches;
  int pre_sketch_size = pre_sketches.size();
  final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
  final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
  vector<KssdSketchInfo>().swap(pre_sketches);
  vector<KssdSketchInfo>().swap(append_sketches);
  string new_folder_path = append_folder_path;
  if(!no_save){
    string command = "mkdir -p " + new_folder_path;
    system(command.c_str());
    saveKssdSketches(final_sketches, pre_info, new_folder_path, sketch_by_file);
  }
  transSketches(final_sketches, append_info, new_folder_path, threads);

  int ** pre_dense_arr;
  uint64_t* pre_ani_arr;
  int pre_dense_span;
  int pre_genome_number;
  if(!no_dense){
    loadDense(pre_dense_arr, folder_path, pre_dense_span, pre_genome_number);
  }

  int ** dense_arr;
  int dense_span = DENSE_SPAN;
  uint64_t* ani_arr;
  vector<EdgeInfo> append_mst = compute_kssd_mst(final_sketches,0,0, append_info, new_folder_path, no_dense, isContainment, threads, dense_arr, dense_span, ani_arr, threshold);
  vector<EdgeInfo> final_graph;
  final_graph.insert(final_graph.end(), pre_mst.begin(), pre_mst.end());
  final_graph.insert(final_graph.end(), append_mst.begin(), append_mst.end());
  vector<EdgeInfo>().swap(pre_mst);
  vector<EdgeInfo>().swap(append_mst);
  sort(final_graph.begin(), final_graph.end(), cmpEdge);
  vector<EdgeInfo> final_mst = kruskalAlgorithm(final_graph, final_sketches.size());
  vector<EdgeInfo>().swap(final_graph);
  if(is_newick_tree){
    string output_newick_file = output_file + ".newick.tree";
    print_kssd_newick_tree(final_sketches, final_mst, pre_sketch_by_file, output_newick_file);
    cerr << "-----write the newick tree into: " << output_newick_file << endl;
  }

  vector<EdgeInfo> forest = generateForest(final_mst, threshold);
  vector<vector<int>> tmpClust = generateClusterWithBfs(forest, final_sketches.size());
  printKssdResult(tmpClust, final_sketches, pre_sketch_by_file, output_file);
  cerr << "-----write the cluster result into: " << output_file << endl;
  cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
  if(!no_save){
    saveKssdMST(final_sketches, final_mst, new_folder_path, sketch_by_file);
  }

  if(!no_dense){
    loadANI(folder_path, pre_ani_arr, sketch_func_id_0);
    for(int i = 0; i < 101; i++)
      ani_arr[i] += pre_ani_arr[i];
    for(int i = 0; i < pre_dense_span; i++){
      for(int j = 0; j < pre_genome_number; j++){
        dense_arr[i][j] += pre_dense_arr[i][j];
      }
    }

    if(!no_save){
      saveANI(new_folder_path, ani_arr, sketch_func_id_0);
      saveDense(new_folder_path, dense_arr, dense_span, final_sketches.size());
    }

    int alpha = 2;
    int denseIndex = threshold / 0.01;
    vector<int> totalNoiseArr;
    for(int i = 0; i < tmpClust.size(); i++){
      if(tmpClust[i].size() == 1) continue;
      vector<PairInt> curDenseArr;
      set<int> denseSet;
      for(int j = 0; j < tmpClust[i].size(); j++){
        int element = tmpClust[i][j];
        PairInt p(element, dense_arr[denseIndex][element]);
        denseSet.insert(dense_arr[denseIndex][element]);
        curDenseArr.push_back(p);
      }
      vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
      totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
    }
    cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
    forest = modifyForest(forest, totalNoiseArr, threads);
    vector<vector<int>> cluster = generateClusterWithBfs(forest, final_sketches.size());
    string outputFileNew = output_file + ".removeNoise";
    printKssdResult(cluster, final_sketches, pre_sketch_by_file, outputFileNew);
    cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
    cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
  }

}

void append_clust_mst(string folder_path, string input_file, string output_file, bool is_newick_tree, bool no_dense, bool sketch_by_file, int min_len, bool no_save, double threshold, int threads){
  int sketch_func_id_0; 
  vector<SketchInfo> pre_sketches; 
  bool pre_sketch_by_file = loadSketches(folder_path, threads, pre_sketches, sketch_func_id_0); 
  if(pre_sketch_by_file != sketch_by_file){
    cerr << "Warning: append_clust_mst(), the input format of append genomes and pre-sketched genome is not same (single input genome vs. genome list)" << endl;
    cerr << "the output cluster file may not have the genome file name" << endl;
  }
  vector<EdgeInfo> pre_mst;
  load_MST(folder_path, pre_mst);
  int sketch_func_id_1, kmer_size, contain_compress, sketch_size, half_k, half_subk, drlevel;
  bool is_containment;
  read_sketch_parameters(folder_path, sketch_func_id_1, kmer_size, is_containment, contain_compress, sketch_size, half_k, half_subk, drlevel);
  assert(sketch_func_id_0 == sketch_func_id_1);
  cerr << "-----use the same sketch parameters with pre-generated sketches" << endl;
  if(sketch_func_id_0 == 0){
    cerr << "---the kmer size is: " << kmer_size << endl;
    if(is_containment)
      cerr << "---use the AAF distance (variable-sketch-size), the sketch size is in proportion with 1/" << contain_compress << endl;
    else 
      cerr << "---use the Mash distance (fixed-sketch-size), the sketch size is: " << sketch_size << endl;
  }
  else if(sketch_func_id_0 == 1){
    cerr << "---use the KSSD sketches" << endl;
    cerr << "---the half_k is: " << half_k << endl;
    cerr << "---the half_subk is: " << half_subk << endl;
    cerr << "---the drlevel is: " << drlevel << endl;
  }
  cerr << "---the thread number is: " << threads << endl;
  cerr << "---the threshold is: " << threshold << endl;
  string sketch_func;
  if(sketch_func_id_0 == 0) 
    sketch_func = "MinHash";
  else if(sketch_func_id_0 == 1)
    sketch_func = "KSSD";
  vector<SketchInfo> append_sketches;
  string append_folder_path;
  compute_sketches(append_sketches, input_file, append_folder_path, sketch_by_file, min_len, kmer_size, sketch_size, sketch_func, is_containment, contain_compress, false, threads);

  vector<SketchInfo> final_sketches;
  int pre_sketch_size = pre_sketches.size();
  final_sketches.insert(final_sketches.end(), pre_sketches.begin(), pre_sketches.end());
  final_sketches.insert(final_sketches.end(), append_sketches.begin(), append_sketches.end());
  vector<SketchInfo>().swap(pre_sketches);
  vector<SketchInfo>().swap(append_sketches);
  string new_folder_path = currentDataTime();
  if(!no_save){
    string command = "mkdir -p " + new_folder_path;
    system(command.c_str());
    saveSketches(final_sketches, new_folder_path, sketch_by_file, sketch_func, is_containment, contain_compress, sketch_size, kmer_size);
  }

  int ** pre_dense_arr;
  uint64_t* pre_ani_arr;
  int pre_dense_span;
  int pre_genome_number;
  if(!no_dense){
    loadDense(pre_dense_arr, folder_path, pre_dense_span, pre_genome_number);
  }

  int ** dense_arr;
  int dense_span = DENSE_SPAN;
  uint64_t* ani_arr;
  vector<EdgeInfo> append_mst = modifyMST(final_sketches, pre_sketch_size, sketch_func_id_0, threads, no_dense, dense_arr, dense_span, ani_arr);
  vector<EdgeInfo> final_graph;
  final_graph.insert(final_graph.end(), pre_mst.begin(), pre_mst.end());
  final_graph.insert(final_graph.end(), append_mst.begin(), append_mst.end());
  vector<EdgeInfo>().swap(pre_mst);
  vector<EdgeInfo>().swap(append_mst);
  sort(final_graph.begin(), final_graph.end(), cmpEdge);
  vector<EdgeInfo> final_mst = kruskalAlgorithm(final_graph, final_sketches.size());
  vector<EdgeInfo>().swap(final_graph);
  if(is_newick_tree){
    string output_newick_file = output_file + ".newick.tree";
    print_newick_tree(final_sketches, final_mst, pre_sketch_by_file, output_newick_file);
    cerr << "-----write the newick tree into: " << output_newick_file << endl;

  }

  vector<EdgeInfo> forest = generateForest(final_mst, threshold);
  vector<vector<int>> tmpClust = generateClusterWithBfs(forest, final_sketches.size());
  printResult(tmpClust, final_sketches, pre_sketch_by_file, output_file);
  cerr << "-----write the cluster result into: " << output_file << endl;
  cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
  if(!no_save){
    saveMST(final_sketches, final_mst, new_folder_path, sketch_by_file);
  }

  if(!no_dense){
    loadANI(folder_path, pre_ani_arr, sketch_func_id_0);
    for(int i = 0; i < 101; i++)
      ani_arr[i] += pre_ani_arr[i];
    for(int i = 0; i < pre_dense_span; i++){
      for(int j = 0; j < pre_genome_number; j++){
        dense_arr[i][j] += pre_dense_arr[i][j];
      }
    }
    if(!no_save){
      saveANI(new_folder_path, ani_arr, sketch_func_id_0);
      saveDense(new_folder_path, dense_arr, dense_span, final_sketches.size());
    }

    int alpha = 2;
    int denseIndex = threshold / 0.01;
    vector<int> totalNoiseArr;
    for(int i = 0; i < tmpClust.size(); i++){
      if(tmpClust[i].size() == 1) continue;
      vector<PairInt> curDenseArr;
      set<int> denseSet;
      for(int j = 0; j < tmpClust[i].size(); j++){
        int element = tmpClust[i][j];
        PairInt p(element, dense_arr[denseIndex][element]);
        denseSet.insert(dense_arr[denseIndex][element]);
        curDenseArr.push_back(p);
      }
      vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
      totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
    }
    cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
    forest = modifyForest(forest, totalNoiseArr, threads);
    vector<vector<int>> cluster = generateClusterWithBfs(forest, final_sketches.size());
    string outputFileNew = output_file + ".removeNoise";
    printResult(cluster, final_sketches, pre_sketch_by_file, outputFileNew);
    cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
    cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
  }

}
void clust_from_mst_fast(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads){
  vector<KssdSketchInfo> sketches;
  vector<EdgeInfo> mst;
  vector<vector<int>> cluster;
  bool sketchByFile = load_kssd_genome_info(folder_path, "mst", sketches);
  load_MST(folder_path, mst);

  if(is_newick_tree){
    string output_newick_file = outputFile + ".newick.tree";
    print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
    cerr << "-----write the newick tree into: " << output_newick_file << endl;
  }

  vector<EdgeInfo> forest = generateForest(mst, threshold);
  vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
  printKssdResult(tmpClust, sketches, sketchByFile, outputFile);
  cerr << "-----write the cluster result into: " << outputFile << endl;
  cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;
  if(!no_dense){
    int **denseArr;
    int genome_number = sketches.size();
    int denseSpan = DENSE_SPAN;
    loadDense(denseArr, folder_path, denseSpan, genome_number);
    int alpha = 2;
    int denseIndex = threshold / 0.01;
    vector<int> totalNoiseArr;
    for(int i = 0; i < tmpClust.size(); i++){
      if(tmpClust[i].size() == 1) continue;
      vector<PairInt> curDenseArr;
      set<int> denseSet;
      for(int j = 0; j < tmpClust[i].size(); j++){
        int element = tmpClust[i][j];
        PairInt p(element, denseArr[denseIndex][element]);
        denseSet.insert(denseArr[denseIndex][element]);
        curDenseArr.push_back(p);
      }
      vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
      totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
    }
    cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
    forest = modifyForest(forest, totalNoiseArr, threads);
    cluster = generateClusterWithBfs(forest, sketches.size());
    string outputFileNew = outputFile + ".removeNoise";
    printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
    cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
    cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
  }
}

void clust_from_mst(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads){
  vector<SketchInfo> sketches;
  vector<EdgeInfo> mst;
  vector<vector<int>> cluster;
  bool sketchByFile = load_genome_info(folder_path, "mst", sketches);
  load_MST(folder_path, mst);

  if(is_newick_tree){
    string output_newick_file = outputFile + ".newick.tree";
    print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
    cerr << "-----write the newick tree into: " << output_newick_file << endl;
  }

  vector<EdgeInfo> forest = generateForest(mst, threshold);
  vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
  printResult(tmpClust, sketches, sketchByFile, outputFile);
  cerr << "-----write the cluster result into: " << outputFile << endl;
  cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

  if(!no_dense){
    int **denseArr;
    int genome_number = sketches.size();
    int denseSpan = DENSE_SPAN;
    loadDense(denseArr, folder_path, denseSpan, genome_number);
    int alpha = 2;
    int denseIndex = threshold / 0.01;
    vector<int> totalNoiseArr;
    for(int i = 0; i < tmpClust.size(); i++){
      if(tmpClust[i].size() == 1) continue;
      vector<PairInt> curDenseArr;
      set<int> denseSet;
      for(int j = 0; j < tmpClust[i].size(); j++){
        int element = tmpClust[i][j];
        PairInt p(element, denseArr[denseIndex][element]);
        denseSet.insert(denseArr[denseIndex][element]);
        curDenseArr.push_back(p);
      }
      vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
      totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
    }
    cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
    forest = modifyForest(forest, totalNoiseArr, threads);
    cluster = generateClusterWithBfs(forest, sketches.size());
    string outputFileNew = outputFile + ".removeNoise";
    printResult(cluster, sketches, sketchByFile, outputFileNew);
    cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
    cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
  }
}
#endif

void clust_from_genome_fast(const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads){
  bool isSave = !noSave;
  vector<KssdSketchInfo> sketches;
  KssdParameters info;
  compute_kssd_sketches(sketches, info, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);
#ifndef GREEDY_CLUST	
  transSketches(sketches, info, folder_path, threads);
#else
#endif
  compute_kssd_clusters(sketches, info, sketchByFile, no_dense, isContainment, folder_path, outputFile, is_newick_tree, threshold, isSave, threads);

}

void compute_kssd_clusters(vector<KssdSketchInfo>& sketches, const KssdParameters info, bool sketchByFile, bool no_dense, bool isContainment, const string folder_path, string outputFile, bool is_newick_tree, double threshold, bool isSave, int threads){
  vector<vector<int>> cluster;
  double t2 = get_sec();

#ifdef GREEDY_CLUST
  //======clust-greedy====================================================================
  int sketch_func_id = 0;
  cluster = KssdgreedyCluster(sketches, sketch_func_id, threshold, threads);
  printKssdResult(cluster, sketches, sketchByFile, outputFile);
  cerr << "-----write the cluster result into: " << outputFile << endl;
  cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
  double t3 = get_sec();
#ifdef Timer
  cerr << "========time of greedyCluster is: " << t3 - t2 << "========" << endl;
#endif
  //======clust-greedy====================================================================
#else



  //======clust-mst=======================================================================
  int **denseArr;
  uint64_t* aniArr; //= new uint64_t[101];
  int denseSpan = DENSE_SPAN;
  vector<EdgeInfo> mst = compute_kssd_mst(sketches, 0,0,info, folder_path,  no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold);
  double t3 = get_sec();
#ifdef Timer
  cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
#endif
  //string sketch_func_id = "fast_kssd_sketch";
  if(isSave){
    if(!no_dense){
      saveANI(folder_path, aniArr, 1);
      saveDense(folder_path, denseArr, denseSpan, sketches.size());
    }
    saveKssdMST(sketches, mst, folder_path, sketchByFile);
  }
  double t4 = get_sec();
#ifdef Timer
  cerr << "========time of saveMST is: " << t4 - t3 << "========" << endl;
#endif

  //generate the Newick tree format
  if(is_newick_tree){
    string output_newick_file = outputFile + ".newick.tree";
    print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
    cerr << "-----write the newick tree into: " << output_newick_file << endl;
  }

  vector<EdgeInfo> forest = generateForest(mst, threshold);
  vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
  printKssdResult(tmpClust, sketches, sketchByFile, outputFile);
  cerr << "-----write the cluster result into: " << outputFile << endl;
  cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;
  //tune cluster by noise cluster
  if(!no_dense){
    int alpha = 2;
    int denseIndex = threshold / 0.01;
    vector<int> totalNoiseArr;
    for(int i = 0; i < tmpClust.size(); i++){
      if(tmpClust[i].size() == 1) continue;
      vector<PairInt> curDenseArr;
      set<int> denseSet;
      for(int j = 0; j < tmpClust[i].size(); j++){
        int element = tmpClust[i][j];
        PairInt p(element, denseArr[denseIndex][element]);
        denseSet.insert(denseArr[denseIndex][element]);
        curDenseArr.push_back(p);
      }
      vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
      totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
    }
    cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
    forest = modifyForest(forest, totalNoiseArr, threads);
    cluster = generateClusterWithBfs(forest, sketches.size());
    string outputFileNew = outputFile + ".removeNoise";
    printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
    cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
    cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
    double t5 = get_sec();
#ifdef Timer
    cerr << "========time of tuning cluster is: " << t5 - t4 << "========" << endl;
#endif
  }
  //======clust-mst=======================================================================
#endif//endif GREEDY_CLUST
}

void compute_kssd_sketches(vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads){
  double t0 = get_sec();
  if(sketchByFile){
    //if(!sketchFiles(inputFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){

    fprintf(stderr, "-----input fileList, sketch by file\n");
    ifstream ifs(inputFile);
    if(!ifs){
      fprintf(stderr, "error open the inputFile: %s\n", inputFile.c_str());
      exit(1);
    }
    vector<string> fileList;
    string fileName;
    while(getline(ifs, fileName)){
      fileList.push_back(fileName);
    }


    if(!sketchFileWithKssd(fileList, minLen, kmerSize, drlevel, sketches, info, threads)){
      cerr << "ERROR: sketchFileWithKssd(), cannot finish the sketch generation by genome files" << endl;
      exit(1);
    }
  }//end sketch by sequence
    else{
      cerr << "use the sketch seuqnce with kssd " << endl;
      if(!sketchSequencesWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads)){
        cerr << "ERROR: sketchSequencesWithKssd (), cannot finish the sketch generation by genome sequences" << endl;
        exit(1);
      }
    }//end sketch by file
    cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;
    double t1 = get_sec();
#ifdef Timer
    cerr << "========time of computing sketch is: " << t1 - t0 << "========" << endl;
#endif
    folder_path = currentDataTime();
    if(isSave){
      string command = "mkdir -p " + folder_path;
      system(command.c_str());
      saveKssdSketches(sketches, info, folder_path, sketchByFile);
      double t2 = get_sec();
#ifdef Timer
      cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
    }

  }

  void clust_from_genomes(string inputFile, string outputFile, bool is_newick_tree, bool sketchByFile, bool no_dense, int kmerSize, int sketchSize, double threshold, string sketchFunc, bool isContainment, int containCompress, int minLen, string folder_path, bool noSave, int threads){
    bool isSave = !noSave;
    vector<SketchInfo> sketches;
    int sketch_func_id;
    if(sketchFunc == "MinHash")	sketch_func_id = 0;
    else if(sketchFunc == "KSSD") sketch_func_id = 1;

    compute_sketches(sketches, inputFile, folder_path, sketchByFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, isSave, threads);

    compute_clusters(sketches, sketchByFile, outputFile, is_newick_tree, no_dense, folder_path, sketch_func_id, threshold, isSave, threads);
  }

  bool tune_kssd_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, int& kmerSize, double& threshold, int &drlevel){
    uint64_t maxSize, minSize, averageSize;
    calSize(sketchByFile, inputFile, threads, minLen, maxSize, minSize, averageSize);

    //======tune the sketch_size===============
    int compression = 1 << (4 * drlevel);
    int sketchSize = averageSize / compression;
    //=====tune the kmer_size===============
    double warning_rate = 0.01;
    double recommend_rate = 0.0001;
    int alphabetSize = 4;//for "AGCT"
    int recommendedKmerSize = ceil(log(maxSize * (1 - recommend_rate) / recommend_rate) / log(4));
    int warningKmerSize = ceil(log(maxSize * (1 - warning_rate) / warning_rate) / log(4));
    if(!isSetKmer){
      kmerSize = recommendedKmerSize;
    }
    else{
      if(kmerSize < warningKmerSize){
        cerr << "the kmerSize " << kmerSize << " is too small for the maximum genome size of " << maxSize << endl;
        cerr << "replace the kmerSize to the: " << recommendedKmerSize << " for reducing the random collision of kmers" << endl;
        kmerSize = recommendedKmerSize;
      }
      else if(kmerSize > recommendedKmerSize + 3){
        cerr << "the kmerSize " << kmerSize << " maybe too large for the maximum genome size of " << maxSize << endl;
        cerr << "replace the kmerSize to the " << recommendedKmerSize << " for increasing the sensitivity of genome comparison" << endl;
        kmerSize = recommendedKmerSize;
      }
    }

    //=====tune the distance threshold===============
    double minJaccard = 0.001;
    if(!isContainment){
      minJaccard = 1.0 / sketchSize;
    }
    else{
      //minJaccard = 1.0 / (averageSize / containCompress);
      minJaccard = 1.0 / (minSize / compression);
    }

    double maxDist;
    if(minJaccard >= 1.0)	
      maxDist = 1.0;
    else
      maxDist = -1.0 / kmerSize * log(2*minJaccard / (1.0 + minJaccard));
    cerr << "-----the max recommand distance threshold is: " << maxDist << endl;
    if(threshold > maxDist){
      cerr << "ERROR: tune_parameters(), the threshold: " << threshold << " is out of the valid distance range estimated by Mash distance or AAF distance" << endl;
      cerr << "Please set a distance threshold with -d option" << endl;
      return false;
    }

#ifdef DEBUG
    if(sketchByFile) cerr << "-----sketch by file!" << endl;
    else cerr << "-----sketch by sequence!" << endl;
    cerr << "-----the kmerSize is: " << kmerSize << endl;
    cerr << "-----the thread number is: " << threads << endl;
    cerr << "-----the threshold is: " << threshold << endl;
    if(isContainment)
      cerr << "-----use the AAF distance, the sketchSize is approximately in proportion with 1/" << compression << endl;
    else
      cerr << "-----use the Mash distance, the sketchSize is about: " << sketchSize << endl;
#endif

    return true;
  }

  bool tune_parameters(bool sketchByFile, bool isSetKmer, string inputFile, int threads, int minLen, bool& isContainment, bool& isJaccard, int& kmerSize, double& threshold, int& containCompress, int& sketchSize){
    uint64_t maxSize, minSize, averageSize;
    calSize(sketchByFile, inputFile, threads, minLen, maxSize, minSize, averageSize);

    //======tune the sketch_size===============
    if(isContainment && isJaccard){
      cerr << "ERROR: tune_parameters(), conflict distance measurement of Mash distance (fixed-sketch-size) and AAF distance (variable-sketch-size) " << endl;
      return false;
    }
#ifdef GREEDY_CLUST
    //======clust-greedy====================================================================
    if(!isContainment && !isJaccard){
      containCompress = averageSize / 1000;
      isContainment = true;
    }
    else if(!isContainment && isJaccard){
      //do nothing
    }
    else{
      if(averageSize / containCompress < 10){
        cerr << "the containCompress " << containCompress << " is too large and the sketch size is too small" << endl;
        containCompress = averageSize / 1000;
        cerr << "set the containCompress to: " << containCompress << endl;
      }
    }
    //=======clust-greedy===================================================================
#endif
    //=====tune the kmer_size===============
    double warning_rate = 0.01;
    double recommend_rate = 0.0001;
    int alphabetSize = 4;//for "AGCT"
    int recommendedKmerSize = ceil(log(maxSize * (1 - recommend_rate) / recommend_rate) / log(4));
    int warningKmerSize = ceil(log(maxSize * (1 - warning_rate) / warning_rate) / log(4));
    if(!isSetKmer){
      kmerSize = recommendedKmerSize;
    }
    else{
      if(kmerSize < warningKmerSize){
        cerr << "the kmerSize " << kmerSize << " is too small for the maximum genome size of " << maxSize << endl;
        cerr << "replace the kmerSize to the: " << recommendedKmerSize << " for reducing the random collision of kmers" << endl;
        kmerSize = recommendedKmerSize;
      }
      else if(kmerSize > recommendedKmerSize + 3){
        cerr << "the kmerSize " << kmerSize << " maybe too large for the maximum genome size of " << maxSize << endl;
        cerr << "replace the kmerSize to the " << recommendedKmerSize << " for increasing the sensitivity of genome comparison" << endl;
        kmerSize = recommendedKmerSize;
      }
    }

    //=====tune the distance threshold===============
    double minJaccard = 0.001;
    if(!isContainment){
      minJaccard = 1.0 / sketchSize;
    }
    else{
      //minJaccard = 1.0 / (averageSize / containCompress);
      minJaccard = 1.0 / (minSize / containCompress);
    }

    double maxDist;
    if(minJaccard >= 1.0)	
      maxDist = 1.0;
    else
      maxDist = -1.0 / kmerSize * log(2*minJaccard / (1.0 + minJaccard));
    cerr << "-----the max recommand distance threshold is: " << maxDist << endl;
    if(threshold > maxDist){
      cerr << "ERROR: tune_parameters(), the threshold: " << threshold << " is out of the valid distance range estimated by Mash distance or AAF distance" << endl;
      cerr << "Please set a distance threshold with -d option" << endl;
      return false;
    }

#ifdef DEBUG
    if(sketchByFile) cerr << "-----sketch by file!" << endl;
    else cerr << "-----sketch by sequence!" << endl;
    cerr << "-----the kmerSize is: " << kmerSize << endl;
    cerr << "-----the thread number is: " << threads << endl;
    cerr << "-----the threshold is: " << threshold << endl;
    if(isContainment)
      cerr << "-----use the AAF distance (variable-sketch-size), the sketchSize is in proportion with 1/" << containCompress << endl;
    else
      cerr << "-----use the Mash distance (fixed-sketch-size), the sketchSize is: " << sketchSize << endl;
#endif

    return true;
  }

  void clust_from_sketch_fast(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, bool isContainment, double threshold, int threads){
    vector<KssdSketchInfo> sketches;
    vector<vector<int>> cluster;
    bool sketchByFile;
    KssdParameters info;
#ifdef GREEDY_CLUST
    //======clust-greedy====================================================================
    double time0 = get_sec();
    sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
    cerr << "-----the size of sketches is: " << sketches.size() << endl;
    double time1 = get_sec();
#ifdef Timer
    cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
    cluster = KssdgreedyCluster(sketches, 0, threshold, threads);
    printKssdResult(cluster, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
    double time2 = get_sec();
#ifdef Timer
    cerr << "========time of greedy incremental cluster is: " << time2 - time1 << endl;
#endif
    //======clust-greedy====================================================================
#else     


    //======clust-mst=======================================================================
    double time0 = get_sec();
    //sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
    sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);

    cerr << "-----the size of sketches is: " << sketches.size() << endl;
    double time1 = get_sec();
#ifdef Timer
    cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
    int** denseArr;
    uint64_t* aniArr; //= new uint64_t[101];
    int denseSpan = DENSE_SPAN;
    //vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, denseArr, denseSpan, aniArr);
    vector<EdgeInfo> mst = compute_kssd_mst(sketches, 0,0,info, folder_path,  no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold);
    double time2 = get_sec();
#ifdef Timer
    cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
#endif
    if(is_newick_tree){
      string output_newick_file = outputFile + ".newick.tree";
      print_kssd_newick_tree(sketches, mst, sketchByFile, output_newick_file);
      cerr << "-----write the newick tree into: " << output_newick_file << endl;
    }
    vector<EdgeInfo> forest = generateForest(mst, threshold);
    vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
    printKssdResult(tmpClust, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

    if(!no_dense){
      int alpha = 2;
      int denseIndex = threshold / 0.01;
      vector<int> totalNoiseArr;
      for(int i = 0; i < tmpClust.size(); i++){
        if(tmpClust[i].size() == 1) continue;
        vector<PairInt> curDenseArr;
        set<int> denseSet;
        for(int j = 0; j < tmpClust[i].size(); j++){
          int element = tmpClust[i][j];
          PairInt p(element, denseArr[denseIndex][element]);
          denseSet.insert(denseArr[denseIndex][element]);
          curDenseArr.push_back(p);
        }
        vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
        totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
      }
      cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
      forest = modifyForest(forest, totalNoiseArr, threads);
      cluster = generateClusterWithBfs(forest, sketches.size());
      string outputFileNew = outputFile + ".removeNoise";
      printKssdResult(cluster, sketches, sketchByFile, outputFileNew);
      cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
      cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
    }
    double time3 = get_sec();
#ifdef Timer
    cerr << "========time of generator forest and cluster is: " << time3 - time2 << "========" << endl;
#endif
#endif
  }

  void clust_from_sketches(string folder_path, string outputFile, bool is_newick_tree, bool no_dense, double threshold, int threads){
    vector<SketchInfo> sketches;
    vector<vector<int>> cluster;
    int sketch_func_id;
    bool sketchByFile;
#ifdef GREEDY_CLUST
    //======clust-greedy====================================================================
    double time0 = get_sec();
    sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
    cerr << "-----the size of sketches is: " << sketches.size() << endl;
    double time1 = get_sec();
#ifdef Timer
    cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
    cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
    printResult(cluster, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
    double time2 = get_sec();
#ifdef Timer
    cerr << "========time of greedy incremental cluster is: " << time2 - time1 << endl;
#endif
    //======clust-greedy====================================================================
#else
    //======clust-mst=======================================================================
    double time0 = get_sec();
    sketchByFile = loadSketches(folder_path, threads, sketches, sketch_func_id);
    cerr << "-----the size of sketches is: " << sketches.size() << endl;
    double time1 = get_sec();
#ifdef Timer
    cerr << "========time of load genome Infos and sketch Infos is: " << time1 - time0 << endl;
#endif
    int** denseArr;
    uint64_t* aniArr; //= new uint64_t[101];
    int denseSpan = DENSE_SPAN;
    vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, no_dense, denseArr, denseSpan, aniArr);
    double time2 = get_sec();
#ifdef Timer
    cerr << "========time of generateMST is: " << time2 - time1 << "========" << endl;
#endif
    if(is_newick_tree){
      string output_newick_file = outputFile + ".newick.tree";
      print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
      cerr << "-----write the newick tree into: " << output_newick_file << endl;
    }
    vector<EdgeInfo> forest = generateForest(mst, threshold);
    vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
    printResult(tmpClust, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

    if(!no_dense){
      int alpha = 2;
      int denseIndex = threshold / 0.01;
      vector<int> totalNoiseArr;
      for(int i = 0; i < tmpClust.size(); i++){
        if(tmpClust[i].size() == 1) continue;
        vector<PairInt> curDenseArr;
        set<int> denseSet;
        for(int j = 0; j < tmpClust[i].size(); j++){
          int element = tmpClust[i][j];
          PairInt p(element, denseArr[denseIndex][element]);
          denseSet.insert(denseArr[denseIndex][element]);
          curDenseArr.push_back(p);
        }
        vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
        totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
      }
      cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
      forest = modifyForest(forest, totalNoiseArr, threads);
      cluster = generateClusterWithBfs(forest, sketches.size());
      string outputFileNew = outputFile + ".removeNoise";
      printResult(cluster, sketches, sketchByFile, outputFileNew);
      cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
      cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
    }
    double time3 = get_sec();
#ifdef Timer
    cerr << "========time of generator forest and cluster is: " << time3 - time2 << "========" << endl;
#endif
    //=======clust-mst======================================================================
#endif
  }

  void compute_sketches(vector<SketchInfo>& sketches, string inputFile, string& folder_path, bool sketchByFile, int minLen, int kmerSize, int sketchSize, string sketchFunc, bool isContainment, int containCompress,  bool isSave, int threads){
    double t0 = get_sec();
    if(sketchByFile){
      if(!sketchFiles(inputFile, minLen, kmerSize, sketchSize, sketchFunc, isContainment, containCompress, sketches, threads)){
        cerr << "ERROR: generate_sketches(), cannot finish the sketch generation by genome files" << endl;
        exit(1);
      }
    }//end sketch by sequence
    else{
      if(!sketchSequences(inputFile, kmerSize, sketchSize, minLen, sketchFunc, isContainment, containCompress, sketches, threads)){
        cerr << "ERROR: generate_sketches(), cannot finish the sketch generation by genome sequences" << endl;
        exit(1);
      }
    }//end sketch by file
    cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;
    double t1 = get_sec();
#ifdef Timer
    cerr << "========time of computing sketch is: " << t1 - t0 << "========" << endl;
#endif
    folder_path = currentDataTime();
    if(isSave){
      string command = "mkdir -p " + folder_path;
      system(command.c_str());
      saveSketches(sketches, folder_path, sketchByFile, sketchFunc, isContainment, containCompress, sketchSize, kmerSize);
      double t2 = get_sec();
#ifdef Timer
      cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
    }
  }

  void compute_clusters(vector<SketchInfo>& sketches, bool sketchByFile, string outputFile, bool is_newick_tree, bool no_dense, string folder_path, int sketch_func_id, double threshold, bool isSave, int threads){
    vector<vector<int>> cluster;
    double t2 = get_sec();
#ifdef GREEDY_CLUST
    //======clust-greedy====================================================================
    cluster = greedyCluster(sketches, sketch_func_id, threshold, threads);
    printResult(cluster, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of " << outputFile << " is: " << cluster.size() << endl;
    double t3 = get_sec();
#ifdef Timer
    cerr << "========time of greedyCluster is: " << t3 - t2 << "========" << endl;
#endif
    //======clust-greedy====================================================================
#else
    //======clust-mst=======================================================================
    int **denseArr;
    uint64_t* aniArr; //= new uint64_t[101];
    int denseSpan = DENSE_SPAN;
    vector<EdgeInfo> mst = modifyMST(sketches, 0, sketch_func_id, threads, no_dense, denseArr, denseSpan, aniArr);
    double t3 = get_sec();
#ifdef Timer
    cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
#endif
    if(isSave){
      saveMST(sketches, mst, folder_path, sketchByFile);
    }
    double t4 = get_sec();
#ifdef Timer
    cerr << "========time of saveMST is: " << t4 - t3 << "========" << endl;
#endif

    //generate the Newick tree format
    if(is_newick_tree){
      string output_newick_file = outputFile + ".newick.tree";
      print_newick_tree(sketches, mst, sketchByFile, output_newick_file);
      cerr << "-----write the newick tree into: " << output_newick_file << endl;
    }

    vector<EdgeInfo> forest = generateForest(mst, threshold);
    vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
    printResult(tmpClust, sketches, sketchByFile, outputFile);
    cerr << "-----write the cluster result into: " << outputFile << endl;
    cerr << "-----the cluster number of: " << outputFile << " is: " << tmpClust.size() << endl;

    //tune cluster by noise cluster
    if(!no_dense){
      if(isSave){
        saveANI(folder_path, aniArr, sketch_func_id);
        saveDense(folder_path, denseArr, denseSpan, sketches.size());
      }
      int alpha = 2;
      int denseIndex = threshold / 0.01;
      vector<int> totalNoiseArr;
      for(int i = 0; i < tmpClust.size(); i++){
        if(tmpClust[i].size() == 1) continue;
        vector<PairInt> curDenseArr;
        set<int> denseSet;
        for(int j = 0; j < tmpClust[i].size(); j++){
          int element = tmpClust[i][j];
          PairInt p(element, denseArr[denseIndex][element]);
          denseSet.insert(denseArr[denseIndex][element]);
          curDenseArr.push_back(p);
        }
        vector<int> curNoiseArr = getNoiseNode(curDenseArr, alpha);
        totalNoiseArr.insert(totalNoiseArr.end(), curNoiseArr.begin(), curNoiseArr.end());
      }
      cerr << "-----the total noiseArr size is: " << totalNoiseArr.size() << endl;
      forest = modifyForest(forest, totalNoiseArr, threads);
      cluster = generateClusterWithBfs(forest, sketches.size());
      string outputFileNew = outputFile + ".removeNoise";
      printResult(cluster, sketches, sketchByFile, outputFileNew);
      cerr << "-----write the cluster without noise into: " << outputFileNew << endl;
      cerr << "-----the cluster number of: " << outputFileNew << " is: " << cluster.size() << endl;
      double t5 = get_sec();
#ifdef Timer
      cerr << "========time of tuning cluster is: " << t5 - t4 << "========" << endl;
#endif
    }

    //======clust-mst=======================================================================
#endif//endif GREEDY_CLUST
  }



  void compute_kssd_sketches_mpi(int my_rank, vector<KssdSketchInfo>& sketches, KssdParameters& info, bool isSave, const string inputFile, string& folder_path, bool sketchByFile, const int minLen, const int kmerSize, const int drlevel, int threads) {
    double t0 = get_sec();

    if (sketchByFile) {
      vector<string> file_list;
      ifstream ifs(inputFile);
      if (!ifs.is_open()) {
        cerr << "ERROR: cannot open input file list: " << inputFile << endl;

      }
      string line;
      while (getline(ifs, line)) {
        if (!line.empty()) file_list.push_back(line);
      }
      ifs.close();

      int comm_sz;
      MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
      size_t total = file_list.size();
      size_t chunk_size = (total + comm_sz - 1) / comm_sz;
      size_t start = my_rank * chunk_size;
      size_t end = min(start + chunk_size, total);
      vector<string> local_files(file_list.begin() + start, file_list.begin() + end);
      cerr << "Rank " << my_rank << " is processing files from index " << start << " to " << end - 1 << " (total files: " << total << ")" << endl;
      if (!sketchFileWithKssd(local_files,  minLen, kmerSize, drlevel, sketches, info, threads)) {
        cerr << "ERROR: compute_sketches(), cannot finish sketch generation by genome files" << endl;
        exit(1);
      }

    } else {
      if (sketchSequencesWithKssd(inputFile, minLen, kmerSize, drlevel, sketches, info, threads)) {
        cerr << "ERROR: compute_sketches(), cannot finish sketch generation by genome sequences" << endl;
        exit(1);
      }
    }

    cerr << "-----the size of sketches (number of genomes or sequences) is: " << sketches.size() << endl;

    double t1 = get_sec();
#ifdef Timer
    cerr << "========time of computing sketch is: " << t1 - t0 << "========" << endl;
#endif

    folder_path = currentDataTime()+ to_string(my_rank);
    isSave = true;
    if (isSave) {
      string command = "mkdir -p " + folder_path;
      system(command.c_str());
      //saveSketches(sketches, folder_path, sketchByFile, sketchFunc, isContainment, containCompress, sketchSize, kmerSize);
      saveKssdSketches(sketches, info, folder_path, sketchByFile);

      double t2 = get_sec();
#ifdef Timer
      cerr << "========time of saveSketches is: " << t2 - t1 << "========" << endl;
#endif
    }
  }



  void distribute_compute_clusters(int my_rank, int comm_sz, vector<KssdSketchInfo>& sketches, const KssdParameters info,int start_index, int end_index, bool sketchByFile, string output_file, bool is_newick_tree, string folder_path, double threshold, bool isSave, int threads, bool no_dense, bool isContainment){
    vector<vector<int>> cluster;
    double t2 = get_sec();
    int **denseArr;
    uint64_t* aniArr; //= new uint64_t[101];
    int denseSpan = DENSE_SPAN;
    //======clust-mst=======================================================================
    vector<EdgeInfo> my_mst = compute_kssd_mst(sketches, start_index, end_index, info, folder_path, no_dense, isContainment, threads, denseArr, denseSpan, aniArr, threshold);
    double t3 = get_sec();
#ifdef Timer
    cerr << "========time of generateMST is: " << t3 - t2 << "========" << endl;
#endif
    //string tmp_recv_mst_path = folder_path + "_mst_" + to_string(my_rank);
    //string cmd0 = "mkdir -p " + tmp_recv_mst_path;
    //system(cmd0.c_str());
    isSave = true;
    if(isSave){
      saveKssdMST(sketches, my_mst, folder_path, sketchByFile);
    }
    double t4 = get_sec();
#ifdef Timer
    cerr << "========time of saveMST is: " << t4 - t3 << "========" << endl;
#endif

    string info_file = folder_path + '/' + "kssd.info.mst";
    string edge_file = folder_path + '/' + "edge.mst";
    char * info_buffer;
    char * edge_buffer;
    size_t info_size, edge_size;
    build_message(info_buffer, info_size, info_file);
    build_message(edge_buffer, edge_size, edge_file);
    cerr << "finish the build_message in distribute_build_MST " << my_rank << endl;

    string tmp_path = "tmp";
    string cmd1 = "mkdir -p " + tmp_path;
    system(cmd1.c_str());

    if(my_rank != 0){
      MPI_Send(&info_size, 1, MPI_UNSIGNED_LONG, 0, my_rank, MPI_COMM_WORLD);
      MPI_Send(&edge_size, 1, MPI_UNSIGNED_LONG, 0, my_rank+comm_sz, MPI_COMM_WORLD);
      MPI_Send(info_buffer, info_size, MPI_CHAR, 0, my_rank+comm_sz*2, MPI_COMM_WORLD);
      MPI_Send(edge_buffer, edge_size, MPI_CHAR, 0, my_rank+comm_sz*3, MPI_COMM_WORLD);
    }
    else{
      size_t recv_info_size, recv_edge_size;
      vector<EdgeInfo> sum_mst;
      sum_mst.insert(sum_mst.end(), my_mst.begin(), my_mst.end());
      for(int id = 1; id < comm_sz; id++){
        MPI_Recv(&recv_info_size, 1, MPI_UNSIGNED_LONG, id, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&recv_edge_size, 1, MPI_UNSIGNED_LONG, id, id+comm_sz, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        char * recv_info_buffer = new char[recv_info_size];
        char * recv_edge_buffer = new char[recv_edge_size];
        MPI_Recv(recv_info_buffer, recv_info_size, MPI_CHAR, id, id+comm_sz*2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recv_edge_buffer, recv_edge_size, MPI_CHAR, id, id+comm_sz*3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        vector<EdgeInfo> cur_MST;
        format_MST(id, recv_edge_buffer, recv_edge_size, tmp_path, cur_MST);
        sum_mst.insert(sum_mst.end(), cur_MST.begin(), cur_MST.end());
        vector<EdgeInfo>().swap(cur_MST);
      }
      sort(sum_mst.begin(), sum_mst.end(), cmpEdge);
      vector<EdgeInfo> final_mst = kruskalAlgorithm(sum_mst, sketches.size());
      vector<EdgeInfo>().swap(sum_mst);
      if(is_newick_tree){
        string output_newick_file = output_file + ".newick.tree";
        print_kssd_newick_tree(sketches, final_mst, sketchByFile, output_newick_file);
        cerr << "-----write the newick tree into: " << output_newick_file << endl;
      }
      vector<EdgeInfo> forest = generateForest(final_mst, threshold);
      vector<vector<int>> tmpClust = generateClusterWithBfs(forest, sketches.size());
      printKssdResult(tmpClust, sketches, sketchByFile, output_file);
      cerr << "-----write the cluster result into: " << output_file << endl;
      cerr << "-----the cluster number of: " << output_file << " is: " << tmpClust.size() << endl;
    }
    cerr << "finish the distribute_build_MST " << endl;

  }




  void clust_from_genomes_fast_MPI(int my_rank, int comm_sz, const string inputFile, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment, const int kmerSize, const double threshold, const int drlevel, const int minLen, bool noSave, int threads){
    bool isSave = !noSave;
    vector<KssdSketchInfo> sketches;
    int sketch_func_id;
    string sketchFunc;
    sketchFunc = "MinHash";
    KssdParameters info;
    if(sketchFunc == "MinHash")	sketch_func_id = 0;
    else if(sketchFunc == "KSSD") sketch_func_id = 1;
    compute_kssd_sketches_mpi(my_rank, sketches, info, isSave, inputFile, folder_path, sketchByFile, minLen, kmerSize, drlevel, threads);

    size_t* info_size_arr = new size_t[comm_sz];
    size_t* hash_size_arr = new size_t[comm_sz];
    size_t* index_size_arr = new size_t[comm_sz];
    string info_file = folder_path + '/' + "kssd.info.sketch";
    string hash_file = folder_path + '/' + "kssd.hash.sketch";
    string index_file = folder_path + '/' + "kssd.sketch.index";
    char * info_buffer;
    char * hash_buffer;
    string final_folder;
    size_t info_size, hash_size, index_size;
    build_message(info_buffer, info_size, info_file);
    build_message(hash_buffer, hash_size, hash_file);
    cerr << "=====================finished the build_message " << my_rank << endl;
    info_size_arr[my_rank] = info_size;
    hash_size_arr[my_rank] = hash_size;
    if(my_rank != 0){
      //MPI_Send(info_buffer, info_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      //MPI_Send(hash_buffer, hash_size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&info_size, 1, MPI_UNSIGNED_LONG, 0, my_rank, MPI_COMM_WORLD);
      MPI_Send(&hash_size, 1, MPI_UNSIGNED_LONG, 0, my_rank+comm_sz, MPI_COMM_WORLD);
      MPI_Send(info_buffer, info_size, MPI_CHAR, 0, my_rank+comm_sz*2, MPI_COMM_WORLD);
      MPI_Send(hash_buffer, hash_size, MPI_CHAR, 0, my_rank+comm_sz*3, MPI_COMM_WORLD);

      size_t sum_info_size, sum_hash_size, sum_index_size, sum_dict_size;
      MPI_Recv(&sum_info_size, 1, MPI_UNSIGNED_LONG, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&sum_index_size, 1, MPI_UNSIGNED_LONG, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&sum_hash_size, 1, MPI_UNSIGNED_LONG, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&sum_dict_size, 1, MPI_UNSIGNED_LONG, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      char * sum_info_buffer = new char[sum_info_size];
      char * sum_hash_buffer = new char[sum_hash_size];
      char * sum_index_buffer = new char[sum_index_size];
      char * sum_dict_buffer = new char[sum_dict_size];
      MPI_Recv(sum_info_buffer, sum_info_size, MPI_CHAR, 0, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(sum_hash_buffer, sum_hash_size, MPI_CHAR, 0, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(sum_index_buffer, sum_index_size, MPI_CHAR, 0, 203, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(sum_dict_buffer, sum_dict_size, MPI_CHAR, 0, 204, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //string sum_folder = folder_path + '_' + "sum_" + to_string(my_rank);
      //string cmd1 = "mkdir -p " + sum_folder;
      //system(cmd1.c_str());
      vector<KssdSketchInfo>().swap(sketches);
      format_sketches_index(sum_info_buffer, sum_info_size, sum_hash_buffer, sum_hash_size, sum_index_buffer, sum_index_size, sum_dict_buffer, sum_dict_size, folder_path, sketches, sketchByFile, threads);
    }
    else{
      //string tmp_recv_folder_path = folder_path + '_' + to_string(my_rank);
      //string cmd0 = "mkdir -p " + tmp_recv_folder_path;
      //system(cmd0.c_str());
      for(int id = 1; id < comm_sz; id++){
        MPI_Recv(&info_size_arr[id], 1, MPI_UNSIGNED_LONG, id, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&hash_size_arr[id], 1, MPI_UNSIGNED_LONG, id, id+comm_sz, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //cout << info_size_arr[id] << '\t' << hash_size_arr[id] << endl;
        char * recv_info_buffer = new char[info_size_arr[id]];
        char * recv_hash_buffer = new char[hash_size_arr[id]];
        MPI_Recv(recv_info_buffer, info_size_arr[id], MPI_CHAR, id, id+comm_sz*2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recv_hash_buffer, hash_size_arr[id], MPI_CHAR, id, id+comm_sz*3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        vector<KssdSketchInfo> cur_sketches;
        format_sketches(recv_info_buffer, info_size_arr[id], recv_hash_buffer, hash_size_arr[id], folder_path, cur_sketches, sketchByFile, threads);
        sketches.insert(sketches.end(), cur_sketches.begin(), cur_sketches.end());
        vector<KssdSketchInfo>().swap(cur_sketches);
      }
      //string sum_folder = folder_path + '_' + "sum";
      //final_folder = sum_folder;
      //string cmd1 = "mkdir -p " + sum_folder;
      //system(cmd1.c_str());
      //saveSketches(sketches, sum_folder, sketchByFile, sketchFunc, isContainment, containCompress, sketchSize, kmerSize);
      saveKssdSketches(sketches, info, folder_path, sketchByFile);
      transSketches(sketches, info, folder_path, threads);
      string sum_info_file = folder_path + '/' + "kssd.info.sketch";
      string sum_hash_file = folder_path + '/' + "kssd.hash.sketch";
      string sum_index_file = folder_path + '/' + "kssd.sketch.index";
      string sum_dict_file = folder_path + '/' + "kssd.sketch.dict";
      char * sum_info_buffer;
      char * sum_hash_buffer;
      char * sum_index_buffer;
      char * sum_dict_buffer;
      size_t sum_info_size, sum_hash_size, sum_index_size, sum_dict_size;
      build_message(sum_info_buffer, sum_info_size, sum_info_file);
      build_message(sum_index_buffer, sum_index_size, sum_index_file);
      build_message(sum_hash_buffer, sum_hash_size, sum_hash_file);
      build_message(sum_dict_buffer, sum_dict_size, sum_dict_file);
      for(int id = 1; id < comm_sz; id++){
        MPI_Send(&sum_info_size, 1, MPI_UNSIGNED_LONG, id, 101, MPI_COMM_WORLD);
        MPI_Send(&sum_index_size, 1, MPI_UNSIGNED_LONG, id, 102, MPI_COMM_WORLD);
        MPI_Send(&sum_hash_size, 1, MPI_UNSIGNED_LONG, id, 103, MPI_COMM_WORLD);
        MPI_Send(&sum_dict_size, 1, MPI_UNSIGNED_LONG, id, 104, MPI_COMM_WORLD);
        MPI_Send(sum_info_buffer, sum_info_size, MPI_CHAR, id, 201, MPI_COMM_WORLD);
        MPI_Send(sum_hash_buffer, sum_hash_size, MPI_CHAR, id, 202, MPI_COMM_WORLD);
        MPI_Send(sum_index_buffer, sum_index_size, MPI_CHAR, id, 203, MPI_COMM_WORLD);
        MPI_Send(sum_dict_buffer, sum_dict_size, MPI_CHAR, id, 204, MPI_COMM_WORLD);
      }

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    cerr << "the sketch size of rank " << my_rank << " is: " << sketches.size() << endl;





    int start_index, end_index;
    if(my_rank != 0){
      MPI_Send(&threads, 1, MPI_INT, 0, my_rank+comm_sz*4, MPI_COMM_WORLD);
      MPI_Recv(&start_index, 1, MPI_UNSIGNED_LONG, 0, my_rank+comm_sz*5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&end_index, 1, MPI_UNSIGNED_LONG, 0, my_rank+comm_sz*6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
      int* thread_arr = new int[comm_sz];
      thread_arr[0] = threads;
      int total_threads = 0;
      total_threads += threads;
      for(int id = 1; id < comm_sz; id++){
        MPI_Recv(&thread_arr[id], 1, MPI_INT, id, id+comm_sz*4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        total_threads += thread_arr[id];
      }
      cerr << "the total threads is: " << total_threads << endl;
      size_t N = sketches.size();
      size_t total_distance_time = N*(N-1)/2;
      //size_t total_distance_time = N;
      size_t* task_arr = new size_t[comm_sz];
      for(int i = 0; i < comm_sz; i++){
        task_arr[i] = total_distance_time * thread_arr[i] / total_threads;
      }
      size_t* start_index_arr = new size_t[comm_sz];
      size_t* end_index_arr = new size_t[comm_sz];
      int cur_start_index = 0;
      int cur_end_index = N-1;
      size_t global_index = 0;
      for(int i = 0; i < comm_sz; i++){
        start_index_arr[i] = global_index;
        size_t accumulate_time = 0;
        while(accumulate_time < task_arr[i] && global_index < N){
          //accumulate_time += N - global_index;
          accumulate_time += global_index + 1;
          global_index++;
        }
        end_index_arr[i] = global_index;
      }
      end_index_arr[comm_sz-1] = max(end_index_arr[comm_sz-1], N);
      start_index = start_index_arr[0];
      end_index = end_index_arr[0];
      for(int i = 1; i < comm_sz; i++){
        cerr << "start_index: " << start_index_arr[i] << '\t' << "end_index: " << end_index_arr[i] << endl;
        MPI_Send(&start_index_arr[i], 1, MPI_UNSIGNED_LONG, i, i+comm_sz*5, MPI_COMM_WORLD);
        MPI_Send(&end_index_arr[i], 1, MPI_UNSIGNED_LONG, i, i+comm_sz*6, MPI_COMM_WORLD);
      }
    }
    string tmp_str = "++++++++++++my rank: " + to_string(my_rank) + ", start_index: " + to_string(start_index) + ", end_index: " + to_string(end_index);
    cout << tmp_str << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    //compute_clusters(sketches, sketchByFile, outputFile, is_newick_tree, folder_path, sketch_func_id, threshold, isSave, threads);
    distribute_compute_clusters(my_rank, comm_sz, sketches, info, start_index, end_index, sketchByFile, outputFile, is_newick_tree, folder_path, threshold, isSave, threads, no_dense, isContainment);
    MPI_Finalize();
  }


  void clust_from_sketches_fast_MPI(int my_rank, int comm_sz, int half_k, int drlevel, string outputFile, string folder_path, bool is_newick_tree, bool no_dense, bool sketchByFile, bool isContainment,const double threshold, bool noSave, int threads)
  {
    bool isSave = !noSave;
    vector<KssdSketchInfo> sketches;
    int sketch_func_id;
    string sketchFunc;
    sketchFunc = "MinHash";
    KssdParameters info;
    if (sketchFunc == "MinHash") sketch_func_id = 0;
    size_t* info_size_arr = new size_t[comm_sz];
    size_t* hash_size_arr = new size_t[comm_sz];
    size_t* index_size_arr = new size_t[comm_sz];
    string info_file = folder_path + '/' + "kssd.info.sketch";
    string hash_file = folder_path + '/' + "kssd.hash.sketch";
    string index_file = folder_path + '/' + "kssd.sketch.index";
    char* info_buffer;
    char* hash_buffer;
    string final_folder;
    size_t info_size, hash_size, index_size;
    info.half_k = 10;
    info.drlevel = 3;


    build_message(info_buffer, info_size, info_file);
    build_message(hash_buffer, hash_size, hash_file);
    cerr << "=====================finished the build_message " << my_rank << endl;
    info_size_arr[my_rank] = info_size;
    hash_size_arr[my_rank] = hash_size;

    // 释放这部分不再需要的内存
    delete[] info_buffer;
    delete[] hash_buffer;
    info_buffer = nullptr;
    hash_buffer = nullptr;


    // Broadcast Logic
    size_t sum_info_size = 0, sum_hash_size = 0, sum_index_size = 0, sum_dict_size = 0;
    char* sum_info_buffer = nullptr;
    char* sum_hash_buffer = nullptr;
    char* sum_index_buffer = nullptr;
    char* sum_dict_buffer = nullptr;

    if (my_rank == 0) {
      sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
      transSketches(sketches, info, folder_path, threads);

      string sum_info_file = folder_path + '/' + "kssd.info.sketch";
      string sum_hash_file = folder_path + '/' + "kssd.hash.sketch";
      string sum_index_file = folder_path + '/' + "kssd.sketch.index";
      string sum_dict_file = folder_path + '/' + "kssd.sketch.dict";

      // Rank 0 reads files to determine sizes and get content
      build_message(sum_info_buffer, sum_info_size, sum_info_file);
      build_message(sum_index_buffer, sum_index_size, sum_index_file);
      build_message(sum_hash_buffer, sum_hash_size, sum_hash_file);
      build_message(sum_dict_buffer, sum_dict_size, sum_dict_file);
    }

    // Step 1: Rank 0 broadcasts the sizes of the buffers to all other processes.
    MPI_Bcast(&sum_info_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum_index_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum_hash_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum_dict_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Step 2: All processes (including rank 0) allocate memory based on the received sizes.
    if (my_rank != 0) {
      sum_info_buffer = new char[sum_info_size];
      sum_index_buffer = new char[sum_index_size];
      sum_hash_buffer = new char[sum_hash_size];
      sum_dict_buffer = new char[sum_dict_size];
    }

    // Step 3: Rank 0 broadcasts the actual data to all other processes.
    MPI_Bcast(sum_info_buffer, sum_info_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(sum_index_buffer, sum_index_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(sum_hash_buffer, sum_hash_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(sum_dict_buffer, sum_dict_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    char sketch_by_file_char = sketchByFile;
    MPI_Bcast(&sketch_by_file_char, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    sketchByFile = sketch_by_file_char;

    // Now, all processes have a copy of the data.
    // Non-zero ranks can proceed with their work.
    if (my_rank != 0) {
      vector<KssdSketchInfo>().swap(sketches);
      format_sketches_index(sum_info_buffer, sum_info_size,
          sum_hash_buffer, sum_hash_size,
          sum_index_buffer, sum_index_size,
          sum_dict_buffer, sum_dict_size,
          folder_path, sketches, sketchByFile, threads);
    }

    // CRITICAL: Clean up all allocated memory before exiting.
    delete[] sum_info_buffer;
    delete[] sum_index_buffer;
    delete[] sum_hash_buffer;
    delete[] sum_dict_buffer;
    delete[] index_size_arr; // Don't forget this one from the beginning 


    //build_message(info_buffer, info_size, info_file);
    //build_message(hash_buffer, hash_size, hash_file);
    //cerr << "=====================finished the build_message " << my_rank << endl;
    //info_size_arr[my_rank] = info_size;
    //hash_size_arr[my_rank] = hash_size;

    //if (my_rank != 0) {
    //    size_t sum_info_size, sum_hash_size, sum_index_size, sum_dict_size;
    //    MPI_Recv(&sum_info_size, 1, MPI_UNSIGNED_LONG, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(&sum_index_size, 1, MPI_UNSIGNED_LONG, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(&sum_hash_size, 1, MPI_UNSIGNED_LONG, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(&sum_dict_size, 1, MPI_UNSIGNED_LONG, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    char* sum_info_buffer = new char[sum_info_size];
    //    char* sum_hash_buffer = new char[sum_hash_size];
    //    char* sum_index_buffer = new char[sum_index_size];
    //    char* sum_dict_buffer = new char[sum_dict_size];
    //    MPI_Recv(sum_info_buffer, sum_info_size, MPI_CHAR, 0, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(sum_hash_buffer, sum_hash_size, MPI_CHAR, 0, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(sum_index_buffer, sum_index_size, MPI_CHAR, 0, 203, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(sum_dict_buffer, sum_dict_size, MPI_CHAR, 0, 204, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //    vector<KssdSketchInfo>().swap(sketches);
    //    format_sketches_index(sum_info_buffer, sum_info_size,
    //                          sum_hash_buffer, sum_hash_size,
    //                          sum_index_buffer, sum_index_size,
    //                          sum_dict_buffer, sum_dict_size,
    //                          folder_path, sketches, sketchByFile, threads);
    //}
    //else {
    //    sketchByFile = loadKssdSketches(folder_path, threads, sketches, info);
    //    transSketches(sketches, info, folder_path, threads);
    //    string sum_info_file = folder_path + '/' + "kssd.info.sketch";
    //    string sum_hash_file = folder_path + '/' + "kssd.hash.sketch";
    //    string sum_index_file = folder_path + '/' + "kssd.sketch.index";
    //    string sum_dict_file = folder_path + '/' + "kssd.sketch.dict";
    //    char* sum_info_buffer;
    //    char* sum_hash_buffer;
    //    char* sum_index_buffer;
    //    char* sum_dict_buffer;
    //    size_t sum_info_size, sum_hash_size, sum_index_size, sum_dict_size;
    //    build_message(sum_info_buffer, sum_info_size, sum_info_file);
    //    build_message(sum_index_buffer, sum_index_size, sum_index_file);
    //    build_message(sum_hash_buffer, sum_hash_size, sum_hash_file);
    //    build_message(sum_dict_buffer, sum_dict_size, sum_dict_file);

    //    for (int id = 1; id < comm_sz; id++) {
    //        MPI_Send(&sum_info_size, 1, MPI_UNSIGNED_LONG, id, 101, MPI_COMM_WORLD);
    //        MPI_Send(&sum_index_size, 1, MPI_UNSIGNED_LONG, id, 102, MPI_COMM_WORLD);
    //        MPI_Send(&sum_hash_size, 1, MPI_UNSIGNED_LONG, id, 103, MPI_COMM_WORLD);
    //        MPI_Send(&sum_dict_size, 1, MPI_UNSIGNED_LONG, id, 104, MPI_COMM_WORLD);
    //        MPI_Send(sum_info_buffer, sum_info_size, MPI_CHAR, id, 201, MPI_COMM_WORLD);
    //        MPI_Send(sum_hash_buffer, sum_hash_size, MPI_CHAR, id, 202, MPI_COMM_WORLD);
    //        MPI_Send(sum_index_buffer, sum_index_size, MPI_CHAR, id, 203, MPI_COMM_WORLD);
    //        MPI_Send(sum_dict_buffer, sum_dict_size, MPI_CHAR, id, 204, MPI_COMM_WORLD);
    //    }
    //}

    //cerr << "the sketch size of rank " << my_rank << " is: " << sketches.size() << endl;

    int start_index, end_index;
    if (my_rank != 0) {
      MPI_Send(&threads, 1, MPI_INT, 0, my_rank + comm_sz * 4, MPI_COMM_WORLD);
      MPI_Recv(&start_index, 1, MPI_UNSIGNED_LONG, 0, my_rank + comm_sz * 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&end_index, 1, MPI_UNSIGNED_LONG, 0, my_rank + comm_sz * 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
      int* thread_arr = new int[comm_sz];
      thread_arr[0] = threads;
      int total_threads = 0;
      total_threads += threads;
      for (int id = 1; id < comm_sz; id++) {
        MPI_Recv(&thread_arr[id], 1, MPI_INT, id, id + comm_sz * 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        total_threads += thread_arr[id];
      }
      cerr << "the total threads is: " << total_threads << endl;
      size_t N = sketches.size();
      size_t total_distance_time = N * (N - 1) / 2;
      size_t* task_arr = new size_t[comm_sz];
      for (int i = 0; i < comm_sz; i++) {
        task_arr[i] = total_distance_time * thread_arr[i] / total_threads;
      }
      size_t* start_index_arr = new size_t[comm_sz];
      size_t* end_index_arr = new size_t[comm_sz];
      int cur_start_index = 0;
      int cur_end_index = N - 1;
      size_t global_index = 0;
      for (int i = 0; i < comm_sz; i++) {
        start_index_arr[i] = global_index;
        size_t accumulate_time = 0;
        while (accumulate_time < task_arr[i] && global_index < N) {
          accumulate_time += global_index + 1;
          global_index++;
        }
        end_index_arr[i] = global_index;
      }
      end_index_arr[comm_sz - 1] = max(end_index_arr[comm_sz - 1], N);
      start_index = start_index_arr[0];
      end_index = end_index_arr[0];
      for (int i = 1; i < comm_sz; i++) {
        cerr << "start_index: " << start_index_arr[i] << '\t' << "end_index: " << end_index_arr[i] << endl;
        MPI_Send(&start_index_arr[i], 1, MPI_UNSIGNED_LONG, i, i + comm_sz * 5, MPI_COMM_WORLD);
        MPI_Send(&end_index_arr[i], 1, MPI_UNSIGNED_LONG, i, i + comm_sz * 6, MPI_COMM_WORLD);
      }
    }

    string tmp_str = "++++++++++++my rank: " + to_string(my_rank) + ", start_index: " + to_string(start_index) + ", end_index: " + to_string(end_index);
    cout << tmp_str << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    distribute_compute_clusters(my_rank, comm_sz, sketches, info,
        start_index, end_index,
        sketchByFile, outputFile, is_newick_tree, folder_path,
        threshold, isSave, threads, no_dense, isContainment);

    MPI_Finalize();
  }

