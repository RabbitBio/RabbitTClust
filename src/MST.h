#ifndef H_MST_GRAPH
#define H_MST_GRAPH

#include <iostream>
#include <vector>
#include <unordered_set>
#include "SketchInfo.h"
struct NeighborNode{
	int id;
	double distance;
	NeighborNode(int i, double d){
		id = i;
		distance = d;
	}
};

struct EdgeInfo{
	int preNode;
	int sufNode;
	double dist;
};


struct Graph{
	int node;
	int curNeighbor;
	std::vector<NeighborNode> neighbor;

};


struct MST{
	std::unordered_set<int> nodes;
	std::vector<EdgeInfo> edges;

};

bool cmpEdge(EdgeInfo e1, EdgeInfo e2);

bool cmpNeighbor(NeighborNode n1, NeighborNode n2);

std::vector<EdgeInfo> kruskalAlgorithm(std::vector<EdgeInfo>graph, int vertices);

vector<EdgeInfo> generateMST(vector<SketchInfo>& sketches, string sketchFunc, int threads);

vector<EdgeInfo> append_MST(vector<SketchInfo>& pre_sketches, vector<SketchInfo>& append_sketches, int sketch_func_id, int threads, int ** &denseArr, int denseSpan, uint64_t* &aniArr);

vector<EdgeInfo> modifyMST(vector<SketchInfo>& sketches, int start_index, int sketch_func_id, int threads, bool no_dense, int** &denseArr, int denseSpan, uint64_t* &aniArr);

vector<EdgeInfo> compute_kssd_mst(vector<KssdSketchInfo>& sketches, KssdParameters info, const string folder_path, int start_index, bool no_dense, bool isContainment, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr, double threshold, KssdInvertedIndex* inverted_index = nullptr);
vector<EdgeInfo> compute_minhash_mst(vector<SketchInfo>& sketches, int start_index, bool no_dense, bool isContainment, int threads, int** &denseArr, int denseSpan, uint64_t* &aniArr, double threshold, int kmerSize, MinHashInvertedIndex* inverted_index = nullptr);

std::vector<EdgeInfo> generateForest(std::vector <EdgeInfo> mst, double threshhold);

std::vector <std::vector<int> > generateCluster(std::vector<EdgeInfo> forest, int vertices);

vector<vector<int>> generateClusterWithBfs(vector<EdgeInfo> forest, int vertices);

vector<EdgeInfo> modifyForest(vector<EdgeInfo> forset, vector<int> noiseArr, int threads);

typedef pair<int, int> PairInt;
vector<int> getNoiseNode(vector<PairInt> densePairArr, int alpha);

struct LinkageRow{
	int c1;      // first cluster id
	int c2;      // second cluster id
	double dist; // distance at merge
	int size;    // size of new cluster
};

std::vector<LinkageRow> get_linkage_from_mst(int N, const std::vector<EdgeInfo>& mst);

string get_newick_tree(const vector<SketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file);
string get_kssd_newick_tree(const vector<KssdSketchInfo>& sketches, const vector<EdgeInfo>& mst, bool sketch_by_file);

// Automatic threshold selection based on MST edge length distribution
struct ThresholdCandidate{
	double threshold;      // Candidate threshold value
	double gap_score;      // Gap size (difference between adjacent edge lengths)
	int edge_index;        // Index of the edge in sorted MST
	double confidence;     // Confidence score (0-1)
	std::string level;     // Suggested taxonomic level (e.g., "species", "genus")
	double stability_score; // Overall stability score (0-1), computed when --stability is enabled
	double stability_split; // Split sensitivity (stability when threshold decreases), computed when --stability is enabled
	double stability_merge; // Merge sensitivity (stability when threshold increases), computed when --stability is enabled
	int cluster_count;      // Number of clusters at this threshold
	int near_edge_count;    // Number of MST edges near threshold used for stability evaluation
};

struct EdgeLengthStats{
	double min_dist;
	double max_dist;
	double median_dist;
	double q1_dist;        // First quartile
	double q3_dist;        // Third quartile
	double mean_dist;
	double std_dev;
	std::vector<double> sorted_distances;
};

// Analyze MST edge length distribution and find natural breakpoints
EdgeLengthStats analyzeEdgeLengthDistribution(const std::vector<EdgeInfo>& mst);

// Find candidate thresholds based on gap detection in edge length distribution
std::vector<ThresholdCandidate> findThresholdCandidates(const std::vector<EdgeInfo>& mst, 
                                                         int max_candidates = 5,
                                                         double min_gap_ratio = 0.1,
                                                         bool enable_stability = false,
                                                         int num_vertices = 0);

// Stability evaluation result
struct StabilityResult{
	double overall;        // Overall stability (average of split and merge)
	double split;          // Split sensitivity (when threshold decreases)
	double merge;          // Merge sensitivity (when threshold increases)
	int near_edge_count;   // Number of edges near threshold used for evaluation
};

// Compute stability score for a threshold by measuring edge flip rate under perturbations
// Returns split/merge sensitivity separately for better interpretability
StabilityResult computeThresholdStability(const std::vector<EdgeInfo>& mst, 
                                         double threshold, 
                                         int num_vertices,
                                         double epsilon = 0.01,
                                         int num_samples = 5,
                                         int min_near_edges = 100);

// Select optimal threshold from candidates (can be extended with taxonomic validation)
ThresholdCandidate selectOptimalThreshold(const std::vector<ThresholdCandidate>& candidates,
                                          const std::vector<EdgeInfo>& mst);

// Print threshold analysis results to file
void printThresholdAnalysis(const std::vector<EdgeInfo>& mst, 
                           const EdgeLengthStats& stats,
                           const std::vector<ThresholdCandidate>& candidates,
                           const ThresholdCandidate& optimal,
                           const std::string& output_file);

#endif
