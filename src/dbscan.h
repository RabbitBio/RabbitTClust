#ifndef H_DBSCAN
#define H_DBSCAN

/**
 * @file dbscan.h
 * @brief DBSCAN (Density-Based Clustering) implementation for genome sequences
 * @author Auto-generated
 * @date 2026-02-02
 * 
 * This module provides DBSCAN clustering functionality for sketch-based genome sequences.
 * DBSCAN can automatically discover clusters and identify noise points (outliers).
 */

#include "SketchInfo.h"
#include <vector>
#include <string>

using std::vector;
using std::string;

/**
 * @brief DBSCAN clustering result structure
 */
struct DBSCANResult {
    vector<vector<int>> clusters;  // Clusters (each cluster is a vector of point indices)
    vector<int> noise;              // Noise points (outliers)
    int num_clusters;               // Number of clusters found
    int num_noise;                  // Number of noise points
};

/**
 * @brief DBSCAN clustering for KSSD sketches
 * 
 * @param sketches Input KSSD sketches
 * @param eps Epsilon (distance threshold) for neighborhood
 * @param minPts Minimum number of points required to form a cluster
 * @param kmer_size K-mer size for distance calculation
 * @param threads Number of threads for parallel computation
 * @param knn_k k-NN parameter: keep only k nearest neighbors per point (0 to disable, default 0)
 * @return DBSCANResult containing clusters and noise points
 */
DBSCANResult KssdDBSCAN(
    vector<KssdSketchInfo>& sketches,
    double eps,
    int minPts,
    int kmer_size,
    int threads,
    int knn_k = 0
);

/**
 * @brief DBSCAN clustering for MinHash sketches
 * 
 * @param sketches Input MinHash sketches
 * @param eps Epsilon (distance threshold) for neighborhood
 * @param minPts Minimum number of points required to form a cluster
 * @param sketch_func_id Sketch function ID (0 for MinHash)
 * @param threads Number of threads for parallel computation
 * @return DBSCANResult containing clusters and noise points
 */
DBSCANResult MinHashDBSCAN(
    vector<SketchInfo>& sketches,
    double eps,
    int minPts,
    int sketch_func_id,
    int threads
);

/**
 * @brief Print DBSCAN results to file
 * 
 * @param result DBSCAN clustering result
 * @param sketches Sketch information (for file names)
 * @param sketch_by_file Whether sketches are organized by file
 * @param output_file Output file path
 * @param eps Epsilon parameter used
 * @param minPts MinPts parameter used
 */
void printDBSCANResult(
    const DBSCANResult& result,
    const vector<SketchInfo>& sketches,
    bool sketch_by_file,
    const string& output_file,
    double eps,
    int minPts
);

/**
 * @brief Print DBSCAN results for KSSD sketches
 */
void printKssdDBSCANResult(
    const DBSCANResult& result,
    const vector<KssdSketchInfo>& sketches,
    bool sketch_by_file,
    const string& output_file,
    double eps,
    int minPts
);

#endif // H_DBSCAN

