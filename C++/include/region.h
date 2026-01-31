#pragma once
#include "class.h"
#include <pcl/kdtree/kdtree_flann.h>

MyMatrixXf getKnnFace(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& PFneighbor, const MyMatrixXf& face, const MyMatrixXf& line, const MyMatrixXf& sign, int K, double interval);

MyMatrixXf getBaseRegion(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& PFneighbor, const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& line, const MyMatrixXf& sign, int K, double interval);

MyMatrixXf isIntersectant(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& knnsearch_face_idx, const MyMatrixXf& line, double interval);
MyMatrixXf isIntersectant2(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& knnsearch_face_idx, const MyMatrixXf& line, double interval, double dist_threshold);


MyMatrixXf holeFilling(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor, int max_loop = 30);
MyMatrixXf holeFilling(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor, const MyMatrixXf& knn_face_idx);
bool holeDetection(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor);

MyMatrixXf shapeRepairing(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor);

MyMatrixXf regionGrowth(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor, int k);

MyMatrixXf regionConnect(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor, int max_loop = 30);

MyMatrixXf regionExtend(const MyMatrixXf& face, const MyMatrixXf& knn_face_idx, const MyMatrixXf& region_face_idx, const MyMatrixXf& normal, const MyMatrixXf& dist, const MyMatrixXf& PFneighbor, float threshold_1, float threshold_2, int max_loop = 200);

MyMatrixXf regionSegment(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& base_region_face_idx, const MyMatrixXf& line, const MyMatrixXf& constrained_edge, float interval);