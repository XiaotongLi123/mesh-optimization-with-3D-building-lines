#pragma once
#include "class.h"
#include "file.h"

struct faceSet {
	MyMatrixXf face_set;
};

extern int globalTest;
extern int dp_count;
extern clock_t myTime;

bool regionRemeshing(const MyMatrixXf& point, MyMatrixXf& face, MyMatrixXf& sign, const MyMatrixXf& region_face_idx, const MyMatrixXf& node_idx, const MyMatrixXf& PFneighbor, MyMatrixXf& constrained_line, float w_1 = 0.5, float w_2 = 0.5, int line_sign = 0);

MyMatrixXf DP_remeshing(const MyMatrixXf& point, const MyMatrixXf& node, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, const float weight_1, const float weight_2, float& opt_val);

MyMatrixXf DP_polygonRemeshing(const MyMatrixXf& point, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, float weight_1, float weight_2, float& opt_val, MyMatrixXf& distribution);
