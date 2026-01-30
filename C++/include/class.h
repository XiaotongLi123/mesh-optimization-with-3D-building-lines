//类定义
#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
//#include <pcl/point_cloud.h> // 包含了点云类等定义
//#include <pcl/kdtree/kdtree_flann.h> // 包含了kdtree等

using namespace std;
using Eigen::Vector3f;

typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MyMatrixXf;

struct pointProjection {
	float dist = 0.0;
	MyMatrixXf nearest_point = MyMatrixXf::Zero(1, 3);
};

