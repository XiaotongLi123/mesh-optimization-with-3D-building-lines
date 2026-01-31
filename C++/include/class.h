// Class definitions
#pragma once
#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using Eigen::Vector3f;

typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MyMatrixXf;

struct pointProjection {
	float dist = 0.0;
	MyMatrixXf nearest_point = MyMatrixXf::Zero(1, 3);
};

