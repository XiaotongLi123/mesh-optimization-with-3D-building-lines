#pragma once
#include "class.h"

// Calculate initial normal vector using RANSAC
MyMatrixXf ransacGetAngleDistribution(const MyMatrixXf& point, const MyMatrixXf& line);
// Calculate the center position of the normal using RANSAC
MyMatrixXf ransacFindCenter(const MyMatrixXf& data, float sigma);