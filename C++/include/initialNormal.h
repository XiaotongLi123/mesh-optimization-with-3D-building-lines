#pragma once
#include "class.h"

//ransac计算初始法向量
MyMatrixXf ransacGetAngleDistribution(const MyMatrixXf& point, const MyMatrixXf& line);
//ransac计算法向中心位置
MyMatrixXf ransacFindCenter(const MyMatrixXf& data, float sigma);