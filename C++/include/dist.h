//计算几何元素之间的距离函数的声明
#pragma once
#include "class.h"

float point2lineDist(const MyMatrixXf& point, const MyMatrixXf& line);

float point2pointDist(const MyMatrixXf& point_1, const MyMatrixXf& point_2);

float point2faceDist(const MyMatrixXf& point, const MyMatrixXf& triangle, MyMatrixXf& nearest_point);
float point2faceDist(const MyMatrixXf& point, const MyMatrixXf& triangle);

float point2segmentDist(const MyMatrixXf& point, const MyMatrixXf& line, MyMatrixXf& nearest_point);
float point2segmentDist(const MyMatrixXf& point, const MyMatrixXf& line);

float point2meshDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& this_point, MyMatrixXf& nearest_face_idx);

float point2meshDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& this_point);

float line2meshMedianDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval);
float line2meshMeanDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval);
float line2meshMaxDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval);

float point2planeDist(const MyMatrixXf& point, const MyMatrixXf& plane);

float line2faceDist(const MyMatrixXf& tri, const MyMatrixXf& line, float interval);

float line2faceDist2(const MyMatrixXf& tri, const MyMatrixXf& line);
