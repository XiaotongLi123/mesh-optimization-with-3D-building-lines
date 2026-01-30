#pragma once
#include "class.h"
extern int globalTest;

//创建点面邻接表
void createPFneighbor(const MyMatrixXf& point, const MyMatrixXf& face, MyMatrixXf& PFneighbor);
//求线段支撑域相邻面的法向量夹角分布，以角度为单位
MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, MyMatrixXf& region_edge_inside);
MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, int sign = 0);
//判断点是否在点集中(点坐标)
bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point, MyMatrixXf& exist_point_id);
bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point);
//判断线是否在线集中
bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line, MyMatrixXf& exist_line_id);
bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line);
//判断面是否在面集中
bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face, MyMatrixXf& exist_face_id);
bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face);
//判断线段支撑域中是否存在约束线段
bool haveConstrainedEdge(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& constrained_edge);
//提取区域中的点和边，并对点和边进行标记（边界点和边标为0，内部点和边标为1，若出现非流形结构，边标记一般大于1）
void getRegionPE(const MyMatrixXf& region_face, MyMatrixXf& region_point, MyMatrixXf& region_edge, int sign = 0);
//求线段支撑域的面积（所有face的面积和）
float getRegionArea(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx);
//查询与点集相邻接的面(输入的是点的id)
MyMatrixXf getPointNeighFace(const MyMatrixXf& point, const MyMatrixXf& PFneighbor);
//提取影响域的边界
void getOutlineSort(const MyMatrixXf& region_face, MyMatrixXf& outline_point_sort, MyMatrixXf& outline_sort);
//将边界排序
MyMatrixXf regionEdgeSort(MyMatrixXf outline);
//求每个face的法向量
MyMatrixXf getFacesNormal(const MyMatrixXf& point, const MyMatrixXf& face);
MyMatrixXf getPerFaceNormal(const MyMatrixXf& tri);
//寻找与当前点相邻接的三角面（输入点的id）
MyMatrixXf findPointNeighFace(const MyMatrixXf& face, float this_point);
//寻找与当前线段相邻接的三角面
MyMatrixXf findEdgeNeighFace(const MyMatrixXf& face, const MyMatrixXf& this_edge);
//寻找与当前面某一边相邻的另一个面
float findCertainNeighFace(const MyMatrixXf& face, const MyMatrixXf& face_idx, const MyMatrixXf& line, float this_face_idx);
//基于map修改面表对应的点号
MyMatrixXf facePointUpdate(const MyMatrixXf& old_face, const MyMatrixXf& map);
//将线段支撑域按照是否连通划分为不同的区域
MyMatrixXf findConnectRegion(const MyMatrixXf& face, const MyMatrixXf& region_face_idx);
//寻找边的连通成分
MyMatrixXf findConnectEdge(const MyMatrixXf& edge);
//取三角面的边
MyMatrixXf getFaceEdge(const MyMatrixXf& face);
