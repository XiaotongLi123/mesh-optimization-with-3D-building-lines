//文件读写声明
#pragma once
#include "class.h"
#include "MyMath.h"


//读取网格obj文件
bool readMeshOBJ(MyMatrixXf &point,MyMatrixXf &face, const char* filename);
//读取线obj文件
bool readLineOBJ(MyMatrixXf &line_set, const char* filename);
//保存网格至obj文件
bool saveMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf sign, const char* filename);
//保存线至obj文件
bool saveLineOBJ(MyMatrixXf line_set, const char* filename);
//读取约束线段txt文件
bool readConstrainedEdgeTXT(MyMatrixXf& constrained_edge, const char* filename);
//保存局部网格至obj文件
bool saveSubMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename);
bool saveSubMeshOBJ2(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename);