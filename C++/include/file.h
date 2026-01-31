// File read/write declarations
#pragma once
#pragma once
#include "class.h"
#include "MyMath.h"


// Read mesh OBJ file
bool readMeshOBJ(MyMatrixXf &point,MyMatrixXf &face, const char* filename);
// Read line OBJ file
bool readLineOBJ(MyMatrixXf &line_set, const char* filename);
// Save mesh to OBJ file
bool saveMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf sign, const char* filename);
// Save line to OBJ file
bool saveLineOBJ(MyMatrixXf line_set, const char* filename);
// Read constrained line segment TXT file
bool readConstrainedEdgeTXT(MyMatrixXf& constrained_edge, const char* filename);
// Save local mesh to OBJ file
bool saveSubMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename);
bool saveSubMeshOBJ2(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename);