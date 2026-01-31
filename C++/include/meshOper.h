#pragma once
#include "class.h"
extern int globalTest;

// Create point-face adjacency table
void createPFneighbor(const MyMatrixXf& point, const MyMatrixXf& face, MyMatrixXf& PFneighbor);
// Calculate the angle distribution (in degrees) between normals of adjacent faces in the supporting region of a line segment
MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, MyMatrixXf& region_edge_inside);
MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, int sign = 0);
// Determine whether a point exists in the point set (by coordinates)
bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point, MyMatrixXf& exist_point_id);
bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point);
// Determine whether a line exists in the line set
bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line, MyMatrixXf& exist_line_id);
bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line);
// Determine whether a face exists in the face set
bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face, MyMatrixXf& exist_face_id);
bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face);
// Determine whether there are constrained edges in the supporting region of a line segment
bool haveConstrainedEdge(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& constrained_edge);
// Extract points and edges in the region, and mark them (boundary points and edges are marked as 0, internal as 1, if non-manifold structure appears, edge mark is usually greater than 1)
void getRegionPE(const MyMatrixXf& region_face, MyMatrixXf& region_point, MyMatrixXf& region_edge, int sign = 0);
// Calculate the area of the supporting region of a line segment (sum of all face areas)
float getRegionArea(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx);
// Query faces adjacent to the point set (input is point id)
MyMatrixXf getPointNeighFace(const MyMatrixXf& point, const MyMatrixXf& PFneighbor);
// Extract the boundary of the influence region
void getOutlineSort(const MyMatrixXf& region_face, MyMatrixXf& outline_point_sort, MyMatrixXf& outline_sort);
// Sort the boundary
MyMatrixXf regionEdgeSort(MyMatrixXf outline);
// Calculate the normal vector of each face
MyMatrixXf getFacesNormal(const MyMatrixXf& point, const MyMatrixXf& face);
MyMatrixXf getPerFaceNormal(const MyMatrixXf& tri);
// Find triangles adjacent to the current point (input is point id)
MyMatrixXf findPointNeighFace(const MyMatrixXf& face, float this_point);
// Find triangles adjacent to the current line segment
MyMatrixXf findEdgeNeighFace(const MyMatrixXf& face, const MyMatrixXf& this_edge);
// Find the other face adjacent to a certain edge of the current face
float findCertainNeighFace(const MyMatrixXf& face, const MyMatrixXf& face_idx, const MyMatrixXf& line, float this_face_idx);
// Modify the point indices in the face table based on a map
MyMatrixXf facePointUpdate(const MyMatrixXf& old_face, const MyMatrixXf& map);
// Divide the supporting region of the line segment into different regions according to connectivity
MyMatrixXf findConnectRegion(const MyMatrixXf& face, const MyMatrixXf& region_face_idx);
// Find the connected components of edges
MyMatrixXf findConnectEdge(const MyMatrixXf& edge);
// Get the edges of a triangle face
MyMatrixXf getFaceEdge(const MyMatrixXf& face);
