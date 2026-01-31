//���õ���ѧ����������
#pragma once
#include "class.h"
#include <pcl/kdtree/kdtree_flann.h>

extern int globalTest;

// Get the equation of a line given two points
MyMatrixXf lineEquation(const MyMatrixXf& point_1, const MyMatrixXf& point_2);
// If there are multiple maximum (or minimum) values, only the first one is returned (the others are ignored, which is a bug)
template <typename T>
T getMax(T* data, int data_num, int& max_index) {
    max_index = 0;
    T max = data[0];
    for (int i = 0; i < data_num; i++) {
        if (max < data[i]) {
            max_index = i;
            max = data[i];
        }
    }
    return max;
}
template <typename T>
T getMax(T* data, int data_num) {
    T max = data[0];
    for (int i = 0; i < data_num; i++) {
        if (max < data[i]) {
            max = data[i];
        }
    }
    return max;
}
float getMax(const MyMatrixXf& data, MyMatrixXf& max_index);
// Get the minimum value and its index in the data array
template <typename T>
T getMin(T* data, int data_num, int& min_index) {
    min_index = 0;
    T min = data[0];
    for (int i = 0; i < data_num; i++) {
        if (min > data[i]) {
            min_index = i;
            min = data[i];
        }
    }
    return min;
}
template <typename T>
T getMin(T* data, int data_num) {
    T min = data[0];
    for (int i = 0; i < data_num; i++) {
        if (min > data[i]) {
            min = data[i];
        }
    }
    return min;
}
float getMin(const MyMatrixXf& data, MyMatrixXf& min_index);
// Calculate the norm of a vector
float norm(const MyMatrixXf& vec);
// Calculate the dot product of two vectors
float dot(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
// Calculate the angle between two vectors
float getVecAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
// Calculate the angle between two lines
float getLineAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
// Calculate the cross product of two vectors
MyMatrixXf crossProduct(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
// Get indices of elements in a vector that satisfy a certain condition
MyMatrixXf getNeighKidx(const MyMatrixXf& idx, int this_idx, int neigh_k);
// Determine whether a point is inside a triangle
bool isInTri(const MyMatrixXf& triangle, const MyMatrixXf& point);
// Calculate the area of a triangle
float getTriArea(const MyMatrixXf& triangle);
// Create pcl kdtree
pcl::KdTreeFLANN<pcl::PointXYZ> createKdtree(const MyMatrixXf& points);
// Get K nearest neighbors of a search point
vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K, const MyMatrixXf& sign);
vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K);
// Sort indices based on corresponding values in a vector
template <typename T>
vector<size_t> sort_indexes(const vector<T>& v,int oper) {
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    if (oper == 0)
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
    else
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
    return idx;
}
// Calculate the mean value
template <typename T>
T getMean(T* data, int num) {
    T sum = 0;
    for (int i = 0; i < num; i++) {
        sum += data[i];
    }
    return sum / num;
}
// Calculate the median value
template <typename T>
T getMedian(T* data, int num) {
    sort(data, data + num, greater<T>());
    if (num % 2 == 0) return (data[num / 2] + data[num / 2 - 1]) / 2;
    else return data[num / 2];
}
// Calculate the quantile value
template <typename T>
T getQuantile(T* data, int num, double q) {
    // q is in the range [0, 1], e.g., 0.5 represents the median, 0.25 represents the 25th percentile
    if (num <= 0 || q < 0 || q > 1) return T();  // boundary check
    std::sort(data, data + num);
    int idx = static_cast<int>(std::ceil(num * q)) - 1;
    if (idx < 0) idx = 0;
    if (idx >= num) idx = num - 1;

    return data[idx];
}
// Parameters of the least squares method
MyMatrixXf getPlanePara(const MyMatrixXf& tri);
// Determine whether a certain element exists in the set
bool anyIsmember(const MyMatrixXf& data_set, float this_data);
bool anyIsmember(const MyMatrixXf& data_set, const MyMatrixXf& target);
// Get certain rows of the matrix
MyMatrixXf getSubMat_Rows(const MyMatrixXf& mat, const MyMatrixXf& rows);
// Get certain columns of the matrix
MyMatrixXf getSubMat_Cols(const MyMatrixXf& mat, const MyMatrixXf& cols);

MyMatrixXf getSubMat(const MyMatrixXf& mat, const MyMatrixXf& rows, const MyMatrixXf& cols);
// Assign value to submatrix
void valueSubMat(MyMatrixXf& mat, const MyMatrixXf& rows, const MyMatrixXf& cols, float value);
// Get indices of elements in a vector that meet a certain condition
MyMatrixXf getIndex(const MyMatrixXf& vec, float value, int oper);
// Convert matrix to vector
vector<float> mat2vec(const MyMatrixXf& mat);
// Convert vector to matrix
template <class T>
MyMatrixXf vec2mat(vector<T> vec) {
    MyMatrixXf mat = MyMatrixXf::Zero(vec.size(), 1);
    for (int i = 0; i < vec.size(); i++) {
        mat(i, 0) = (float)vec[i];
    }
    return mat;
}
// Flip the matrix
MyMatrixXf matFlip(const MyMatrixXf& mat);
// Set operations
MyMatrixXf MyUnion(const MyMatrixXf& A, const MyMatrixXf& B);
MyMatrixXf MyIntersection(const MyMatrixXf& A, const MyMatrixXf& B);
MyMatrixXf MyDifference(const MyMatrixXf& A, const MyMatrixXf& B);