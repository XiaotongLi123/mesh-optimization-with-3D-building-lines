//常用的数学函数的声明
#pragma once
#include "class.h"
#include <pcl/kdtree/kdtree_flann.h>

extern int globalTest;

//由线段两端点求解直线方程
MyMatrixXf lineEquation(const MyMatrixXf& point_1, const MyMatrixXf& point_2);
//求数组最大值，并返回最大值索引（如有相同只返回第一个的索引）（考虑改成模板）
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
//求数组最小值，并返回最小值索引（如有相同只返回第一个的索引）
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
//求向量的模长
float norm(const MyMatrixXf& vec);
//求两个向量的内积
float dot(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
//求两个向量的夹角
float getVecAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
//求两个直线的夹角，输入直线的方向向量
float getLineAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
//求两个向量的叉乘
MyMatrixXf crossProduct(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2);
//求某个索引两侧k邻近的索引
MyMatrixXf getNeighKidx(const MyMatrixXf& idx, int this_idx, int neigh_k);
//判断三维空间点是否在三角平面内
bool isInTri(const MyMatrixXf& triangle, const MyMatrixXf& point);
//计算三角面的面积
float getTriArea(const MyMatrixXf& triangle);
//创建k-d树
//kdtree<float, 3> createKdtree(MyMatrixXf points, MyMatrixXf sign);
pcl::KdTreeFLANN<pcl::PointXYZ> createKdtree(const MyMatrixXf& points);
//进行knn搜索
//int getNearestNeighbor(kdtree<float, 3> tree, MyMatrixXf search_point);
vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K, const MyMatrixXf& sign);
vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K);
//对数组的索引进行排序（降序,oper==0,升序,oper==1）
template <typename T>
vector<size_t> sort_indexes(const vector<T>& v,int oper) {
    // 初始化索引向量
    vector<size_t> idx(v.size());
    //使用iota对向量赋0-?的连续值
    iota(idx.begin(), idx.end(), 0);

    // 通过比较v的值对索引idx进行排序
    if (oper == 0)
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
    else
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
    return idx;
}
//求数组的平均值
template <typename T>
T getMean(T* data, int num) {
    T sum = 0;
    for (int i = 0; i < num; i++) {
        sum += data[i];
    }
    return sum / num;
}
//求数组的中值
template <typename T>
T getMedian(T* data, int num) {
    sort(data, data + num, greater<T>());
    if (num % 2 == 0) return (data[num / 2] + data[num / 2 - 1]) / 2;
    else return data[num / 2];
}
//求数组的q百分位数位置
template <typename T>
T getQuantile(T* data, int num, double q) {
    // q 为 [0, 1] 之间的小数，如 0.5 表示中位数，0.25 表示 25% 分位点
    if (num <= 0 || q < 0 || q > 1) return T();  // 边界检查

    // 升序排序（与 MATLAB 的 sort(data, 'ascend') 一致）
    std::sort(data, data + num);

    // 百分位位置（ceil，与 MATLAB 一致）
    int idx = static_cast<int>(std::ceil(num * q)) - 1;  // C++ 下标从 0 开始
    if (idx < 0) idx = 0;
    if (idx >= num) idx = num - 1;

    return data[idx];
}
//根据三点求平面方程的参数
MyMatrixXf getPlanePara(const MyMatrixXf& tri);
//判断某元素是否在矩阵内
bool anyIsmember(const MyMatrixXf& data_set, float this_data);
bool anyIsmember(const MyMatrixXf& data_set, const MyMatrixXf& target);
//取矩阵的某些行
MyMatrixXf getSubMat_Rows(const MyMatrixXf& mat, const MyMatrixXf& rows);
//取矩阵的某些列
MyMatrixXf getSubMat_Cols(const MyMatrixXf& mat, const MyMatrixXf& cols);

MyMatrixXf getSubMat(const MyMatrixXf& mat, const MyMatrixXf& rows, const MyMatrixXf& cols);
//给子矩阵赋值
void valueSubMat(MyMatrixXf& mat, const MyMatrixXf& rows, const MyMatrixXf& cols, float value);
//获得满足条件的向量索引
MyMatrixXf getIndex(const MyMatrixXf& vec, float value, int oper);
//矩阵转向量
vector<float> mat2vec(const MyMatrixXf& mat);
//向量转矩阵
template <class T>
MyMatrixXf vec2mat(vector<T> vec) {
    MyMatrixXf mat = MyMatrixXf::Zero(vec.size(), 1);
    for (int i = 0; i < vec.size(); i++) {
        mat(i, 0) = (float)vec[i];
    }
    return mat;
}
//矩阵翻转
MyMatrixXf matFlip(const MyMatrixXf& mat);
//求并集
MyMatrixXf MyUnion(const MyMatrixXf& A, const MyMatrixXf& B);
//求交集
MyMatrixXf MyIntersection(const MyMatrixXf& A, const MyMatrixXf& B);
//求差集
MyMatrixXf MyDifference(const MyMatrixXf& A, const MyMatrixXf& B);