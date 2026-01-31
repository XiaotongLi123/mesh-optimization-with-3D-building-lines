// Common mathematical function definitions
#include "myMath.h"
#include "dist.h"


float getMin(const MyMatrixXf& data, MyMatrixXf& min_index) {
    float min = data(0, 0);
    min_index = MyMatrixXf::Zero(1, 2);
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            if (min > data(i, j)) {
                min = data(i, j);
                min_index(0, 0) = i;
                min_index(0, 1) = j;
            }
        }
    }
    return min;
}

float getMax(const MyMatrixXf& data, MyMatrixXf& max_index) {
    float max = data(0, 0);
    max_index = MyMatrixXf::Zero(1, 2);
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            if (max < data(i, j)) {
                max = data(i, j);
                max_index(0, 0) = i;
                max_index(0, 1) = j;
            }
        }
    }
    return max;
}

MyMatrixXf lineEquation(const MyMatrixXf& point_1, const MyMatrixXf& point_2) {
    MyMatrixXf line_para = MyMatrixXf::Zero(1, 3);
    float x1 = point_1(0, 0);
    float x2 = point_2(0, 0);
    float y1 = point_1(0, 1);
    float y2 = point_2(0, 1);
    float z1 = point_1(0, 2);
    float z2 = point_2(0, 2);
    if (x1 == x2 && y1 == y2 && z1 == z2) {
        return line_para;
    }
    line_para(0, 0) = (x2 - x1) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    line_para(0, 1) = (y2 - y1) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    line_para(0, 2) = (z2 - z1) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    return line_para;
}

float norm(const MyMatrixXf& vec) {
    return sqrt(vec(0, 0) * vec(0, 0) + vec(0, 1) * vec(0, 1) + vec(0, 2) * vec(0, 2));
}

float dot(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2) {
    return vec_1(0, 0) * vec_2(0, 0) + vec_1(0, 1) * vec_2(0, 1) + vec_1(0, 2) * vec_2(0, 2);
}

float getVecAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2) {
    float temp = dot(vec_1, vec_2) / (norm(vec_1) * norm(vec_2));
    if (temp >= 1) {
        return acos(1);
    }
    if (temp <= -1) {
        return acos(-1);
    }
    return acos(temp);
}

float getLineAngle(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2) {
    float angle = getVecAngle(vec_1, vec_2);
    if (angle > M_PI / 2)
        return M_PI - angle;
    else
        return angle;
}

MyMatrixXf crossProduct(const MyMatrixXf& vec_1, const MyMatrixXf& vec_2) {
    MyMatrixXf res_vec = MyMatrixXf::Zero(1, 3);
    res_vec(0, 0) = vec_1(0, 1) * vec_2(0, 2) - vec_2(0, 1) * vec_1(0, 2);
    res_vec(0, 1) = vec_1(0, 2) * vec_2(0, 0) - vec_2(0, 2) * vec_1(0, 0);
    res_vec(0, 2) = vec_1(0, 0) * vec_2(0, 1) - vec_2(0, 0) * vec_1(0, 1);
    return res_vec;
}

MyMatrixXf getNeighKidx(const MyMatrixXf& idx, int this_idx, int neigh_k) {
    int idx_num = idx.rows();
    MyMatrixXf temp = MyMatrixXf::Zero(3 * idx_num, 1);
    temp << idx,
        idx,
        idx;
    return temp.block(idx_num + this_idx - neigh_k, 0, 2 * neigh_k + 1, 1);
}

bool isInTri(const MyMatrixXf& triangle, const MyMatrixXf& point) {
    bool is = 0;
    MyMatrixXf L1 = MyMatrixXf::Zero(2, 3);
    MyMatrixXf L2 = MyMatrixXf::Zero(2, 3);
    MyMatrixXf L3 = MyMatrixXf::Zero(2, 3);
    L1 << triangle.block(1, 0, 1, 3),
        triangle.block(2, 0, 1, 3);
    L2 << triangle.block(2, 0, 1, 3),
        triangle.block(0, 0, 1, 3);
    L3 << triangle.block(0, 0, 1, 3),
        triangle.block(1, 0, 1, 3);

    // Compute distances from the point to the three edges
    float D1 = point2lineDist(point, L1);
    float D2 = point2lineDist(point, L2);
    float D3 = point2lineDist(point, L3);
    float D = point2lineDist(triangle.block(0, 0, 1, 3), L1);
    // Compute lengths of the three edges
    float length_1 = point2pointDist(triangle.block(1, 0, 1, 3), triangle.block(2, 0, 1, 3));
    float length_2 = point2pointDist(triangle.block(2, 0, 1, 3), triangle.block(0, 0, 1, 3));
    float length_3 = point2pointDist(triangle.block(0, 0, 1, 3), triangle.block(1, 0, 1, 3));
    // Compute areas
    float area = D * length_1 / 2;
    float area_1 = D1 * length_1 / 2;
    float area_2 = D2 * length_2 / 2;
    float area_3 = D3 * length_3 / 2;
    float error = area - (area_1 + area_2 + area_3);

    if (abs(error) < area * 1e-7) {
        is = 1;
    }
    return is;
}

float getTriArea(const MyMatrixXf& triangle) {
    MyMatrixXf V1 = triangle.block(0, 0, 1, 3);
    MyMatrixXf V2 = triangle.block(1, 0, 1, 3);
    MyMatrixXf V3 = triangle.block(2, 0, 1, 3);
    float a = point2pointDist(V2, V3);
    float b = point2pointDist(V1, V3);
    float c = point2pointDist(V1, V2);
    float p = (a + b + c) / 2;
    return sqrt(p * (p - a) * (p - b) * (p - c));
}


pcl::KdTreeFLANN<pcl::PointXYZ> createKdtree(const MyMatrixXf& points) {
    // [1] Create a point cloud pointer
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    // [2] Generate pointcloud data
    cloud->width = points.rows();
    cloud->height = 1;
    cloud->points.resize(cloud->width * cloud->height);

    for (std::size_t i = 0; i < cloud->size(); ++i)
    {
        (*cloud)[i].x = points(i, 0);
        (*cloud)[i].y = points(i, 1);
        (*cloud)[i].z = points(i, 2);
    }

    // [3] Create KD-Tree object
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    // [4] Pass data into KDTREE, i.e., set point cloud data as KD-Tree structure
    kdtree.setInputCloud(cloud);
    return kdtree;
}

vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K, const MyMatrixXf& sign) {
    // K-nearest neighbor search
    std::vector<int> pointIdxKNNSearch(2 * K);
    // Set search distance
    std::vector<float> pointKNNSquaredDistance(500);
    pcl::PointXYZ pcl_search_point;
    pcl_search_point.x = search_point(0, 0);
    pcl_search_point.y = search_point(0, 1);
    pcl_search_point.z = search_point(0, 2);
    int count = 0; // Statistics on the number of neighboring points searched so far
    vector<int> k_neigh(K, -1);
    int sum = 0; // Statistics on the number of valid neighboring points found so far
    int this_search_num = 2 * K;
    while (count < kdtree.getInputCloud()->points.size()) {
        kdtree.nearestKSearch(pcl_search_point, this_search_num, pointIdxKNNSearch, pointKNNSquaredDistance);
        for (int i = this_search_num - 2 * K; i < this_search_num; i++) {
            if (sign(pointIdxKNNSearch[i], 0) == 1) {
                k_neigh[sum++] = pointIdxKNNSearch[i];
            }
            if (sum == K)return k_neigh;
        }
        count += 2 * K;
        this_search_num += 2 * K;
        if (this_search_num > kdtree.getInputCloud()->points.size())
            this_search_num = kdtree.getInputCloud()->points.size();
    }
    return k_neigh;
}
vector<int> getKneighbor(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& search_point, int K) {
    // K-nearest neighbor search
    std::vector<int> pointIdxKNNSearch(K);
    // Set search distance
    std::vector<float> pointKNNSquaredDistance(100);
    pcl::PointXYZ pcl_search_point;
    pcl_search_point.x = search_point(0, 0);
    pcl_search_point.y = search_point(0, 1);
    pcl_search_point.z = search_point(0, 2);
    kdtree.nearestKSearch(pcl_search_point, K, pointIdxKNNSearch, pointKNNSquaredDistance);

    return pointIdxKNNSearch;
}

MyMatrixXf getPlanePara(const MyMatrixXf& tri) {
    MyMatrixXf S11 = MyMatrixXf::Zero(2, 2);
    MyMatrixXf S12 = MyMatrixXf::Zero(2, 2);
    MyMatrixXf S13 = MyMatrixXf::Zero(2, 2);
    S11 << tri(1, 1) - tri(0, 1), tri(1, 2) - tri(0, 2),
        tri(2, 1) - tri(0, 1), tri(2, 2) - tri(0, 2);
    S12 << tri(1, 0) - tri(0, 0), tri(1, 2) - tri(0, 2),
        tri(2, 0) - tri(0, 0), tri(2, 2) - tri(0, 2);
    S13 << tri(1, 0) - tri(0, 0), tri(1, 1) - tri(0, 1),
        tri(2, 0) - tri(0, 0), tri(2, 1) - tri(0, 1);
    MyMatrixXf plane = MyMatrixXf::Zero(1, 4);
    plane(0, 0) = S11.determinant();
    plane(0, 1) = -S12.determinant();
    plane(0, 2) = S13.determinant();
    plane(0, 3) = -plane(0, 0) * tri(0, 0) - plane(0, 1) * tri(0, 1) - plane(0, 2) * tri(0, 2);
    return plane;
}

bool anyIsmember(const MyMatrixXf& data_set, float this_data) {
    int row = data_set.rows();
    int col = data_set.cols();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (abs(data_set(i, j) - this_data) < 1e-5) return true;
        }
    }
    return false;
}
bool anyIsmember(const MyMatrixXf& data_set, const MyMatrixXf& target) {
    int row = data_set.rows();
    int col = data_set.cols();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            for (int m = 0; m < target.rows(); m++) {
                for (int n = 0; n < target.cols(); n++) {
                    if (abs(data_set(i, j) - target(m, n)) < 1e-10) return true;
                }
            }
        }
    }
    return false;
}

MyMatrixXf getSubMat_Rows(const MyMatrixXf& mat, const MyMatrixXf& rows) {
    int r;
    MyMatrixXf sub_mat;
    if (rows.cols() == 1) {
        sub_mat = MyMatrixXf::Zero(rows.rows(), mat.cols());
        for (int i = 0; i < rows.rows(); i++) {
            r = rows(i, 0);
            sub_mat.row(i) = mat.row(r);
        }
        return sub_mat;
    }
    if (rows.rows() == 1) {
        sub_mat = MyMatrixXf::Zero(rows.cols(), mat.cols());
        for (int i = 0; i < rows.cols(); i++) {
            r = rows(0, i);
            sub_mat.row(i) = mat.row(r);
        }
        return sub_mat;
    }
    cout << "rows is not a vector!" << endl;
    return sub_mat;
}

MyMatrixXf getSubMat_Cols(const MyMatrixXf& mat, const MyMatrixXf& cols) {
    MyMatrixXf sub_mat;
    int c;
    if (cols.cols() == 1) {
        sub_mat = MyMatrixXf::Zero(mat.rows(), cols.rows());
        for (int i = 0; i < cols.rows(); i++) {
            c = cols(i, 0);
            sub_mat.col(i) = mat.col(c);
        }
        return sub_mat;
    }
    if (cols.rows() == 1) {
        sub_mat = MyMatrixXf::Zero(mat.rows(), cols.cols());
        for (int i = 0; i < cols.cols(); i++) {
            c = cols(0, i);
            sub_mat.col(i) = mat.col(c);
        }
        return sub_mat;
    }
    cout << "cols is not a vector!" << endl;
    return sub_mat;
}

MyMatrixXf getSubMat(const MyMatrixXf& mat, const MyMatrixXf& rows, const MyMatrixXf& cols) {
    MyMatrixXf sub_mat_rows = getSubMat_Rows(mat, rows);
    MyMatrixXf sub_mat = getSubMat_Cols(sub_mat_rows, cols);
    return sub_mat;
}

void valueSubMat(MyMatrixXf &mat, const MyMatrixXf& rows, const MyMatrixXf& cols,float value) {
    int r, c;
    for (int i = 0; i < rows.rows(); i++) {
        for (int j = 0; j < cols.cols(); j++) {
            r = rows(i); c = cols(j);
            mat(r, c) = value;
        }
    }
}

MyMatrixXf getIndex(const MyMatrixXf& vec, float value, int oper) {
    MyMatrixXf index = MyMatrixXf::Zero(vec.rows(), vec.cols());
    int count = 0;
    if (vec.rows() == 1) {
        for (int i = 0; i < vec.cols(); i++) {
            switch (oper) {
                 // Equal
            case 0: {
                if (vec(i) == value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Not equal
            case 1: {
                if (vec(i) != value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Larger than
            case 2: {
                if (vec(i) > value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Smaller than
            case 3: {
                if (vec(i) < value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
            default: {
                cout << "Invalid operation!" << endl;
                return index;
            }
            }
        }
        return index.block(0, 0, 1, count);
    }
    if (vec.cols() == 1) {
        for (int i = 0; i < vec.rows(); i++) {
            switch (oper) {
                 // Equal
            case 0: {
                if (vec(i) == value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Not equal
            case 1: {
                if (vec(i) != value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Larger than
            case 2: {
                if (vec(i) > value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
                 // Smaller than
            case 3: {
                if (vec(i) < value) {
                    index(count) = i;
                    count += 1;
                }
                break;
            }
            default: {
                cout << "Invalid operation!" << endl;
                return index;
            }
            }
        }
        return index.block(0, 0, count, 1);
    }
    cout << "Input is not a vector!" << endl;
    return index;
}

vector<float> mat2vec(const MyMatrixXf& mat) {
    vector<float> vec(mat.rows() * mat.cols());
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            vec[i * mat.cols() + j] = mat(i, j);
        }
    }
    return vec;
}

MyMatrixXf matFlip(const MyMatrixXf& mat) {
    MyMatrixXf res_mat = MyMatrixXf::Zero(mat.rows(), mat.cols());
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            res_mat(i, j) = mat(mat.rows() - 1 - i, mat.cols() - 1 - j);
        }
    }
    return res_mat;
}

MyMatrixXf MyUnion(const MyMatrixXf& A, const MyMatrixXf& B) {
    vector<float> a, b, temp;
    a = mat2vec(A);
    b = mat2vec(B);
    temp.resize(a.size() + b.size());
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    auto it = set_union(a.begin(), a.end(), b.begin(), b.end(), temp.begin());
    temp.resize(it - temp.begin());
    MyMatrixXf res = vec2mat(temp);
    return res;
}

MyMatrixXf MyIntersection(const MyMatrixXf& A, const MyMatrixXf& B) {
    vector<float> a, b, temp;
    a = mat2vec(A);
    b = mat2vec(B);
    temp.resize(min(a.size(), b.size()));
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    auto it = set_intersection(a.begin(), a.end(), b.begin(), b.end(), temp.begin());
    temp.resize(it - temp.begin());
    MyMatrixXf res = vec2mat(temp);
    return res;
}

MyMatrixXf MyDifference(const MyMatrixXf& A, const MyMatrixXf& B) {
    vector<float> a, b, temp;
    a = mat2vec(A);
    b = mat2vec(B);
    temp.resize(a.size());
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    auto it = set_difference(a.begin(), a.end(), b.begin(), b.end(), temp.begin());
    temp.resize(it - temp.begin());
    MyMatrixXf res = vec2mat(temp);
    return res;
}

