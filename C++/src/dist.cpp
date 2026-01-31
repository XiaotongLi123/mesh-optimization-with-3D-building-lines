// Definition of functions to calculate distances between geometric elements
#include "dist.h"
#include "myMath.h"
#include "lineProcess.h"
#include "omp.h"

float point2lineDist(const MyMatrixXf& point, const MyMatrixXf& line) {
    MyMatrixXf para = lineEquation(line.block(0, 0, 1, 3), line.block(1, 0, 1, 3));
    return sqrt((para(0, 1) * (point(0, 2) - line(0, 2)) - para(0, 2) * (point(0, 1) - line(0, 1))) * (para(0, 1) * (point(0, 2) - line(0, 2)) - para(0, 2) * (point(0, 1) - line(0, 1)))
        + (para(0, 0) * (point(0, 2) - line(0, 2)) - para(0, 2) * (point(0, 0) - line(0, 0))) * (para(0, 0) * (point(0, 2) - line(0, 2)) - para(0, 2) * (point(0, 0) - line(0, 0)))
        + (para(0, 0) * (point(0, 1) - line(0, 1)) - para(0, 1) * (point(0, 0) - line(0, 0))) * (para(0, 0) * (point(0, 1) - line(0, 1)) - para(0, 1) * (point(0, 0) - line(0, 0))))
        / norm(para);
}

float point2pointDist(const MyMatrixXf& point_1, const MyMatrixXf& point_2) {
    return sqrt((point_1(0, 0) - point_2(0, 0)) * (point_1(0, 0) - point_2(0, 0)) + (point_1(0, 1) - point_2(0, 1)) * (point_1(0, 1) - point_2(0, 1)) + (point_1(0, 2) - point_2(0, 2)) * (point_1(0, 2) - point_2(0, 2)));
}

float point2faceDist(const MyMatrixXf& point, const MyMatrixXf& tri,MyMatrixXf& nearest_point) {
    MyMatrixXf s11(2, 2);
    MyMatrixXf s12(2, 2);
    MyMatrixXf s13(2, 2);
    s11 << tri(1, 1) - tri(0, 1), tri(1, 2) - tri(0, 2),
        tri(2, 1) - tri(0, 1), tri(2, 2) - tri(0, 2);
    s12 << tri(1, 0) - tri(0, 0), tri(1, 2) - tri(0, 2),
        tri(2, 0) - tri(0, 0), tri(2, 2) - tri(0, 2);
    s13 << tri(1, 0) - tri(0, 0), tri(1, 1) - tri(0, 1),
        tri(2, 0) - tri(0, 0), tri(2, 1) - tri(0, 1);

    float a_t = s11.determinant();
    float b_t = -s12.determinant();
    float c_t = s13.determinant();

    float t_0 = -((point(0, 0) - tri(0, 0)) * a_t + (point(0, 1) - tri(0, 1)) * b_t + (point(0, 2) - tri(0, 2)) * c_t) / (a_t * a_t + b_t * b_t + c_t * c_t);

    MyMatrixXf point0 = MyMatrixXf::Zero(1, 3);
    point0(0, 0) = point(0, 0) + a_t * t_0;
    point0(0, 1) = point(0, 1) + b_t * t_0;
    point0(0, 2) = point(0, 2) + c_t * t_0;
    
    if (isInTri(tri,point0)) {
        nearest_point = point0;
        return abs((a_t * point(0, 0) + b_t * point(0, 1) + c_t * point(0, 2) - a_t * tri(0, 0) - b_t * tri(0, 1) - c_t * tri(0, 2)) / sqrt(a_t * a_t + b_t * b_t + c_t * c_t));
    }
    else {
        MyMatrixXf L1 = MyMatrixXf::Zero(2, 3);
        MyMatrixXf L2 = MyMatrixXf::Zero(2, 3);
        MyMatrixXf L3 = MyMatrixXf::Zero(2, 3);
        L1 << tri.row(1),
            tri.row(2);
        L2 << tri.row(2),
            tri.row(0);
        L3 << tri.row(0),
            tri.row(1);
        MyMatrixXf nearest_edge_point = MyMatrixXf::Zero(3, 3);
        MyMatrixXf nearest_1, nearest_2, nearest_3;
        float dist[3] = { point2segmentDist(point, L1,nearest_1) ,point2segmentDist(point, L2,nearest_2) ,point2segmentDist(point, L3,nearest_3) };
        nearest_edge_point << nearest_1,
            nearest_2,
            nearest_3;
        int nearest_id = 0;
        float p2f_dist = getMin(dist, 3, nearest_id);
        nearest_point = nearest_edge_point.row(nearest_id);
        return p2f_dist;
    }
}
float point2faceDist(const MyMatrixXf& point, const MyMatrixXf& triangle) {
    MyMatrixXf V1 = triangle.block(0, 0, 1, 3);
    MyMatrixXf V2 = triangle.block(1, 0, 1, 3);
    MyMatrixXf V3 = triangle.block(2, 0, 1, 3);

    MyMatrixXf s11(2, 2);
    MyMatrixXf s12(2, 2);
    MyMatrixXf s13(2, 2);
    s11 << V2(0, 1) - V1(0, 1), V2(0, 2) - V1(0, 2),
        V3(0, 1) - V1(0, 1), V3(0, 2) - V1(0, 2);
    s12 << V2(0, 0) - V1(0, 0), V2(0, 2) - V1(0, 2),
        V3(0, 0) - V1(0, 0), V3(0, 2) - V1(0, 2);
    s13 << V2(0, 0) - V1(0, 0), V2(0, 1) - V1(0, 1),
        V3(0, 0) - V1(0, 0), V3(0, 1) - V1(0, 1);

    float a_t = s11.determinant();
    float b_t = -s12.determinant();
    float c_t = s13.determinant();

    float t_0 = -((point(0, 0) - V1(0, 0)) * a_t + (point(0, 1) - V1(0, 1)) * b_t + (point(0, 2) - V1(0, 2)) * c_t) / (a_t * a_t + b_t * b_t + c_t * c_t);

    MyMatrixXf point0 = MyMatrixXf::Zero(1, 3);
    point0(0, 0) = point(0, 0) + a_t * t_0;
    point0(0, 1) = point(0, 1) + b_t * t_0;
    point0(0, 2) = point(0, 2) + c_t * t_0;

    if (isInTri(triangle, point0)) {
        return abs((a_t * point(0, 0) + b_t * point(0, 1) + c_t * point(0, 2) - a_t * triangle(0, 0) - b_t * triangle(0, 1) - c_t * triangle(0, 2)) / sqrt(a_t * a_t + b_t * b_t + c_t * c_t));
    }
    else {
        MyMatrixXf L1 = MyMatrixXf::Zero(2, 3);
        MyMatrixXf L2 = MyMatrixXf::Zero(2, 3);
        MyMatrixXf L3 = MyMatrixXf::Zero(2, 3);
        L1 << V2,
            V3;
        L2 << V3,
            V1;
        L3 << V1,
            V2;
        float dist[3] = { point2segmentDist(point, L1) ,point2segmentDist(point, L2) ,point2segmentDist(point, L3) };

        int nearest_id = 0;
        float p2f_dist = getMin(dist, 3, nearest_id);
        return p2f_dist;
    }
}

float point2segmentDist(const MyMatrixXf& point, const MyMatrixXf& line, MyMatrixXf& nearest_point) {
    float dist = 0.0;
    nearest_point = MyMatrixXf::Zero(1, 3);
    MyMatrixXf V1 = line.block(0, 0, 1, 3);
    MyMatrixXf V2 = line.block(1, 0, 1, 3);
    MyMatrixXf line_para = lineEquation(V1, V2);
    float d = -(line_para(0, 0) * point(0, 0) + line_para(0, 1) * point(0, 1) + line_para(0, 2) * point(0, 2));
    float t = -(line_para(0, 0) * V1(0, 0) + line_para(0, 1) * V1(0, 1) + line_para(0, 2) * V1(0, 2) + d) / pow(norm(line_para), 2);
    float x = V1(0, 0) + line_para(0, 0) * t;
    float y = V1(0, 1) + line_para(0, 1) * t;
    float z = V1(0, 2) + line_para(0, 2) * t;

    Vector3f vec_1(point(0, 0) - V1(0, 0), point(0, 1) - V1(0, 1), point(0, 2) - V1(0, 2));
    Vector3f vec_2(x - point(0, 0), y - point(0, 1), z - point(0, 2));
    Vector3f vec_4(V2(0, 0) - V1(0, 0), V2(0, 1) - V1(0, 1), V2(0, 2) - V1(0, 2));
    Vector3f vec_5(V2(0, 0) - point(0, 0), V2(0, 1) - point(0, 1), V2(0, 2) - point(0, 2));
    Vector3f vec_3 = vec_1 + vec_2;

    int m = -1;
    for (int i = 0; i < 3; i++) {
        if (vec_4[i] != 0) {
            m = i;
            break;
        }
    }
    if (m == -1) {
        cout << "The two endpoints of the line coincide!" << endl;
        return dist;
    }
    float lambda = vec_3[m] / vec_4[m];

    if (lambda >= 0 && lambda <= 1) {
        dist = vec_2.norm();
        nearest_point << x, y, z;
    }
    if (lambda < 0) {
        dist = vec_1.norm();
        nearest_point = V1;
    }
    if (lambda > 1) {
        dist = vec_5.norm();
        nearest_point = V2;
    }
    return dist;
}
float point2segmentDist(const MyMatrixXf& point, const MyMatrixXf& line) {
    float dist = 0.0;
    MyMatrixXf V1 = line.block(0, 0, 1, 3);
    MyMatrixXf V2 = line.block(1, 0, 1, 3);
    MyMatrixXf line_para = lineEquation(V1, V2);
    float d = -(line_para(0, 0) * point(0, 0) + line_para(0, 1) * point(0, 1) + line_para(0, 2) * point(0, 2));
    float t = -(line_para(0, 0) * V1(0, 0) + line_para(0, 1) * V1(0, 1) + line_para(0, 2) * V1(0, 2) + d) / pow(norm(line_para), 2);
    float x = V1(0, 0) + line_para(0, 0) * t;
    float y = V1(0, 1) + line_para(0, 1) * t;
    float z = V1(0, 2) + line_para(0, 2) * t;

    Vector3f vec_1(point(0, 0) - V1(0, 0), point(0, 1) - V1(0, 1), point(0, 2) - V1(0, 2));
    Vector3f vec_2(x - point(0, 0), y - point(0, 1), z - point(0, 2));
    Vector3f vec_4(V2(0, 0) - V1(0, 0), V2(0, 1) - V1(0, 1), V2(0, 2) - V1(0, 2));
    Vector3f vec_5(V2(0, 0) - point(0, 0), V2(0, 1) - point(0, 1), V2(0, 2) - point(0, 2));
    Vector3f vec_3 = vec_1 + vec_2;

    int m = -1;
    for (int i = 0; i < 3; i++) {
        if (vec_4[i] != 0) {
            m = i;
            break;
        }
    }
    if (m == -1) {
        cout << "The two endpoints of the line coincide!" << endl;
        return dist;
    }
    float lambda = vec_3[m] / vec_4[m];

    if (lambda >= 0 && lambda <= 1) {
        dist = vec_2.norm();
    }
    if (lambda < 0) {
        dist = vec_1.norm();
    }
    if (lambda > 1) {
        dist = vec_5.norm();
    }
    return dist;
}

float point2meshDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& this_point, MyMatrixXf& nearest_face_idx) {
    int region_face_num = region_face_idx.rows();
    MyMatrixXf p2f_dist = MyMatrixXf::Zero(region_face_num, 1);
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    int r;
    for (int i = 0; i < region_face_num; i++) {
        r = region_face_idx(i, 0);
        for (int j = 0; j < 3; j++) {
            tri.row(j) << point.row(face(r, j));
        }
        p2f_dist(i, 0) = point2faceDist(this_point, tri);
    }
    float p2m_dist = p2f_dist.minCoeff();
    MyMatrixXf nearest_ori = MyMatrixXf::Zero(999, 1);
    int count = 0;
    for (int i = 0; i < region_face_num; i++) {
        if (abs(p2f_dist(i, 0) - p2m_dist) < 1e-7) {
            nearest_ori(count, 0) = i;
            count += 1;
        }
    }
    nearest_face_idx = nearest_ori.block(0, 0, count, 1);
    return p2m_dist;
}

float point2meshDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& this_point) {
    int region_face_num = region_face_idx.rows();
    MyMatrixXf p2f_dist = MyMatrixXf::Zero(region_face_num, 1);
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    int r;
    for (int i = 0; i < region_face_num; i++) {
        r = region_face_idx(i, 0);
        for (int j = 0; j < 3; j++) {
            tri.row(j) << point.row(face(r, j));
        }
        p2f_dist(i, 0) = point2faceDist(this_point, tri);
    }
    return p2f_dist.minCoeff();
}

float line2meshMedianDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval) {
    MyMatrixXf node = getNode(line, interval);
    int node_num = node.rows();

    float* p2m_dist = new float[node_num];
#pragma omp parallel for
    for (int i = 0; i < node_num; i++) {
        p2m_dist[i] = point2meshDist(point, face, region_face_idx, node.row(i));
    }
    float l2m_dist = getMedian(p2m_dist, node_num);
    delete[] p2m_dist;
    return l2m_dist;
}

float line2meshMeanDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval) {
    MyMatrixXf node = getNode(line, interval);
    int node_num = node.rows();

    float* p2m_dist = new float[node_num];
#pragma omp parallel for
    for (int i = 0; i < node_num; i++) {
        p2m_dist[i] = point2meshDist(point, face, region_face_idx, node.row(i));
    }
    float l2m_dist = getMean(p2m_dist, node_num);
    delete[] p2m_dist;
    return l2m_dist;
}

float line2meshMaxDist(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& line, float interval) {
    MyMatrixXf node = getNode(line, interval);
    int node_num = node.rows();

    MyMatrixXf p2m_dist = MyMatrixXf::Zero(node_num, 1);
#pragma omp parallel for
    for (int i = 0; i < node_num; i++) {
        p2m_dist(i, 0) = point2meshDist(point, face, region_face_idx, node.row(i));
    }
    float l2m_dist = p2m_dist.maxCoeff();
    return l2m_dist;
}

float point2planeDist(const MyMatrixXf& point, const MyMatrixXf& plane) {
    return abs(plane(0, 0) * point(0, 0) + plane(0, 1) * point(0, 1) + plane(0, 2) * point(0, 2) + plane(0, 3)) / norm(plane.block(0, 0, 1, 3));
}

float line2faceDist(const MyMatrixXf& tri, const MyMatrixXf& line, float interval) {
    MyMatrixXf node = getNode(line, interval);
    int node_num = node.rows();

    float* p2m_dist = new float[node_num];
#pragma omp parallel for
    for (int i = 0; i < node_num; i++) {
        p2m_dist[i] = point2faceDist(node.row(i), tri);
    }
    float l2f_dist = getMin(p2m_dist, node_num);;
    delete[] p2m_dist;
    return l2f_dist;
}

float line2faceDist2(const MyMatrixXf& tri, const MyMatrixXf& line) {
    MyMatrixXf plane = getPlanePara(tri);
    float p2p_dist_1 = point2planeDist(line.row(0), plane);
    float p2p_dist_2 = point2planeDist(line.row(1), plane);
    return max({ p2p_dist_1,p2p_dist_2 });
}
