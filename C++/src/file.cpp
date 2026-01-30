//文件读写函数定义
//文件读写函数定义

#include "file.h"

bool readMeshOBJ(MyMatrixXf& point, MyMatrixXf& face, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    char attri;
    MyMatrixXf point_ori = MyMatrixXf::Zero(9999999, 3);
    MyMatrixXf face_ori = MyMatrixXf::Zero(9999999, 3);
    int pts_num = 0;
    int face_num = 0;
    char line[256];  // 定义一个字符数组来临时存储注释
    while (!feof(fp)) {
        fscanf(fp, "%c", &attri);
        if (attri == 'v') {
            fscanf(fp, "%f%f%f", &point_ori(pts_num, 0), &point_ori(pts_num, 1), &point_ori(pts_num, 2));
            pts_num += 1;
        }
        if (attri == 'f') {
            fscanf(fp, "%f%f%f", &face_ori(face_num, 0), &face_ori(face_num, 1), &face_ori(face_num, 2));
            face_num += 1;
        }
        if (attri == '#') {
            fgets(line, sizeof(line), fp);
        }
    }
    fclose(fp);
    point = point_ori.block(0, 0, pts_num, 3);
    face = face_ori.block(0, 0, face_num, 3);
    face = face - MyMatrixXf::Ones(face_num, 3);
    return true;
}

bool readLineOBJ(MyMatrixXf& line_set, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    char attri;
    MyMatrixXf line_set_ori = MyMatrixXf::Zero(999999, 3);
    int pts_num = 0;
    while (!feof(fp)) {
        fscanf(fp, "%c", &attri);
        if (attri == 'v') {
            fscanf(fp, "%f%f%f", &line_set_ori(pts_num, 0), &line_set_ori(pts_num, 1), &line_set_ori(pts_num, 2));
            pts_num += 1;
        }
    }
    fclose(fp);
    line_set = line_set_ori.block(0, 0, pts_num, 3);
    return true;
}

bool saveMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf sign, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    for (int i = 0; i < point.rows(); i++) {
        if (sign(i, 0) != 0)
            fprintf(fp, "v %f %f %f\n", point(i, 0), point(i, 1), point(i, 2));
        else
            fprintf(fp, "v 0 0 0\n");
    }
    face = face + MyMatrixXf::Ones(face.rows(), 3);
    for (int i = 0; i < face.rows(); i++) {
        fprintf(fp, "f %d %d %d\n", (int)face(i, 0), (int)face(i, 1), (int)face(i, 2));
    }
    fclose(fp);
    return true;
}

bool saveLineOBJ(MyMatrixXf line_set, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    for (int i = 0; i < line_set.rows(); i++) {
        fprintf(fp, "v %f %f %f\n", line_set(i, 0), line_set(i, 1), line_set(i, 2));
    }
    for (int i = 0; i < line_set.rows() / 2; i++) {
        fprintf(fp, "l %d %d\n", 2 * i + 1, 2 * i + 2);
    }
    fclose(fp);
    return true;
}

bool readConstrainedEdgeTXT(MyMatrixXf& constrained_edge, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    MyMatrixXf constrained_edge_ori = MyMatrixXf::Zero(999999, 2);
    int edge_num = 0;
    while (!feof(fp)) {
        fscanf(fp, "%f%f", &constrained_edge_ori(edge_num, 0), &constrained_edge_ori(edge_num, 1));
        edge_num += 1;
    }
    fclose(fp);
    constrained_edge = constrained_edge_ori.block(0, 0, edge_num, 2);
    constrained_edge = constrained_edge - MyMatrixXf::Ones(edge_num, 2);
    return true;
}

bool saveSubMeshOBJ(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    MyMatrixXf sub_point_sign = MyMatrixXf::Zero(point.rows(), 1);

    MyMatrixXf out_point_ori = MyMatrixXf::Zero(3 * region_face.rows(), 3);
    int out_point_num = 0;
    map<int, int> index_map;
    int point_idx;
    for (int i = 0; i < region_face.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            point_idx = region_face(i, j);
            if (sub_point_sign(point_idx, 0) == 0) {
                out_point_ori.row(out_point_num) = point.row(point_idx);
                index_map.insert(pair<int, int>(point_idx, out_point_num));
                out_point_num += 1;
                sub_point_sign(point_idx, 0) = 1;
            }
            region_face(i, j) = index_map[point_idx];
        }
    }
    MyMatrixXf out_point = out_point_ori.block(0, 0, out_point_num, 3);
    for (int i = 0; i < out_point.rows(); i++) {
        if (sign(i, 0) != 0)
            fprintf(fp, "v %f %f %f\n", out_point(i, 0), out_point(i, 1), out_point(i, 2));
        else
            fprintf(fp, "v 0 0 0\n");
    }
    region_face = region_face + MyMatrixXf::Ones(region_face.rows(), 3);
    for (int i = 0; i < region_face.rows(); i++) {
        fprintf(fp, "f %d %d %d\n", (int)region_face(i, 0), (int)region_face(i, 1), (int)region_face(i, 2));
    }
    fclose(fp);
    return true;
}


bool saveSubMeshOBJ2(MyMatrixXf point, MyMatrixXf face, MyMatrixXf region_face_idx, MyMatrixXf sign, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        cout << "文件打开异常！" << endl;
        return false;
    }
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);

    for (int i = 0; i < point.rows(); i++) {
        if (sign(i, 0) != 0)
            fprintf(fp, "v %f %f %f\n", point(i, 0), point(i, 1), point(i, 2));
        else
            fprintf(fp, "v 0 0 0\n");
    }
    region_face = region_face + MyMatrixXf::Ones(region_face.rows(), 3);
    for (int i = 0; i < region_face.rows(); i++) {
        fprintf(fp, "f %d %d %d\n", (int)region_face(i, 0), (int)region_face(i, 1), (int)region_face(i, 2));
    }
    fclose(fp);
    return true;
}