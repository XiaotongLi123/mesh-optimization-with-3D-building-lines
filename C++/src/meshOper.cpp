#include "meshOper.h"
#include "myMath.h"
#include "dist.h"
#include "omp.h"

void createPFneighbor(const MyMatrixXf& point, const MyMatrixXf& face, MyMatrixXf& PFneighbor) {
    //构建点面邻接表
    int pts_num = point.rows();
    int face_num = face.rows();
    int r;
    PFneighbor = MyMatrixXf::Zero(pts_num, 30);//第一列存储邻接面数量，不超过29
    for (int i = 0; i < face_num; i++) {
        for (int j = 0; j < 3; j++) {
            r = face(i, j);
            int this_PFneigh_num = PFneighbor(r, 0);//该点当前的面邻接数量
            if (this_PFneigh_num >= 29) {
                continue;
            }
            PFneighbor(r, this_PFneigh_num + 1) = i;
            PFneighbor(r, 0) += 1;
        }
    }
}

MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, MyMatrixXf& region_edge_inside) {
    //求线段支撑域相邻面的法向量夹角分布，以角度为单位
    int region_face_num = region_face_idx.rows();
    MyMatrixXf region_face = MyMatrixXf::Zero(region_face_num, 3);
    for (int i = 0; i < region_face_num; i++) {
        region_face.row(i) = face.row(region_face_idx(i, 0));
    }
    MyMatrixXf normal = getFacesNormal(point, region_face);
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);
    int inside_edge_num = 0;
    MyMatrixXf region_edge_inside_ori = MyMatrixXf::Zero(region_edge.rows(), 2);
    for (int i = 0; i < region_edge.rows(); i++) {
        if (region_edge(i, 2) == 1) {
            region_edge_inside_ori.row(inside_edge_num) = region_edge.block(i, 0, 1, 2);
            inside_edge_num += 1;
        }
    }
    region_edge_inside = region_edge_inside_ori.block(0, 0, inside_edge_num, 2);
    MyMatrixXf normal_angle = MyMatrixXf::Zero(inside_edge_num, 1);
    MyMatrixXf neigh_face_pair = MyMatrixXf::Zero(2, 1);
    for (int i = 0; i < inside_edge_num; i++) {
        neigh_face_pair = findEdgeNeighFace(region_face, region_edge_inside.row(i));
        normal_angle(i, 0) = getVecAngle(normal.row(neigh_face_pair(0, 0)), normal.row(neigh_face_pair(1, 0))) * 180 / M_PI;
    }
    return normal_angle;
}
MyMatrixXf regionFaceNormalAngle(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx,int sign) {
    //求线段支撑域相邻面的法向量夹角分布，以角度为单位
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    MyMatrixXf normal = getFacesNormal(point, region_face);
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge, sign);
    int inside_edge_num = 0;
    MyMatrixXf region_edge_inside_ori = MyMatrixXf::Zero(region_edge.rows(), 2);
    for (int i = 0; i < region_edge.rows(); i++) {
        if (region_edge(i, 2) == 1) {
            region_edge_inside_ori.row(inside_edge_num) = region_edge.block(i, 0, 1, 2);
            inside_edge_num += 1;
        }
    }
    MyMatrixXf region_edge_inside = region_edge_inside_ori.block(0, 0, inside_edge_num, 2);
    MyMatrixXf normal_angle = MyMatrixXf::Zero(inside_edge_num, 1);
    MyMatrixXf neigh_face_pair = MyMatrixXf::Zero(2, 1);
    for (int i = 0; i < inside_edge_num; i++) {
        neigh_face_pair = findEdgeNeighFace(region_face, region_edge_inside.row(i));
        normal_angle(i, 0) = getVecAngle(normal.row(neigh_face_pair(0, 0)), normal.row(neigh_face_pair(1, 0))) * 180 / M_PI;
    }
    return normal_angle;
}

bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point, MyMatrixXf& exist_point_id) {
    int pts_num = point_set.rows();
    bool is_exist_point = false;
    MyMatrixXf exist_point_id_ori = MyMatrixXf::Zero(pts_num, 1);
    int count = 0;
    for (int i = 0; i < pts_num; i++) {
        if (point2pointDist(point_set.row(i), this_point) < 1e-3) {
            is_exist_point = true;
            exist_point_id_ori(count, 0) = i;
            count += 1;
        }
    }
    exist_point_id = exist_point_id_ori.block(0, 0, count, 1);
    return is_exist_point;
}
bool isExistPoint(const MyMatrixXf& point_set, const MyMatrixXf& this_point) {
    int pts_num = point_set.rows();
    bool is_exist_point = false;
    for (int i = 0; i < pts_num; i++) {
        if (point2pointDist(point_set.row(i), this_point) < 1e-3) {
            is_exist_point = true;
            break;
        }
    }
    return is_exist_point;
}
bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line, MyMatrixXf& exist_line_id) {
    int line_num = line_set.rows();
    bool is_exist_line = false;
    MyMatrixXf exist_line_id_ori = MyMatrixXf::Zero(line_num, 1);
    int count = 0;
    for (int i = 0; i < line_num; i++) {
        if (anyIsmember(line_set.row(i), this_line(0, 0)) && anyIsmember(line_set.row(i), this_line(0, 1))) {
            is_exist_line = true;
            exist_line_id_ori(count, 0) = i;
            count += 1;
        }
    }
    exist_line_id = exist_line_id_ori.block(0, 0, count, 1);
    return is_exist_line;
}

bool isExistLine(const MyMatrixXf& line_set, const MyMatrixXf& this_line) {
    int line_num = line_set.rows();
    bool is_exist_line = false;
    for (int i = 0; i < line_num; i++) {
        if (anyIsmember(line_set.row(i), this_line(0, 0)) && anyIsmember(line_set.row(i), this_line(0, 1))) {
            is_exist_line = true;
            break;
        }
    }
    return is_exist_line;
}

bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face, MyMatrixXf& exist_face_id) {
    int face_num = face_set.rows();
    bool is_exist_face = false;
    MyMatrixXf exist_face_id_ori = MyMatrixXf::Zero(face_num, 1);
    int count = 0;
    for (int i = 0; i < face_num; i++) {
        if (anyIsmember(face_set.row(i), this_face(0, 0)) && anyIsmember(face_set.row(i), this_face(0, 1)) && anyIsmember(face_set.row(i), this_face(0, 2))) {
            is_exist_face = true;
            exist_face_id_ori(count, 0) = i;
            count += 1;
        }
    }
    exist_face_id = exist_face_id_ori.block(0, 0, count, 1);
    return is_exist_face;
}
bool isExistFace(const MyMatrixXf& face_set, const MyMatrixXf& this_face) {
    int face_num = face_set.rows();
    bool is_exist_face = false;
    for (int i = 0; i < face_num; i++) {
        if (anyIsmember(face_set.row(i), this_face(0, 0)) && anyIsmember(face_set.row(i), this_face(0, 1)) && anyIsmember(face_set.row(i), this_face(0, 2))) {
            is_exist_face = true;
        }
    }
    return is_exist_face;
}

bool haveConstrainedEdge(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& constrained_edge) {
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);

    int constrained_edge_num = constrained_edge.rows();
    if (constrained_edge_num == 0)return false;
    MyMatrixXf index, in_region_edge_ori, in_region_edge;
    for (int i = 0; i < constrained_edge_num; i++) {
        index = getIndex(region_edge.col(2), 1, 0);
        in_region_edge_ori = getSubMat_Rows(region_edge, index);
        in_region_edge = in_region_edge_ori.block(0, 0, in_region_edge_ori.rows(), 2);
        if (isExistLine(in_region_edge, constrained_edge.row(i))) return true;
    }
    return false;
}

void getRegionPE(const MyMatrixXf& region_face, MyMatrixXf& region_point, MyMatrixXf& region_edge,int sign) {
    //提取区域中的点和边，并对点和边进行标记（边界点和边标为0，内部点和边标为1，若出现非流形结构，边标记一般大于1）
    //region_point为2列矩阵，region_face为3列矩阵
    int region_face_num = region_face.rows();
    MyMatrixXf region_point_ori = MyMatrixXf::Zero(3 * region_face_num, 2);
    region_point_ori.col(1) = MyMatrixXf::Ones(3 * region_face_num, 1);
    MyMatrixXf region_edge_ori = MyMatrixXf::Zero(3 * region_face_num, 3);
    int region_edge_num = 0;
    int region_point_num = 0;
    int r;

    //提取区域中的边
    MyMatrixXf this_face = MyMatrixXf::Zero(1, 3);
    MyMatrixXf this_edge = MyMatrixXf::Zero(1, 2);
    MyMatrixXf exist_line_id;
    MyMatrixXf cols = MyMatrixXf::Zero(1, 2);
    cols << 0, 1;
    for (int i = 0; i < region_face_num; i++) {
        this_face = region_face.row(i);
        for (int j = 0; j < 3; j++) {
            if (j == 0) {
                this_edge(0, 0) = this_face(0, 0);
                this_edge(0, 1) = this_face(0, 1); 
            }
            if (j == 1) {
                this_edge(0, 0) = this_face(0, 1);
                this_edge(0, 1) = this_face(0, 2);
            }
            if (j == 2) {
                this_edge(0, 0) = this_face(0, 2);
                this_edge(0, 1) = this_face(0, 0);
            }
            if (isExistLine(region_edge_ori.block(0, 0, region_edge_num, 2), this_edge, exist_line_id)) {
                for (int k = 0; k < exist_line_id.rows(); k++) {
                    r = exist_line_id(k, 0);
                    region_edge_ori(r, 2) += 1;
                }
            }
            else {
                region_edge_ori(region_edge_num, 0) = this_edge(0, 0);
                region_edge_ori(region_edge_num, 1) = this_edge(0, 1);
                region_edge_num += 1;
            }
        }
    }
    region_edge = region_edge_ori.block(0, 0, region_edge_num, 3);
    //提取区域中的端点
    MyMatrixXf row_index;
    MyMatrixXf col_index = MyMatrixXf::Ones(1, 1);
    for (int i = 0; i < region_edge_num; i++) {
        for (int j = 0; j < 2; j++) {
            if (!anyIsmember(region_point_ori.block(0, 0, region_point_num, 1), region_edge(i, j))) {
                region_point_ori(region_point_num, 0) = region_edge(i, j);
                if (region_edge(i, 2) == 0) region_point_ori(region_point_num, 1) = 0;
                region_point_num += 1;
            }
            else {
                if (region_edge(i, 2) == 0) {
                    row_index = getIndex(region_point_ori.col(0), region_edge(i, j), 0);
                    valueSubMat(region_point_ori, row_index, col_index, 0);
                }
            }
        }
    }
    region_point = region_point_ori.block(0, 0, region_point_num, 2);
}

float getRegionArea(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx) {
    int face_num = region_face_idx.rows();
    float area = 0;
//#pragma omp parallel for
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    for (int i = 0; i < face_num; i++) {
        tri = getSubMat_Rows(point, face.row(region_face_idx(i, 0)));
        area += getTriArea(tri);
    }
    return area;
}

MyMatrixXf getPointNeighFace(const MyMatrixXf& point, const MyMatrixXf& PFneighbor) {
    //查询与点集相邻接的面(输入的是点的id)
    int point_num = point.rows();
    MyMatrixXf neigh_face_idx = MyMatrixXf::Zero(29 * point_num, 1);
    int neigh_face_num = 0;
    int r;
    for (int i = 0; i < point_num; i++) {
        r = point(i, 0);
        for (int j = 0; j < PFneighbor(r, 0); j++) {
            if (!anyIsmember(neigh_face_idx.block(0, 0, neigh_face_num, 1), PFneighbor(r, j + 1))) {
                neigh_face_idx(neigh_face_num, 0) = PFneighbor(r, j + 1);
                neigh_face_num += 1;
            }
        }
    }
    return neigh_face_idx.block(0, 0, neigh_face_num, 1);
}

void getOutlineSort(const MyMatrixXf& region_face, MyMatrixXf& outline_point_sort, MyMatrixXf& outline_sort) {
    //提取影响域的边界
    //注意：若面按照逆时针存储（1，2——2，3——3，1），这里提取到的边界也是逆时针顺序
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);

    MyMatrixXf outline = MyMatrixXf::Zero(region_edge.rows() - region_edge.col(2).sum(), 2);
    int outline_num = 0;
    for (int i = 0; i < region_edge.rows(); i++) {
        if (region_edge(i, 2) == 0) {
            outline.row(outline_num) = region_edge.block(i, 0, 1, 2);
            outline_num += 1;
        }
    }

    //提取出按顺序排列的边界
    outline_sort = regionEdgeSort(outline);
    outline_point_sort = outline_sort.col(0);
}

MyMatrixXf regionEdgeSort(MyMatrixXf outline) {
    //将边界排序
    int outline_num = outline.rows();
    MyMatrixXf t = MyMatrixXf::Zero(1, 2);
    for (int i = 0; i < outline_num - 1; i++) {
        if (outline(i, 1) == outline(i + 1, 0))continue;
        for (int j = i + 2; j < outline_num; j++) {
            if (outline(i, 1) == outline(j, 0)) {
                t = outline.row(j);
                outline.row(j) = outline.row(i + 1);
                outline.row(i + 1) = t;
            }
        }
    }
    return outline;
}

MyMatrixXf getFacesNormal(const MyMatrixXf& point, const MyMatrixXf& face) {
    int face_num = face.rows();
    MyMatrixXf normal = MyMatrixXf::Zero(face_num, 3);
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    for (int i = 0; i < face_num; i++) {
        tri = getSubMat_Rows(point, face.row(i));
        normal.row(i) = getPerFaceNormal(tri);
    }
    return normal;
}

MyMatrixXf getPerFaceNormal(const MyMatrixXf& tri) {
    MyMatrixXf normal = crossProduct(tri.row(1) - tri.row(0), tri.row(2) - tri.row(1));
    /*if (globalTest == 12) {
        cout << tri << endl << endl;
        cout << normal << endl << endl;
    }*/
    if (norm(normal) != 0) {
        normal = normal / norm(normal);
    }
    return normal;
}

MyMatrixXf findPointNeighFace(const MyMatrixXf& face, float this_point) {
    //寻找与当前点相邻接的三角面（输入点的id）
    int face_num = face.rows();
    MyMatrixXf neigh_face_id = MyMatrixXf::Zero(face_num, 1);
    int count = 0;
    for (int i = 0; i < face_num; i++) {
        if (anyIsmember(face.row(i), this_point)) {
            neigh_face_id(count, 0) = i;
            count += 1;
        }
    }
    return neigh_face_id.block(0, 0, count, 1);
}

MyMatrixXf findEdgeNeighFace(const MyMatrixXf& face, const MyMatrixXf& this_edge) {
    //寻找与当前线段相邻接的三角面
    int face_num = face.rows();
    MyMatrixXf neigh_face_id = MyMatrixXf::Zero(2, 1);
    int count = 0;
    for (int i = 0; i < face_num; i++) {
        if (anyIsmember(face.row(i), this_edge(0,0))&& anyIsmember(face.row(i), this_edge(0, 1))) {
            neigh_face_id(count, 0) = i;
            count += 1;
        }
        if (count == 2)break;
    }
    return neigh_face_id.block(0, 0, count, 1);
}

float findCertainNeighFace(const MyMatrixXf& face, const MyMatrixXf& face_idx, const MyMatrixXf& line, float this_face_idx) {
    //寻找与当前面某一边相邻的另一个面
    float third_face_id = -1;
    int face_num = face.rows();
    for (int i = 0; i < face_num; i++) {
        if (anyIsmember(face.row(i), line(0, 0)) && anyIsmember(face.row(i), line(0, 1))) {
            if (face_idx(i, 0) != this_face_idx) {
                third_face_id = i;
                break;
            }
        }
    }
    return third_face_id;
}

MyMatrixXf facePointUpdate(const MyMatrixXf& old_face, const MyMatrixXf& map) {
    //基于map修改面表对应的点号
    int face_num = old_face.rows();
    int r;
    MyMatrixXf new_face = MyMatrixXf::Zero(face_num, 3);
    for (int i = 0; i < face_num; i++) {
        for (int j = 0; j < 3; j++) {
            r = old_face(i, j);
            new_face(i, j) = map(r, 0);
        }
    }
    return new_face;
}

MyMatrixXf findConnectRegion(const MyMatrixXf& face, const MyMatrixXf& region_face_idx) {
    //将线段支撑域按照是否连通划分为不同的区域
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    int region_face_num = region_face_idx.rows();
    MyMatrixXf connect_region = MyMatrixXf::Zero(region_face_num, region_face_num + 1);
    int connect_region_count = 0;
    MyMatrixXf region_face_sign = MyMatrixXf::Zero(region_face_num, 1);
    int r;
    MyMatrixXf seed, this_face_idx;
    MyMatrixXf new_face_idx = MyMatrixXf::Zero(1, face.rows());
    MyMatrixXf edge = MyMatrixXf::Zero(3, 2);
    while (region_face_sign.sum() != region_face_num) {
        //选出新的连通区域的种子面，并进行相应的初始化操作
        seed = getIndex(region_face_sign, 0, 0);
        this_face_idx = getSubMat_Rows(region_face_idx, seed.block(0, 0, 1, 1));
        connect_region(connect_region_count, 0) = 1;
        connect_region(connect_region_count, 1) = this_face_idx(0, 0);
        r = seed(0, 0);
        region_face_sign(r, 0) = 1;
        while (this_face_idx.cols() != 0) {        
            int new_face_count = 0;
            for (int i = 0; i < this_face_idx.cols(); i++) {
                r = this_face_idx(0, i);
                edge << face(r, 0), face(r, 1),
                    face(r, 1), face(r, 2),
                    face(r, 2), face(r, 0);
                for (int j = 0; j < 3; j++) {
                    int third_face_id = findCertainNeighFace(region_face, region_face_idx, edge.row(j), r);
                    if (third_face_id == -1)continue;
                    if (region_face_sign(third_face_id, 0) == 0) {
                        new_face_idx(0, new_face_count) = region_face_idx(third_face_id, 0);
                        new_face_count += 1;
                        region_face_sign(third_face_id, 0) = 1;
                    }
                }
            }
            this_face_idx = new_face_idx.block(0, 0, 1, new_face_count);
            connect_region.block(connect_region_count, connect_region(connect_region_count, 0) + 1, 1, new_face_count) = this_face_idx;
            connect_region(connect_region_count, 0) += new_face_count;
        }
        connect_region_count += 1;
    }
    return connect_region.block(0, 0, connect_region_count, region_face_num + 1);
}

MyMatrixXf findConnectEdge(const MyMatrixXf& edge) {
    //寻找边的连通成分
    //若连通成分的首尾点相同，则表示连通成分闭合
    //考虑一点多连的情况（目前最多考虑一点三连）
    int edge_num = edge.rows();
    int max_length = edge_num + 1;
    MyMatrixXf used = MyMatrixXf::Zero(edge_num, 1);
    int connect_edge_num = 0;

    MyMatrixXf connect_edge = MyMatrixXf::Zero(edge_num, max_length + 1);
    MyMatrixXf this_connect = -1 * MyMatrixXf::Ones(1, max_length);
    for (int i = 0; i < edge_num; i++) {
        if (used(i, 0) == 1)continue;
        //初始化
        used(i, 0) = 1;
        this_connect = -1 * MyMatrixXf::Ones(1, max_length);
        int left_neigh = -1;
        this_connect(0, max_length - 1) = edge(i, 0);
        this_connect(0, 0) = edge(i, 1);
        int left_neigh_num = 0;
        
        while (left_neigh != 0) {
            left_neigh = 0;
            for (int j = 0; j < edge_num; j++) {
                if (used(j, 0) == 1)continue;
                for (int k = 0; k < 2; k++) {
                    if (edge(j, k) == this_connect(0, max_length - 1 - left_neigh_num)) {
                        left_neigh_num += 1;
                        this_connect(0, max_length - 1 - left_neigh_num) = edge(j, 1 - k);
                        used(j, 0) = 1;
                        left_neigh = 1;
                        break;
                    }
                }
            }
        }

        int right_neigh = -1;
        int right_neigh_num = 0;
        while (right_neigh != 0) {
            right_neigh = 0;
            for (int j = 0; j < edge_num; j++) {
                if (used(j, 0) == 1)continue;
                for (int k = 0; k < 2; k++) {
                    if (edge(j, k) == this_connect(0, right_neigh_num)) {
                        right_neigh_num += 1;
                        this_connect(0, right_neigh_num) = edge(j, 1 - k);
                        used(j, 0) = 1;
                        right_neigh = 1;
                        break;
                    }
                }
            }
        }
        connect_edge(connect_edge_num, 0) = 2 + left_neigh_num + right_neigh_num;
        connect_edge.block(connect_edge_num, 1, 1, left_neigh_num + 1) = this_connect.block(0, max_length - 1 - left_neigh_num, 1, left_neigh_num + 1);
        connect_edge.block(connect_edge_num, left_neigh_num + 2, 1, right_neigh_num + 1) = this_connect.block(0, 0, 1, right_neigh_num + 1);
        connect_edge_num += 1;
    }
    return connect_edge.block(0, 0, connect_edge_num, max_length + 1);
}

MyMatrixXf getFaceEdge(const MyMatrixXf& face) {
    MyMatrixXf face_edge = MyMatrixXf::Zero(3, 2);
    face_edge << face(0, 0), face(0, 1),
        face(0, 1), face(0, 2),
        face(0, 2), face(0, 0);
    return face_edge;
}