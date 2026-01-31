// main
#include <iostream>
#include "class.h"
#include "file.h"
#include "lineProcess.h"
#include "myMath.h"
#include "dist.h"
#include "meshOper.h"
#include "region.h"
#include "retriangulate.h"
#include "time.h"

bool readTXT(MyMatrixXf& mat, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        cout << "File open error!" << endl;
        return false;
    }
    MyMatrixXf mat_ori = MyMatrixXf::Zero(999999, 1);
    int edge_num = 0;
    while (!feof(fp)) {
        fscanf(fp, "%f", &mat_ori(edge_num, 0));
        if (mat_ori(edge_num, 0) != 0)
            edge_num++;
    }
    fclose(fp);

    mat = mat_ori.block(0, 0, edge_num, 1);
    return true;
}

bool writeTXT(const MyMatrixXf& mat, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        cout << "File open error!" << endl;
        return false;
    }
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++)
            fprintf(fp, "%f ", mat(i, j));
        fprintf(fp, "\n");
    }
    fclose(fp);
    return true;
}

int globalTest = 0;
int dp_count = 0;
clock_t myTime = 0;
clock_t myTime_2 = 0;

int main(int argc, char** argv)
{
    if (argc < 7) {
        cout << "Usage: " << argv[0]
            << " <line_obj> <mesh_obj> <out_dir> <para_1> <para_2> <node_interval>\n";
        cout << "Example: " << argv[0]
            << " mesh_optimization.exe line.obj mesh_5000.obj C:\\mesh_optimization\\C++\\out 30 0.7 0.03\n";
        return -1;
    }

    string line_file = argv[1];
    string mesh_file = argv[2];
    string output_dir = argv[3];

    float para_1 = atof(argv[4]);
    float para_2 = atof(argv[5]);
    float node_interval = atof(argv[6]);
    float sample_interval_1 = node_interval / 20;
    float sample_interval_2 = node_interval * 5;

    cout << "Input parameters:\n"
        << " line_obj: " << line_file << "\n"
        << " mesh_obj: " << mesh_file << "\n"
        << " output_dir: " << output_dir << "\n"
        << " para_1=" << para_1 << " para_2=" << para_2 << " node_interval=" << node_interval << endl;

    clock_t start, end;
    MyMatrixXf point, face, line_set;
    int k_1 = 3;
    int k_2 = 500;
    readLineOBJ(line_set, (output_dir + "/" + line_file).c_str());
    readMeshOBJ(point, face, (output_dir + "/" + mesh_file).c_str());
    int line_num = line_set.rows() / 2;
    int pts_num = point.rows();
    int face_num = face.rows();
    MyMatrixXf constrained_line_set;

    MyMatrixXf para_1_set = para_1 * MyMatrixXf::Ones(line_num, 1);
    MyMatrixXf para_2_set = para_2 * MyMatrixXf::Ones(line_num, 1);

    // Start timing
    start = clock();
    lineSort(line_set);
    MyMatrixXf each_node_num = MyMatrixXf::Zero(line_num, 1);
    MyMatrixXf node_ori = MyMatrixXf::Zero(line_num * 1000, 3);
    int count = 0;
    for (int i = 0; i < line_num; i++) {
        MyMatrixXf temp = getNode(line_set.block(2 * i, 0, 2, 3), node_interval);
        each_node_num(i, 0) = temp.rows();
        node_ori.block(count, 0, each_node_num(i, 0), 3) = temp;
        count += each_node_num(i, 0);
    }
    MyMatrixXf node = node_ori.block(0, 0, count, 3);
    int all_node_num = node.rows();
    MyMatrixXf new_point = MyMatrixXf::Zero(pts_num + all_node_num, 3);
    new_point << point,
        node;
    MyMatrixXf sign = MyMatrixXf::Zero(pts_num + all_node_num, 1);
    for (int i = 0; i < pts_num; i++) {
        if (!(point(i, 0) == 0 && point(i, 1) == 0 && point(i, 2) == 0))
            sign(i, 0) = 1;
    }
    MyMatrixXf idx = MyMatrixXf::Zero(all_node_num, 1);
    for (int i = 0; i < all_node_num; i++) {
        idx(i, 0) = i + pts_num;
    }

    // Build k-d tree
    pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(new_point);

    MyMatrixXf successful = MyMatrixXf::Zero(line_num, 1);
    // Start embedding lines
    for (int m = 0; m < line_num; m++) {
        try {
            // Build point-face adjacency table
            cout << m << endl;
            globalTest = m;
            dp_count = 0;

            MyMatrixXf PFneighbor;
            createPFneighbor(new_point, face, PFneighbor);

            // LSM initialization
            MyMatrixXf knn_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_1, sample_interval_1);
            knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
            knn_face_idx = regionConnect(new_point, face, knn_face_idx, PFneighbor);
            // Line-to-mesh projection
            MyMatrixXf base_face_idx = isIntersectant(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1);
            base_face_idx = holeFilling(base_face_idx, face, PFneighbor);
            base_face_idx = regionConnect(new_point, face, base_face_idx, PFneighbor);
            int base_face_num = base_face_idx.rows();
            cout << "region projected!" << endl;
            // LSM growing
            MyMatrixXf range_face_idx_ori = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_2, sample_interval_2);
            MyMatrixXf range_face_idx = MyUnion(range_face_idx_ori, base_face_idx);
            int range_face_num = range_face_idx.rows();
            MyMatrixXf l2f_dist = MyMatrixXf::Zero(range_face_num, 1);
            for (int i = 0; i < range_face_num; i++) {
                l2f_dist(i, 0) = line2faceDist2(getSubMat_Rows(new_point, face.row(range_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
            }
            MyMatrixXf normal = getFacesNormal(new_point, getSubMat_Rows(face, range_face_idx));
            float* l2f_dist_base = new float[base_face_num];
            for (int i = 0; i < base_face_num; i++) {
                l2f_dist_base[i] = line2faceDist2(getSubMat_Rows(new_point, face.row(base_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
            }
            float extend_threshold = getQuantile(l2f_dist_base, base_face_num, para_2_set(m, 0));
            delete[] l2f_dist_base;

            MyMatrixXf extend_face_idx = regionExtend(face, range_face_idx, base_face_idx, normal, l2f_dist, PFneighbor, para_1_set(m, 0), extend_threshold);
            extend_face_idx = shapeRepairing(new_point, face, extend_face_idx, PFneighbor);
            extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor);
            cout << "region extended!" << endl;
            
            // LSM clipping
            MyMatrixXf segment_face_idx = regionSegment(new_point, face, extend_face_idx, base_face_idx, line_set.block(2 * m, 0, 2, 3), constrained_line_set, sample_interval_1);
            segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor);
            cout << "region segmented!" << endl;
            if (haveConstrainedEdge(face, segment_face_idx, constrained_line_set))
                continue;

            MyMatrixXf node_idx;
            if (m == 0)
                node_idx = idx.block(0, 0, each_node_num(m, 0), 1);
            else
                node_idx = idx.block(each_node_num.block(0, 0, m, 1).sum(), 0, each_node_num(m, 0), 1);
            
            //LSM retriangulation
            bool error = regionRemeshing(new_point, face, sign, segment_face_idx, node_idx, PFneighbor, constrained_line_set, 0.5, 0.5, m);
            if (error == false) {
                cout << "remeshing successfully!" << endl;
                successful(m, 0) = 1;
            }
        }
        catch (int n) {
            n = m;
            cout << "Unknown error occurred during the embedding of line " << n << "!" << endl;
            continue;
        }
    }
    end = clock();
    cout << double(end - start) << endl; 
    cout << double(myTime) << endl;

    string out_file = output_dir + "\\" + mesh_file + "_out.obj";
    saveMeshOBJ(new_point, face, sign, out_file.c_str());
}