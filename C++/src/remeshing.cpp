//
//void coMatClustering(MyMatrixXf& coMat,MyMatrixXf& cluster) {
//    Graph G;
//    for (int i = 0; i < coMat.rows(); i++) {
//        for (int j = i + 1; j < coMat.cols(); j++) {
//            if (coMat(i, j) == 1) {
//                add_edge(i, j, 1, G);
//            }
//        }
//    }
//    MyMatrixXf sign = MyMatrixXf::Zero(coMat.rows(), 1);
//    MyMatrixXf cluster_ori = MyMatrixXf::Zero(coMat.rows(), coMat.rows() + 1);
//    int cluster_num = 0;
//    for (int i = 0; i < coMat.rows(); i++) {
//        if (sign(i, 0) == 1)
//            continue;
//        if (coMat.row(i).sum() == 0)
//            continue;
//        // 初始化可达性标记数组
//        std::vector<bool> reachable(num_vertices(G), false);
//
//        // 执行广度优先搜索，从 start 顶点开始
//        ReachableVisitor vis(reachable);
//        breadth_first_search(G,
//            vertex(i, G),
//            visitor(vis));
//
//        // 输出所有可以到达 start 点的顶点
//        int count = 0;
//        for (Vertex v = 0; v < num_vertices(G); ++v) {
//            if (reachable[v]) {
//                cluster_ori(cluster_num, (count++) + 1) = v;
//                cluster_ori(cluster_num, 0) += 1;
//                sign(v, 0) = 1;
//            }
//        }
//        cluster_num += 1;
//    }
//    cluster = cluster_ori.block(0, 0, cluster_num, cluster_ori.cols());
//}
//struct region {
//    MyMatrixXf region_face_idx;
//};
//
//int main()
//{
//    clock_t start, end;
//    
//    float para_1 = 45;
//    float para_2 = 1.65;
//    float node_interval = 5;
//    float sample_interval_1 = 0.25;
//    float sample_interval_2 = 25;
//    int k_1 = 3;
//    int k_2 = 200;
//    MyMatrixXf point, face, line_set;
//    readLineOBJ(line_set, "F:\\3D_reconstruction_test\\Data\\guangzhou\\line.obj");
//    readMeshOBJ(point, face, "F:\\3D_reconstruction_test\\Data\\guangzhou\\mesh_50000.obj");
//    //readLineOBJ(line_set, "F:\\test\\line_2.obj");
//    //readMeshOBJ(point, face, "F:\\test\\timber_out_1.obj");
//    int line_num = line_set.rows() / 2;
//    int pts_num = point.rows();
//    int face_num = face.rows();
//    MyMatrixXf constrained_line_set;
//    //readConstrainedEdgeTXT(constrained_line_set, "F:\\test\\constrained.txt");
//    //readConstrainedEdgeTXT(constrained_line_set, "F:\\3D_reconstruction_test\\Data\\timber_frame_house\\constrained_1.txt");
//
//    //开始计时
//    start = clock();
//    lineSort(line_set);
//    MyMatrixXf each_node_num = MyMatrixXf::Zero(line_num, 1);
//    MyMatrixXf node_ori = MyMatrixXf::Zero(line_num * 1000, 3);
//    int count = 0;
//    for (int i = 0; i < line_num; i++) {
//        MyMatrixXf temp = getNode(line_set.block(2 * i, 0, 2, 3), node_interval);
//        each_node_num(i, 0) = temp.rows();
//        node_ori.block(count, 0, each_node_num(i, 0), 3) = temp;
//        count += each_node_num(i, 0);
//    }
//    MyMatrixXf node = node_ori.block(0, 0, count, 3);
//    int all_node_num = node.rows();
//    MyMatrixXf new_point = MyMatrixXf::Zero(pts_num + all_node_num, 3);
//    new_point << point,
//        node;
//    MyMatrixXf sign = MyMatrixXf::Zero(pts_num + all_node_num, 1);
//    for (int i = 0; i < pts_num; i++) {
//        if (!(point(i, 0) == 0 && point(i, 1) == 0 && point(i, 2) == 0))
//            sign(i, 0) = 1;
//    }
//    MyMatrixXf idx = MyMatrixXf::Zero(all_node_num, 1);
//    for (int i = 0; i < all_node_num; i++) {
//        idx(i, 0) = i + pts_num;
//    }
//
//    //构建k-d树
//    pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(new_point);
//
//    MyMatrixXf PFneighbor;
//    createPFneighbor(new_point, face, PFneighbor);
//    //构建关联矩阵
//    MyMatrixXf coMat = MyMatrixXf::Zero(line_num, line_num);
//    vector<region> range(line_num);
//    for (int m = 0; m < line_num; m++) {
//        range[m].region_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_2, sample_interval_2);
//    }
//    for (int i = 0; i < line_num; i++) {
//        for (int j = i; j < line_num; j++) {
//            if (i == j) {
//                coMat(i, j) = 1;
//                continue;
//            }
//            MyMatrixXf temp = MyIntersection(range[i].region_face_idx, range[j].region_face_idx);
//            if (temp.rows() != 0) {
//                coMat(i, j) = 1;
//                coMat(j, i) = 1;
//            }
//        }
//    }
//    MyMatrixXf cluster;
//    coMatClustering(coMat, cluster);
//    cout << cluster << endl;
//    //开始嵌入线
//#pragma omp parallel for
//    for (int j = 0; j < cluster.rows(); j++) {
//        for (int k = 0; k < cluster(j, 0); k++) {
//            //构建点面邻接表
//            int m = cluster(j, k + 1);
//            cout << m << endl;
//
//            createPFneighbor(new_point, face, PFneighbor);
//
//            //计算knn支撑域
//            MyMatrixXf knn_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_1, sample_interval_1);
//            knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
//            knn_face_idx = regionConnect(new_point, face, knn_face_idx, PFneighbor);
//
//            //end = clock();
//            //cout << double(end - start) << endl;
//
//            //计算基准支撑域
//            //start = clock();
//            MyMatrixXf base_face_idx = isIntersectant(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1);
//            base_face_idx = holeFilling(base_face_idx, face, PFneighbor);
//            base_face_idx = regionConnect(new_point, face, base_face_idx, PFneighbor);
//            int base_face_num = base_face_idx.rows();
//
//            //计算扩展支撑域
//            //start = clock();
//            MyMatrixXf range_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_2, sample_interval_2);
//            int range_face_num = range_face_idx.rows();
//            MyMatrixXf l2f_dist = MyMatrixXf::Zero(range_face_num, 1);
//
//            for (int i = 0; i < range_face_num; i++) {
//                l2f_dist(i, 0) = line2faceDist2(getSubMat_Rows(new_point, face.row(range_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//            }
//            MyMatrixXf normal = getFacesNormal(new_point, getSubMat_Rows(face, range_face_idx));
//            float* l2f_dist_base = new float[base_face_num];
//            for (int i = 0; i < base_face_num; i++) {
//                l2f_dist_base[i] = line2faceDist2(getSubMat_Rows(new_point, face.row(base_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//            }
//            float max_dist = getMedian(l2f_dist_base, base_face_num);
//            delete[] l2f_dist_base;
//            MyMatrixXf extend_face_idx = regionExtend(face, range_face_idx, base_face_idx, normal, l2f_dist, PFneighbor, para_1, para_2 * max_dist);
//            extend_face_idx = shapeRepairing(new_point, face, extend_face_idx, PFneighbor);
//            extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor);
//
//            MyMatrixXf segment_face_idx = regionSegment(new_point, face, extend_face_idx, base_face_idx, line_set.block(2 * m, 0, 2, 3), constrained_line_set);
//            segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor);
//            if (haveConstrainedEdge(face, segment_face_idx, constrained_line_set))
//                continue;
//
//            //writeTXT(range_face_idx, "F:\\test\\segment.txt");
//            MyMatrixXf node_idx;
//            if (m == 0)
//                node_idx = idx.block(0, 0, each_node_num(m, 0), 1);
//            else
//                node_idx = idx.block(each_node_num.block(0, 0, m, 1).sum(), 0, each_node_num(m, 0), 1);
//            //start = clock();
//            bool error = regionRemeshing(new_point, face, sign, segment_face_idx, node_idx, PFneighbor, constrained_line_set, 0.5, 0.5);
//            //saveMeshOBJ(new_point, face, sign, "F:\\test\\dublin_2_out.obj");
//        }
//    } 
//    end = clock();
//    cout << double(end - start) << endl;
//    saveMeshOBJ(new_point, face, sign, "F:\\test\\guangzhou_out.obj");
//}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件



//// MeshOptimization.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
////
//
//#include <iostream>
//#include "class.h"
//#include "file.h"
//#include "lineProcess.h"
//#include "myMath.h"
//#include "dist.h"
//#include "meshOper.h"
//#include "region.h"
//#include "retriangulate.h"
//#include "time.h"
//
//std::string concatenateStrings(const char* str1, int number, const char* str2) {
//    // 直接用 std::string 操作
//    std::string result = std::string(str1) + std::to_string(number) + str2;
//    return result; // std::string 自动管理内存
//}
//
//bool readTXT(MyMatrixXf& mat, const char* filename) {
//    FILE* fp = fopen(filename, "r");
//    if (fp == NULL) {
//        cout << "文件打开异常！" << endl;
//        return false;
//    }
//    MyMatrixXf mat_ori = MyMatrixXf::Zero(999999, 1);
//    int edge_num = 0;
//    while (!feof(fp)) {
//        fscanf(fp, "%f", &mat_ori(edge_num, 0));
//        if (mat_ori(edge_num, 0) != 0)
//            edge_num++;
//    }
//    fclose(fp);
//
//    mat = mat_ori.block(0, 0, edge_num, 1);
//    return true;
//}
//
//bool writeTXT(const MyMatrixXf& mat, const char* filename) {
//    FILE* fp = fopen(filename, "w");
//    if (fp == NULL) {
//        cout << "文件打开异常！" << endl;
//        return false;
//    }
//    for (int i = 0; i < mat.rows(); i++) {
//        for (int j = 0; j < mat.cols(); j++)
//            fprintf(fp, "%f ", mat(i, j));
//        fprintf(fp, "\n");
//    }
//    fclose(fp);
//    return true;
//}
//
//int globalTest = 0;
//int dp_count = 0;
//clock_t myTime = 0;
//clock_t myTime_2 = 0;
//
//
//int main()
//{
//    clock_t start, end;
//    MyMatrixXf point, face, line_set;
//    float para_1 = 30;
//    float para_2 = 1.25;
//    float node_interval = 25;
//    float sample_interval_1 = 1.25;
//    float sample_interval_2 = 25;
//    int k_1 = 3;
//    int k_2 = 250;
//    readLineOBJ(line_set, "G:\\3D_reconstruction_test_weight\\simulate3\\1\\line_0.obj");
//    int line_num = line_set.rows() / 2;
//    MyMatrixXf para_1_set = para_1 * MyMatrixXf::Ones(line_num, 1);
//    MyMatrixXf para_2_set = para_2 * MyMatrixXf::Ones(line_num, 1);
//
//    for (int m = 0; m < line_num; m++) {
//        readMeshOBJ(point, face, "G:\\3D_reconstruction_test_weight\\simulate3\\1\\input_18.obj");
//
//        int pts_num = point.rows();
//        int face_num = face.rows();
//        MyMatrixXf constrained_line_set;
//
//        //开始计时
//        //start = clock();
//        lineSort(line_set);
//        MyMatrixXf each_node_num = MyMatrixXf::Zero(line_num, 1);
//        MyMatrixXf node_ori = MyMatrixXf::Zero(line_num * 1000, 3);
//        int count = 0;
//        for (int i = 0; i < line_num; i++) {
//            MyMatrixXf temp = getNode(line_set.block(2 * i, 0, 2, 3), node_interval);
//            each_node_num(i, 0) = temp.rows();
//            node_ori.block(count, 0, each_node_num(i, 0), 3) = temp;
//            count += each_node_num(i, 0);
//        }
//        MyMatrixXf node = node_ori.block(0, 0, count, 3);
//        int all_node_num = node.rows();
//        MyMatrixXf new_point = MyMatrixXf::Zero(pts_num + all_node_num, 3);
//        new_point << point,
//            node;
//        MyMatrixXf sign = MyMatrixXf::Zero(pts_num + all_node_num, 1);
//        for (int i = 0; i < pts_num; i++) {
//            if (!(point(i, 0) == 0 && point(i, 1) == 0 && point(i, 2) == 0))
//                sign(i, 0) = 1;
//        }
//        MyMatrixXf idx = MyMatrixXf::Zero(all_node_num, 1);
//        for (int i = 0; i < all_node_num; i++) {
//            idx(i, 0) = i + pts_num;
//        }
//
//        //构建k-d树
//        pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(new_point);
//
//        //开始嵌入线
//            
//        //构建点面邻接表
//        cout << m << endl;
//        //if (m != 0)
//        //    continue;
//        globalTest = m;
//        dp_count = 0;
//
//        MyMatrixXf PFneighbor;
//        createPFneighbor(new_point, face, PFneighbor);
//
//        //计算knn支撑域
//        MyMatrixXf knn_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_1, sample_interval_1);
//        knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
//        knn_face_idx = regionConnect(new_point, face, knn_face_idx, PFneighbor);
//
//        //计算基准支撑域
//        MyMatrixXf base_face_idx = isIntersectant(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1);
//        base_face_idx = holeFilling(base_face_idx, face, PFneighbor);
//        base_face_idx = regionConnect(new_point, face, base_face_idx, PFneighbor);
//        int base_face_num = base_face_idx.rows();
//        cout << "region projected!" << endl;
//
//        //计算扩展支撑域
//        MyMatrixXf range_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_2, sample_interval_2);
//        int range_face_num = range_face_idx.rows();
//        MyMatrixXf l2f_dist = MyMatrixXf::Zero(range_face_num, 1);
//
//        for (int i = 0; i < range_face_num; i++) {
//            l2f_dist(i, 0) = line2faceDist2(getSubMat_Rows(new_point, face.row(range_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//        }
//        MyMatrixXf normal = getFacesNormal(new_point, getSubMat_Rows(face, range_face_idx));
//        float* l2f_dist_base = new float[base_face_num];
//        for (int i = 0; i < base_face_num; i++) {
//            l2f_dist_base[i] = line2faceDist2(getSubMat_Rows(new_point, face.row(base_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//        }
//        float max_dist = getMedian(l2f_dist_base, base_face_num);
//        delete[] l2f_dist_base;
//
//        MyMatrixXf extend_face_idx = regionExtend(face, range_face_idx, base_face_idx, normal, l2f_dist, PFneighbor, para_1_set(m, 0), para_2_set(m, 0) * max_dist);
//        cout << "region extended!" << endl;
//        extend_face_idx = shapeRepairing(new_point, face, extend_face_idx, PFneighbor);
//        extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor, range_face_idx);
//
//
//        MyMatrixXf segment_face_idx = regionSegment(new_point, face, extend_face_idx, base_face_idx, line_set.block(2 * m, 0, 2, 3), constrained_line_set, sample_interval_1);
//        segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor);
//        cout << "region segmented!" << endl;
//        if (haveConstrainedEdge(face, segment_face_idx, constrained_line_set))
//            continue;
//
//        MyMatrixXf node_idx;
//        if (m == 0)
//            node_idx = idx.block(0, 0, each_node_num(m, 0), 1);
//        else
//            node_idx = idx.block(each_node_num.block(0, 0, m, 1).sum(), 0, each_node_num(m, 0), 1);
//
//        bool error = regionRemeshing(new_point, face, sign, segment_face_idx, node_idx, PFneighbor, constrained_line_set, 0.5, 0.5);
//        if (error == false)
//            cout << "remeshing successfully!" << endl;
//
//
//        saveMeshOBJ(new_point, face, sign, concatenateStrings("G:\\3D_reconstruction_test_weight\\simulate3\\1\\PC\\18\\", m, "_p.obj").c_str());
//        cout << double(myTime) << endl;
//    }  
//}










//// MeshOptimization.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
////
//
//#include <iostream>
//#include "class.h"
//#include "file.h"
//#include "lineProcess.h"
//#include "myMath.h"
//#include "dist.h"
//#include "meshOper.h"
//#include "region.h"
//#include "retriangulate.h"
//#include "time.h"
//
//bool readTXT(MyMatrixXf& mat, const char* filename) {
//    FILE* fp = fopen(filename, "r");
//    if (fp == NULL) {
//        cout << "文件打开异常！" << endl;
//        return false;
//    }
//    MyMatrixXf mat_ori = MyMatrixXf::Zero(999999, 1);
//    int edge_num = 0;
//    while (!feof(fp)) {
//        fscanf(fp, "%f", &mat_ori(edge_num, 0));
//        if (mat_ori(edge_num, 0) != 0)
//            edge_num++;
//    }
//    fclose(fp);
//
//    mat = mat_ori.block(0, 0, edge_num, 1);
//    return true;
//}
//
//bool writeTXT(const MyMatrixXf& mat, const char* filename) {
//    FILE* fp = fopen(filename, "w");
//    if (fp == NULL) {
//        cout << "文件打开异常！" << endl;
//        return false;
//    }
//    for (int i = 0; i < mat.rows(); i++) {
//        for (int j = 0; j < mat.cols(); j++)
//            fprintf(fp, "%f ", mat(i, j));
//        fprintf(fp, "\n");
//    }
//    fclose(fp);
//    return true;
//}
//
//int globalTest = 0;
//int dp_count = 0;
//clock_t myTime = 0;
//clock_t myTime_2 = 0;
//
//
//int main()
//{
//    clock_t start, end;
//    MyMatrixXf point, face, line_set;
//    float para_1 = 45;
//    float para_2 = 1.25;
//    //float dist_threshold = 0.006;
//    float node_interval = 0.03;
//    float sample_interval_1 = node_interval / 20;
//    float sample_interval_2 = node_interval * 5;
//    int k_1 = 3;
//    int k_2 = 250;
//    readLineOBJ(line_set, "D:\\C_code\\remeshing\\guangzhou_boundary_2.obj");
//    readMeshOBJ(point, face, "D:\\C_code\\remeshing\\guangzhou_mesh2.obj");
//    //readLineOBJ(line_set, "G:\\3D_reconstruction_test_weight\\Data\\dublin_2\\line.obj");
//    //readMeshOBJ(point, face, "G:\\3D_reconstruction_test_weight\\Data\\dublin_2\\mesh_10000.obj");
//    //readLineOBJ(line_set, "F:\\test\\line_2.obj");
//    //readMeshOBJ(point, face, "F:\\test\\timber_out_1.obj");
//    int line_num = line_set.rows() / 2;
//    int pts_num = point.rows();
//    int face_num = face.rows();
//    MyMatrixXf constrained_line_set;
//
//    MyMatrixXf para_1_set = para_1 * MyMatrixXf::Ones(line_num, 1);
//    MyMatrixXf para_2_set = para_2 * MyMatrixXf::Ones(line_num, 1);
//
//    //MyMatrixXf para_1_set = 15 * MyMatrixXf::Ones(line_num, 1);
//    //MyMatrixXf para_2_set = MyMatrixXf::Ones(line_num, 1);
//    //for (int i = 0; i < line_num; i++) {
//    //    if (i == 0 || i == 1 ||i==2||i==3|| i == 27 || i == 28||i==29||i==30 || (i >= 367 && i <= 397) || i >= 431) {
//    //        para_1_set(i, 0) = para_1;
//    //        para_2_set(i, 0) = para_2;
//    //    }
//    //    if (i >= 398 && i <= 430) {
//    //        para_1_set(i, 0) = 0;
//    //        para_2_set(i, 0) = 0;
//    //    }
//    //}
//
//
//    //开始计时
//    //start = clock();
//    lineSort(line_set);
//    MyMatrixXf each_node_num = MyMatrixXf::Zero(line_num, 1);
//    MyMatrixXf node_ori = MyMatrixXf::Zero(line_num * 1000, 3);
//    int count = 0;
//    for (int i = 0; i < line_num; i++) {
//        MyMatrixXf temp = getNode(line_set.block(2 * i, 0, 2, 3), node_interval);
//        each_node_num(i, 0) = temp.rows();
//        node_ori.block(count, 0, each_node_num(i, 0), 3) = temp;
//        count += each_node_num(i, 0);
//    }
//    MyMatrixXf node = node_ori.block(0, 0, count, 3);
//    int all_node_num = node.rows();
//    MyMatrixXf new_point = MyMatrixXf::Zero(pts_num + all_node_num, 3);
//    new_point << point,
//        node;
//    MyMatrixXf sign = MyMatrixXf::Zero(pts_num + all_node_num, 1);
//    for (int i = 0; i < pts_num; i++) {
//        if (!(point(i, 0) == 0 && point(i, 1) == 0 && point(i, 2) == 0))
//            sign(i, 0) = 1;
//    }
//    MyMatrixXf idx = MyMatrixXf::Zero(all_node_num, 1);
//    for (int i = 0; i < all_node_num; i++) {
//        idx(i, 0) = i + pts_num;
//    }
//
//    //构建k-d树
//    pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(new_point);
//
//    MyMatrixXf successful = MyMatrixXf::Zero(line_num, 1);
//    //开始嵌入线
//    for (int m = 0; m < line_num; m++) {
//        try {
//            //构建点面邻接表
//            cout << m << endl;
//            //if (m == 289)
//            //    continue;
//            globalTest = m;
//            dp_count = 0;
//
//            MyMatrixXf PFneighbor;
//            createPFneighbor(new_point, face, PFneighbor);
//
//            //计算knn支撑域
//            MyMatrixXf knn_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_1, sample_interval_1);
//            //writeTXT(knn_face_idx, "F:\\test\\knn.txt");
//            //string filename_8 = "F:\\test\\region_1\\" + to_string(m) + ".obj";
//            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_8.c_str());
//            knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
//            //string filename_6 = "F:\\test\\region_3\\" + to_string(m) + ".obj";
//            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_6.c_str());
//            knn_face_idx = regionConnect(new_point, face, knn_face_idx, PFneighbor);
//            //string filename_7 = "F:\\test\\region_2\\" + to_string(m) + ".obj";
//            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_7.c_str());
//            //计算基准支撑域
//            MyMatrixXf base_face_idx = isIntersectant(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1);
//            //MyMatrixXf base_face_idx = isIntersectant2(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1, dist_threshold);
//            base_face_idx = holeFilling(base_face_idx, face, PFneighbor); 
//            base_face_idx = regionConnect(new_point, face, base_face_idx, PFneighbor);
//            writeTXT(sign, "F:\\test\\sign.txt");
//            int base_face_num = base_face_idx.rows();
//            cout << "region projected!" << endl;
//            string filename_5 = "F:\\test\\region_3\\" + to_string(m) + ".obj";
//            saveSubMeshOBJ(new_point, face, base_face_idx, sign, filename_5.c_str());
//
//            //计算扩展支撑域
//            MyMatrixXf range_face_idx_ori = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_2, sample_interval_2);
//            MyMatrixXf range_face_idx = MyUnion(range_face_idx_ori, base_face_idx);
//            int range_face_num = range_face_idx.rows();
//            MyMatrixXf l2f_dist = MyMatrixXf::Zero(range_face_num, 1);
//
//            for (int i = 0; i < range_face_num; i++) {
//                l2f_dist(i, 0) = line2faceDist2(getSubMat_Rows(new_point, face.row(range_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//            }
//            MyMatrixXf normal = getFacesNormal(new_point, getSubMat_Rows(face, range_face_idx));
//            float* l2f_dist_base = new float[base_face_num];
//            for (int i = 0; i < base_face_num; i++) {
//                l2f_dist_base[i] = line2faceDist2(getSubMat_Rows(new_point, face.row(base_face_idx(i, 0))), line_set.block(2 * m, 0, 2, 3));
//            }
//            float extend_threshold = getMedian(l2f_dist_base, base_face_num);
//            delete[] l2f_dist_base;
//
//            //start = clock();
//            MyMatrixXf extend_face_idx = regionExtend(face, range_face_idx, base_face_idx, normal, l2f_dist, PFneighbor, para_1_set(m, 0), para_2_set(m, 0) * extend_threshold);
//            //end = clock();
//            //cout << double(end - start) << endl;
//            extend_face_idx = shapeRepairing(new_point, face, extend_face_idx, PFneighbor);
//            extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor);
//            //extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor, range_face_idx);
//            string filename_4 = "F:\\test\\region_2\\" + to_string(m) + ".obj";
//            saveSubMeshOBJ2(new_point, face, extend_face_idx, sign, filename_4.c_str());
//            cout << "region extended!" << endl;
//            //writeTXT(extend_face_idx, "F:\\test\\extend_test.txt");
//
//            //start = clock();
//            MyMatrixXf segment_face_idx = regionSegment(new_point, face, extend_face_idx, base_face_idx, line_set.block(2 * m, 0, 2, 3), constrained_line_set, sample_interval_1);
//            //segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor, range_face_idx);
//            segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor);
//            string filename_3 = "F:\\test\\region_1\\" + to_string(m) + ".obj";
//            saveSubMeshOBJ2(new_point, face, segment_face_idx, sign, filename_3.c_str());
//            //end = clock();
//            //cout << double(end - start) << endl;
//            cout << "region segmented!" << endl;
//            if (haveConstrainedEdge(face, segment_face_idx, constrained_line_set))
//                continue;
//            //if (holeDetection(segment_face_idx, face, PFneighbor))
//            //    continue;
//
//            //string filename = "F:\\test\\region\\" + to_string(m) + ".obj";
//            ////const char* file = nullptr;
//            ////file = filename.c_str();
//            //saveSubMeshOBJ(new_point, face, segment_face_idx, sign, filename.c_str());
//
//            //MyMatrixXf segment_face_idx;
//            //readTXT(segment_face_idx, "G:\\3D_reconstruction_test_weight\\simulate3\\1\\LSM\\36_1.txt");
//            //segment_face_idx = segment_face_idx - MyMatrixXf::Ones(segment_face_idx.rows(),1);
//            //writeTXT(segment_face_idx, "F:\\test\\segment_test.txt");
//            MyMatrixXf node_idx;
//            if (m == 0)
//                node_idx = idx.block(0, 0, each_node_num(m, 0), 1);
//            else
//                node_idx = idx.block(each_node_num.block(0, 0, m, 1).sum(), 0, each_node_num(m, 0), 1);
//            //start = clock();
//            bool error = regionRemeshing(new_point, face, sign, segment_face_idx, node_idx, PFneighbor, constrained_line_set, 0.5, 0.5, m);
//            if (error == false) {
//                cout << "remeshing successfully!" << endl;
//                successful(m, 0) = 1;
//            }
//                
//            //writeTXT(constrained_line_set, "F:\\test\\constrained_line_set.txt");
//            //end = clock();
//            //cout << double(end - start) << endl << endl;
//
//            string filename_2 = "F:\\test\\all\\" + to_string(m) + ".obj";
//            saveMeshOBJ(new_point, face, sign, filename_2.c_str());
//            cout << double(myTime) << endl;
//        }
//        catch (int n){
//            n = m;
//            cout << "第" << n << "条线的嵌入存在未知错误！" << endl;
//            continue;
//        }
//        //cout << double(myTime_2) << endl;
//    }
//    //end = clock();
//    //cout << double(end - start) << endl; 
//    //cout << double(myTime_2) << endl;
//    cout << double(myTime) << endl;
//    //writeTXT(constrained_line_set, "constrained_line_set.txt");
//    //writeTXT(successful, "successful.txt");
//    //saveMeshOBJ(new_point, face, sign, "ours_mesh.obj");
//}



// MeshOptimization.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

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
        cout << "文件打开异常！" << endl;
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
        cout << "文件打开异常！" << endl;
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
        cout << "用法: " << argv[0]
            << " <line_obj> <mesh_obj> <out_dir> <para_1> <para_2> <node_interval>\n";
        cout << "示例: " << argv[0]
            << " remeshing.exe line.obj mesh_5000.obj G:\\MeshOptimization\\ISPRS\\ISPRS_data\\mesh_simplification 30 0.7 0.05\n";
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

    cout << "输入参数:\n"
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

    //开始计时
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

    //构建k-d树
    pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(new_point);

    MyMatrixXf successful = MyMatrixXf::Zero(line_num, 1);
    //开始嵌入线
    for (int m = 0; m < line_num; m++) {
        try {
            //构建点面邻接表
            cout << m << endl;
            //if (m == 289)
            //    continue;
            globalTest = m;
            dp_count = 0;

            MyMatrixXf PFneighbor;
            createPFneighbor(new_point, face, PFneighbor);

            //计算knn支撑域
            MyMatrixXf knn_face_idx = getKnnFace(tree, PFneighbor, face, line_set.block(2 * m, 0, 2, 3), sign, k_1, sample_interval_1);
            //writeTXT(knn_face_idx, "F:\\test\\knn.txt");
            //string filename_8 = "F:\\test\\region_1\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_8.c_str());
            knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
            //string filename_6 = "F:\\test\\region_3\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_6.c_str());
            knn_face_idx = regionConnect(new_point, face, knn_face_idx, PFneighbor);
            //string filename_7 = "F:\\test\\region_2\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ(new_point, face, knn_face_idx, sign, filename_7.c_str());
            //计算基准支撑域
            MyMatrixXf base_face_idx = isIntersectant(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1);
            //MyMatrixXf base_face_idx = isIntersectant2(new_point, face, knn_face_idx, line_set.block(2 * m, 0, 2, 3), sample_interval_1, dist_threshold);
            base_face_idx = holeFilling(base_face_idx, face, PFneighbor);
            base_face_idx = regionConnect(new_point, face, base_face_idx, PFneighbor);
            //writeTXT(sign, "F:\\test\\sign.txt");
            int base_face_num = base_face_idx.rows();
            cout << "region projected!" << endl;
            //string filename_5 = "F:\\test\\region_3\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ(new_point, face, base_face_idx, sign, filename_5.c_str());

            //计算扩展支撑域
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

            //start = clock();
            MyMatrixXf extend_face_idx = regionExtend(face, range_face_idx, base_face_idx, normal, l2f_dist, PFneighbor, para_1_set(m, 0), extend_threshold);
            //end = clock();
            //cout << double(end - start) << endl;
            extend_face_idx = shapeRepairing(new_point, face, extend_face_idx, PFneighbor);
            extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor);
            //extend_face_idx = holeFilling(extend_face_idx, face, PFneighbor, range_face_idx);
            //string filename_4 = "F:\\test\\region_2\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ2(new_point, face, extend_face_idx, sign, filename_4.c_str());
            //cout << "region extended!" << endl;
            //writeTXT(extend_face_idx, "F:\\test\\extend_test.txt");

            //start = clock();
            MyMatrixXf segment_face_idx = regionSegment(new_point, face, extend_face_idx, base_face_idx, line_set.block(2 * m, 0, 2, 3), constrained_line_set, sample_interval_1);
            //segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor, range_face_idx);
            segment_face_idx = holeFilling(segment_face_idx, face, PFneighbor);
            //string filename_3 = "F:\\test\\region_1\\" + to_string(m) + ".obj";
            //saveSubMeshOBJ2(new_point, face, segment_face_idx, sign, filename_3.c_str());
            //end = clock();
            //cout << double(end - start) << endl;
            cout << "region segmented!" << endl;
            if (haveConstrainedEdge(face, segment_face_idx, constrained_line_set))
                continue;
            //if (holeDetection(segment_face_idx, face, PFneighbor))
            //    continue;

            //string filename = "F:\\test\\region\\" + to_string(m) + ".obj";
            ////const char* file = nullptr;
            ////file = filename.c_str();
            //saveSubMeshOBJ(new_point, face, segment_face_idx, sign, filename.c_str());

            //MyMatrixXf segment_face_idx;
            //readTXT(segment_face_idx, "G:\\3D_reconstruction_test_weight\\simulate3\\1\\LSM\\36_1.txt");
            //segment_face_idx = segment_face_idx - MyMatrixXf::Ones(segment_face_idx.rows(),1);
            //writeTXT(segment_face_idx, "F:\\test\\segment_test.txt");
            MyMatrixXf node_idx;
            if (m == 0)
                node_idx = idx.block(0, 0, each_node_num(m, 0), 1);
            else
                node_idx = idx.block(each_node_num.block(0, 0, m, 1).sum(), 0, each_node_num(m, 0), 1);
            //start = clock();
            bool error = regionRemeshing(new_point, face, sign, segment_face_idx, node_idx, PFneighbor, constrained_line_set, 0.5, 0.5, m);
            if (error == false) {
                cout << "remeshing successfully!" << endl;
                successful(m, 0) = 1;
            }

            //writeTXT(constrained_line_set, "F:\\test\\constrained_line_set.txt");
            //end = clock();
            //cout << double(end - start) << endl << endl;

            //string out_file = output_dir + "\\" + mesh_file + "_out.obj";
            //saveMeshOBJ(new_point, face, sign, out_file.c_str());
            //cout << double(myTime) << endl;
        }
        catch (int n) {
            n = m;
            cout << "第" << n << "条线的嵌入存在未知错误！" << endl;
            continue;
        }
        //cout << double(myTime_2) << endl;
    }
    end = clock();
    cout << double(end - start) << endl; 
    //cout << double(myTime_2) << endl;
    cout << double(myTime) << endl;
    //writeTXT(constrained_line_set, "constrained_line_set.txt");
    //writeTXT(successful, "successful.txt");

    string out_file = output_dir + "\\" + mesh_file + "_out.obj";
    saveMeshOBJ(new_point, face, sign, out_file.c_str());
}





