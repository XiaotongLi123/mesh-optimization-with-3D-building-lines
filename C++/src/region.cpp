#include "region.h"
#include "lineProcess.h"
#include "myMath.h"
#include "dist.h"
#include "meshOper.h"
#include <vector>
#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "time.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/property_map/property_map.hpp>
#include "omp.h"

// Define the graph type
using namespace boost;
typedef adjacency_list<vecS, vecS, undirectedS, no_property,
    property<edge_weight_t, float>> Graph;
typedef property_map<Graph, edge_weight_t>::type WeightMap;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

// Custom breadth-first search visitor class
class ReachableVisitor : public default_bfs_visitor {
public:
    ReachableVisitor(std::vector<bool>& reachable) : reachable(reachable) {}

    // Called when visiting each vertex
    void examine_vertex(Vertex u, const Graph&) {
        reachable[u] = true;
    }

private:
    std::vector<bool>& reachable;
};

void faceNearestPointPos(const MyMatrixXf& proj_point, const MyMatrixXf& tri, MyMatrixXf& is_near);

float getVertexSumAngle(const MyMatrixXf& point, const MyMatrixXf& face, int this_vertex);

MyMatrixXf shortestPath(Graph G, vector<vector<size_t>>& path, size_t start_v, const MyMatrixXf& end_v);

void removeEdge(Graph& G, size_t start_v, size_t end_v);
void removeEdge(Graph& G, size_t v);

MyMatrixXf myMin(const MyMatrixXf& mat, MyMatrixXf& min_id, float error_rate);

void coMatClustering(MyMatrixXf& coMat);

bool segmentRegionGrowth(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& extend_region_constrained_edge, MyMatrixXf& region_segment);

MyMatrixXf extendSegmentRegion(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& base_region_face_idx, const MyMatrixXf& segment_region_face_idx, const MyMatrixXf& line, const MyMatrixXf& constrained_edge, float interval);

MyMatrixXf getKnnFace(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& PFneighbor, const MyMatrixXf& face, const MyMatrixXf& line, const MyMatrixXf& sign, int K, double interval) {
    MyMatrixXf node = getNode(line, interval);
    int face_num = face.rows();
    MyMatrixXf inface_idx = MyMatrixXf::Zero(999999, 1);
    MyMatrixXf face_label = MyMatrixXf::Zero(face_num, 1);
    int inface_count = 0;
    int r;
    vector<int> neigh_idx(K);
    for (int i = 0; i < node.rows(); i++) {
        neigh_idx = getKneighbor(kdtree, node.row(i), K, sign);
        for (int j = 0; j < K; j++) {
            for (int k = 0; k < PFneighbor(neigh_idx[j], 0); k++) {
                r = PFneighbor(neigh_idx[j], k + 1);
                if (face_label(r, 0) != 1) {
                    inface_idx(inface_count, 0) = PFneighbor(neigh_idx[j], k + 1);
                    face_label(r, 0) = 1;
                    inface_count += 1;
                }
            }
        }
    }
    return inface_idx.block(0, 0, inface_count, 1);
}

MyMatrixXf getBaseRegion(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const MyMatrixXf& PFneighbor, const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& line, const MyMatrixXf& sign, int K, double interval) {
    MyMatrixXf node = getNode(line, interval);
    int face_num = face.rows();
    MyMatrixXf inface_idx = MyMatrixXf::Zero(999999, 1);
    MyMatrixXf face_label = MyMatrixXf::Zero(face_num, 1);
    int inface_count = 0;
    MyMatrixXf neigh_face_idx = MyMatrixXf::Zero(999, 1);
    MyMatrixXf neigh_face_label = MyMatrixXf::Zero(face_num, 1);
    MyMatrixXf nearest_face_idx;
    int neigh_face_count;
    vector<int> neigh_idx(K);
    int r;
    for (int i = 0; i < node.rows(); i++) {
        neigh_face_count = 0;
        neigh_face_label = MyMatrixXf::Zero(face_num, 1);
        neigh_idx = getKneighbor(kdtree, node.row(i), K, sign);
        for (int j = 0; j < K; j++) {
            for (int k = 0; k < PFneighbor(neigh_idx[j], 0); k++) {
                r = PFneighbor(neigh_idx[j], k + 1);
                if (neigh_face_label(r, 0) != 1) {
                    neigh_face_idx(neigh_face_count++, 0) = PFneighbor(neigh_idx[j], k + 1);
                    neigh_face_label(r, 0) = 1;
                }
            }
        }
        point2meshDist(point, face, neigh_face_idx.block(0, 0, neigh_face_count, 1), node.row(i), nearest_face_idx);
        r = nearest_face_idx(0, 0); r = neigh_face_idx(r, 0);
        if (face_label(r, 0) != 1) {
            r = nearest_face_idx(0, 0);
            inface_idx(inface_count++, 0) = neigh_face_idx(r, 0);
            r = neigh_face_idx(r, 0);
            face_label(r, 0) = 1;
        }
    }
    return inface_idx.block(0, 0, inface_count, 1);
}

MyMatrixXf isIntersectant(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& knnsearch_face_idx, const MyMatrixXf& line, double interval) {
    MyMatrixXf node = getNode(line, interval);
    MyMatrixXf knn_face = getSubMat_Rows(face, knnsearch_face_idx);
    MyMatrixXf is_region_face = MyMatrixXf::Zero(knn_face.rows(), 1);
    MyMatrixXf nearest_point, nearest_face_idx;
    MyMatrixXf neigh_face_id, near_pt;
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    MyMatrixXf is_near = MyMatrixXf::Zero(3, 1);
    int r;
    for (int n = 0; n < node.rows(); n++) { 
        point2meshDist(point, face, knnsearch_face_idx, node.row(n), nearest_face_idx);
        for (int i = 0; i < nearest_face_idx.rows(); i++) {
            r = nearest_face_idx(i, 0);
            tri = getSubMat_Rows(point, face.row(knnsearch_face_idx(r, 0)));
            point2faceDist(node.row(n), tri, nearest_point);
            faceNearestPointPos(nearest_point, tri, is_near);
            int near_edge_num = is_near.sum();
            switch (near_edge_num) {
            case 0: {
                is_region_face(r, 0) = 1;
                break;
            }
            case 1: {
                neigh_face_id = findEdgeNeighFace(knn_face, getSubMat_Cols(face.row(knnsearch_face_idx(r, 0)), getIndex(is_near, 0, 0)));
                valueSubMat(is_region_face, neigh_face_id, MyMatrixXf::Zero(1, 1), 1);
                break;
            }
            case 2: {
                near_pt = getSubMat_Cols(face.row(knnsearch_face_idx(r, 0)), getIndex(is_near, 0, 0));
                neigh_face_id = findPointNeighFace(knn_face, near_pt(0, 0));
                valueSubMat(is_region_face, neigh_face_id, MyMatrixXf::Zero(1, 1), 1);
                break;
            }
            }
        }
    }

    return getSubMat_Rows(knnsearch_face_idx, getIndex(is_region_face, 1, 0));
}

MyMatrixXf isIntersectant2(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& knnsearch_face_idx, const MyMatrixXf& line, double interval, double dist_threshold) {
    MyMatrixXf node = getNode(line, interval);
    MyMatrixXf knn_face = getSubMat_Rows(face, knnsearch_face_idx);
    MyMatrixXf is_region_face = MyMatrixXf::Zero(knn_face.rows(), 1);
    MyMatrixXf nearest_point, nearest_face_idx;
    MyMatrixXf neigh_face_id, near_pt;
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    MyMatrixXf is_near = MyMatrixXf::Zero(3, 1);
    int r;
    for (int n = 0; n < node.rows(); n++) {
        point2meshDist(point, face, knnsearch_face_idx, node.row(n), nearest_face_idx);
        for (int i = 0; i < nearest_face_idx.rows(); i++) {
            r = nearest_face_idx(i, 0);
            tri = getSubMat_Rows(point, face.row(knnsearch_face_idx(r, 0)));
            point2faceDist(node.row(n), tri, nearest_point);
            faceNearestPointPos(nearest_point, tri, is_near);
            int near_edge_num = is_near.sum();
            switch (near_edge_num) {
            case 0: {
                is_region_face(r, 0) = 1;
                break;
            }
            case 1: {
                neigh_face_id = findEdgeNeighFace(knn_face, getSubMat_Cols(face.row(knnsearch_face_idx(r, 0)), getIndex(is_near, 0, 0)));
                valueSubMat(is_region_face, neigh_face_id, MyMatrixXf::Zero(1, 1), 1);
                break;
            }
            case 2: {
                near_pt = getSubMat_Cols(face.row(knnsearch_face_idx(r, 0)), getIndex(is_near, 0, 0));
                neigh_face_id = findPointNeighFace(knn_face, near_pt(0, 0));
                valueSubMat(is_region_face, neigh_face_id, MyMatrixXf::Zero(1, 1), 1);
                break;
            }
            }
        }
    }

    for (int i = 0; i < knn_face.rows(); i++) {
        if (is_region_face(i, 0) == 1)
            continue;
        float dist = line2faceDist2(getSubMat_Rows(point, knn_face.row(i)), line);
        if (dist <= dist_threshold) {
            is_region_face(i, 0) = 1;
        }
    }
    return getSubMat_Rows(knnsearch_face_idx, getIndex(is_region_face, 1, 0));
}

void faceNearestPointPos(const MyMatrixXf& proj_point, const MyMatrixXf& tri,MyMatrixXf& is_near) {
    MyMatrixXf line = MyMatrixXf::Zero(6, 3);
    line << tri.row(1),
        tri.row(2),
        tri.row(2),
        tri.row(0),
        tri.row(0),
        tri.row(1);
    float dist_threshold = 0.0;
    is_near = MyMatrixXf::Zero(3, 1);
    for (int i = 0; i < 3; i++) {
        dist_threshold = 0.25 * point2lineDist(tri.row(i), line.block(2 * i, 0, 2, 3));
        if (point2lineDist(proj_point, line.block(2 * i, 0, 2, 3)) <= dist_threshold) is_near(i, 0) = 1;
    }
}

MyMatrixXf holeFilling(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor, int max_loop) {
    // Hole filling algorithm
    // Step 1: Grow until there are no holes inside
    // Step 2: Fill holes with non-maximum connected components
    MyMatrixXf region_point, region_edge;
    getRegionPE(getSubMat_Rows(face, region_face_idx), region_point, region_edge);
    // Region growing
    int loop = 0;
    MyMatrixXf enlarge_face_idx, connect_outline;
    MyMatrixXf grow_face_idx = MyMatrixXf::Zero(face.rows(), 1);
    while (loop < max_loop) {
        enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
        getRegionPE(getSubMat_Rows(face, enlarge_face_idx), region_point, region_edge);
        connect_outline = findConnectEdge(getSubMat_Rows(region_edge.block(0, 0, region_edge.rows(), 2), getIndex(region_edge.col(2), 0, 0)));
        if (connect_outline.rows() == 1)
            break;
        loop++;
    }
    if (loop == max_loop) {
        cout << "Non-manifold or other errors exist in the mesh, unable to perform hole filling!" << endl;
        return region_face_idx;
    }
    // Extract the largest connected component
    grow_face_idx.block(0, 0, enlarge_face_idx.rows(), 1) = MyMatrixXf::Zero(enlarge_face_idx.rows(), 1);
    int grow_face_num = 0;
    for (int i = 0; i < enlarge_face_idx.rows(); i++) {
        if (!anyIsmember(region_face_idx, enlarge_face_idx(i, 0))) {
            grow_face_idx(grow_face_num++, 0) = enlarge_face_idx(i, 0);
        }
    }
    MyMatrixXf grow_connect_region = findConnectRegion(face, grow_face_idx.block(0, 0, grow_face_num, 1));
    MyMatrixXf max_grow_connect_id;
    getMax(grow_connect_region.col(0), max_grow_connect_id);
    int region_face_num = region_face_idx.rows();
    MyMatrixXf repaired_region_face_idx = MyMatrixXf::Zero(face.rows(), 1);
    repaired_region_face_idx.block(0, 0, region_face_idx.rows(), 1) << region_face_idx;
    for (int i = 0; i < grow_connect_region.rows(); i++) {
        if (i == max_grow_connect_id(0, 0))continue;
        repaired_region_face_idx.block(region_face_num, 0, grow_connect_region(i, 0), 1) = grow_connect_region.block(i, 1, 1, grow_connect_region(i, 0)).transpose();
        region_face_num += grow_connect_region(i, 0);
    }
    return repaired_region_face_idx.block(0, 0, region_face_num, 1);
}
MyMatrixXf holeFilling(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor, const MyMatrixXf& knn_face_idx) {
    // Hole filling algorithm
    // Step 1: Extract the boundary of knn_region
    // Step 2: Perform connectivity analysis on the region to be grown
    // Step 3: If the boundary of each connected region does not contain the boundary of knn_region, add it in
    MyMatrixXf region_point, region_edge;
    getRegionPE(getSubMat_Rows(face, region_face_idx), region_point, region_edge);
    MyMatrixXf connect_outline = findConnectEdge(getSubMat_Rows(region_edge.block(0, 0, region_edge.rows(), 2), getIndex(region_edge.col(2), 0, 0)));
    if (connect_outline.rows() == 1)
        return region_face_idx;
    // Step 1: Extract the boundary of knn_region
    MyMatrixXf knn_face = getSubMat_Rows(face, knn_face_idx);
    MyMatrixXf knn_outline, knn_outline_point;
    getOutlineSort(knn_face, knn_outline_point, knn_outline);

    // Step 2: Perform connectivity analysis on the region to be grown
    MyMatrixXf grow_face_idx = MyDifference(knn_face_idx, region_face_idx);
    MyMatrixXf grow_connect_region = findConnectRegion(face, grow_face_idx);
    if (grow_connect_region.rows() == 1) {
        return region_face_idx;
    }
    // Step 3: If the boundary of each connected region does not contain the boundary of knn_region, add it in
    int region_face_num = region_face_idx.rows();
    MyMatrixXf repaired_region_face_idx = MyMatrixXf::Zero(face.rows(), 1);
    repaired_region_face_idx.block(0, 0, region_face_idx.rows(), 1) << region_face_idx;
    int sign = 1;
    for (int i = 0; i < grow_connect_region.rows(); i++) {
        sign = 1;
        MyMatrixXf grow_outline_point, grow_outline;
        MyMatrixXf this_grow_face_idx = grow_connect_region.block(i, 1, 1, grow_connect_region(i, 0)).transpose();
        getOutlineSort(getSubMat_Rows(face, this_grow_face_idx), grow_outline_point, grow_outline);
        for (int j = 0; j < grow_outline.rows(); j++) {
            if (isExistLine(knn_outline, grow_outline.row(j))) {
                sign = 0; 
                break;
            }
        }
        if (sign == 0)continue;
        repaired_region_face_idx.block(region_face_num, 0, grow_connect_region(i, 0), 1) = grow_connect_region.block(i, 1, 1, grow_connect_region(i, 0)).transpose();
        region_face_num += grow_connect_region(i, 0);
    }
    return repaired_region_face_idx.block(0, 0, region_face_num, 1);
}

bool holeDetection(const MyMatrixXf& region_face_idx, const MyMatrixXf& face, const MyMatrixXf& PFneighbor) {
    // Hole filling algorithm
    // Step 1: Grow until there are no holes inside
    // Step 2: Fill holes with non-maximum connected components
    MyMatrixXf region_point, region_edge;
    getRegionPE(getSubMat_Rows(face, region_face_idx), region_point, region_edge);
    // Region growing
    int loop = 0;
    MyMatrixXf enlarge_face_idx, connect_outline;
    enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
    getRegionPE(getSubMat_Rows(face, enlarge_face_idx), region_point, region_edge);
    connect_outline = findConnectEdge(getSubMat_Rows(region_edge.block(0, 0, region_edge.rows(), 2), getIndex(region_edge.col(2), 0, 0)));
    if (connect_outline.rows() == 1)
        return false;
    else
        return true;

}

MyMatrixXf shapeRepairing(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor) {
    // Shape repairing
    // One region growing: check if all three points of the grown triangle are in the current region, if so, add the triangle
    // Loop until no more triangles meet the condition
    // Loop: check if the sum of the interior angles of all triangles adjacent to boundary points is less than 240бу, otherwise expand the vertex
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);
    // Region growing
    int sign = 1;
    int loop = 0;
    int count = region_face_idx.rows();
    MyMatrixXf repaired_region_face_idx = MyMatrixXf::Zero(face.rows(), 1);
    repaired_region_face_idx.block(0, 0, region_face_idx.rows(), 1) << region_face_idx;
    MyMatrixXf enlarge_face_idx, repaired_region_face;
    MyMatrixXf grow_face_idx = MyMatrixXf::Zero(face.rows(), 1);
    MyMatrixXf is_extend = MyMatrixXf::Zero(face.rows(), 1);
    MyMatrixXf this_face = MyMatrixXf::Zero(1, 3);
    int r;

    while (sign == 1 && loop < 30) {
        enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
        //cout << enlarge_face_idx << endl << endl;

        // Extract the grown faces
        grow_face_idx.block(0,0, enlarge_face_idx.rows(),1) << MyMatrixXf::Zero(enlarge_face_idx.rows(), 1);
        int grow_face_num = 0;
        for (int i = 0; i < enlarge_face_idx.rows(); i++) {
            if (!anyIsmember(repaired_region_face_idx.block(0, 0, count, 1), enlarge_face_idx(i, 0))) {
                grow_face_idx(grow_face_num++, 0) = enlarge_face_idx(i, 0);
            }
        }
        // Determine whether each triangle meets the condition
        is_extend.block(0, 0, grow_face_num, 1) = MyMatrixXf::Zero(grow_face_num, 1);
        for (int i = 0; i < grow_face_num; i++) {
            r = grow_face_idx(i, 0);
            this_face << face.row(r);
            if (anyIsmember(region_point.col(0), this_face(0, 0)) && anyIsmember(region_point.col(0), this_face(0, 1)) && anyIsmember(region_point.col(0), this_face(0, 2)))
                is_extend(i, 0) = 1;
        }

        // Growing
        loop++;
        if (is_extend.block(0, 0, grow_face_num, 1).sum() == 0)
            sign = 0;
        else {
            repaired_region_face_idx.block(count, 0, is_extend.block(0, 0, grow_face_num, 1).sum(), 1) = getSubMat_Rows(grow_face_idx.block(0, 0, grow_face_num, 1), getIndex(is_extend.block(0, 0, grow_face_num, 1), 1, 0));
            count += is_extend.block(0, 0, grow_face_num, 1).sum();
            repaired_region_face = getSubMat_Rows(face, repaired_region_face_idx.block(0, 0, count, 1));
            getRegionPE(repaired_region_face, region_point, region_edge);
        }
    }
    MyMatrixXf new_repaired_region_face_idx = repaired_region_face_idx.block(0, 0, count, 1);

    // Vertex angle judgment
    sign = 1;
    loop = 0;
    int old_rows;
    MyMatrixXf this_neigh_face_idx, in_this_neigh_face_idx, in_this_neigh_face, outline_point;
    while (sign == 1 && loop < 30) {
        old_rows = new_repaired_region_face_idx.rows();
        // Extract the boundary points of the region
        getRegionPE(getSubMat_Rows(face, new_repaired_region_face_idx), region_point, region_edge);
        outline_point = getSubMat_Rows(region_point.col(0), getIndex(region_point.col(1), 0, 0));
        int outline_point_num = outline_point.rows();
        for (int i = 0; i < outline_point_num; i++) {
            // Extract adjacent faces
            r = outline_point(i, 0);
            this_neigh_face_idx = PFneighbor.block(r, 1, 1, PFneighbor(r, 0));
            // Find adjacent faces within the region
            in_this_neigh_face_idx = MyIntersection(this_neigh_face_idx, new_repaired_region_face_idx);
            // Calculate the sum of angles
            in_this_neigh_face = getSubMat_Rows(face, in_this_neigh_face_idx);
            float sum_angle = getVertexSumAngle(point, in_this_neigh_face, r);
            if (sum_angle > 240) {
                new_repaired_region_face_idx = MyUnion(this_neigh_face_idx, new_repaired_region_face_idx);
            }     
        }
        if (old_rows == new_repaired_region_face_idx.rows())
            sign = 0;
        loop++;
    }
    return new_repaired_region_face_idx;
}

float getVertexSumAngle(const MyMatrixXf& point, const MyMatrixXf& face, int this_vertex) {
    int face_num = face.rows();
    float sum_angle = 0;
    int r;
    MyMatrixXf not_vertex;
    MyMatrixXf vec_1 = MyMatrixXf::Zero(1, 3);
    MyMatrixXf vec_2 = MyMatrixXf::Zero(1, 3);
    for (int i = 0; i < face_num; i++) {
        not_vertex = getSubMat_Cols(face.row(i), getIndex(face.row(i), this_vertex, 1));
        r = not_vertex(0, 0);
        vec_1 = point.row(r) - point.row(this_vertex);
        r = not_vertex(0, 1);
        vec_2 = point.row(r) - point.row(this_vertex);
        sum_angle += getVecAngle(vec_1, vec_2);
    }
    return sum_angle * 180 / M_PI;
}

MyMatrixXf regionGrowth(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor, int k) {
    MyMatrixXf grow_face_idx = region_face_idx;
    MyMatrixXf region_point, region_edge;
    for (int i = 0; i < k; i++) {
        getRegionPE(getSubMat_Rows(face, grow_face_idx), region_point, region_edge);
        grow_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
    }
    return grow_face_idx;
}

MyMatrixXf regionExtend(const MyMatrixXf& face, const MyMatrixXf& knn_face_idx, const MyMatrixXf& region_face_idx, const MyMatrixXf& normal, const MyMatrixXf& dist, const MyMatrixXf& PFneighbor, float threshold_1, float threshold_2, int max_loop) {
    // Extension of the supporting region of the line segment
    // dist is the distance from each face in knn_face to the constrained line segment
    int r, r_1;
    int sign = 1;
    MyMatrixXf knn_face = getSubMat_Rows(face, knn_face_idx);
    int loop = 0;
    int extend_face_num = region_face_idx.rows();
    MyMatrixXf extend_face_idx = MyMatrixXf::Zero(knn_face_idx.rows(), 1);
    extend_face_idx.block(0, 0, region_face_idx.rows(), 1) << region_face_idx;
    MyMatrixXf region_point, region_edge;
    MyMatrixXf extend_face, enlarge_face_idx, grow_face_idx, grow_face;
    
    MyMatrixXf edge = MyMatrixXf::Zero(3, 2);
    MyMatrixXf index = MyMatrixXf::Zero(1, 2);
    MyMatrixXf index_1 = MyMatrixXf::Zero(1, 2);
    MyMatrixXf index_2 = MyMatrixXf::Zero(1, 2);
    MyMatrixXf neigh_region_face_idx, neigh_face_id, neigh_face_idx, temp;

    MyMatrixXf is_extend = MyMatrixXf::Zero(knn_face_idx.rows(), 1);// Determine whether to extend
    MyMatrixXf angle = MyMatrixXf::Zero(knn_face_idx.rows(), 1);// Angle between adjacent faces of the growing edge
    MyMatrixXf neigh_region_face_idx_ori = MyMatrixXf::Zero(99999, 1);
    while (sign == 1 && loop < max_loop) {
        extend_face = getSubMat_Rows(face, extend_face_idx.block(0, 0, extend_face_num, 1));
        getRegionPE(extend_face, region_point, region_edge);

        enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
        enlarge_face_idx = MyIntersection(enlarge_face_idx, knn_face_idx);
        grow_face_idx = MyDifference(enlarge_face_idx, extend_face_idx.block(0, 0, extend_face_num, 1));
        grow_face = getSubMat_Rows(face, grow_face_idx);

        // Determine whether each face can be regionally expanded
        int grow_face_num = grow_face_idx.rows();
        is_extend.block(0, 0, grow_face_num, 1) = MyMatrixXf::Zero(grow_face_num, 1);

        for (int i = 0; i < grow_face_num; i++) {
            // Find all neighboring faces (adjacent edges) for each grown face
            
            edge << grow_face(i, 1), grow_face(i, 2),
                grow_face(i, 2), grow_face(i, 0),
                grow_face(i, 0), grow_face(i, 1);
            
            int count = 0;
            for (int j = 0; j < 3; j++) {
                neigh_face_id = findEdgeNeighFace(knn_face, edge.row(j));
                neigh_face_idx = getSubMat_Rows(knn_face_idx, neigh_face_id);
                neigh_face_idx = MyDifference(neigh_face_idx, grow_face_idx.row(i));
                temp = MyIntersection(neigh_face_idx, extend_face_idx.block(0, 0, extend_face_num, 1));
                neigh_region_face_idx_ori.block(count, 0, temp.rows(), 1) = temp;
                count += temp.rows();
            }
            neigh_region_face_idx = neigh_region_face_idx_ori.block(0, 0, count, 1);
            // Determine whether to grow based on the distance from the constrained segment to the current face and neighboring faces in the supporting region
            if (neigh_region_face_idx.rows() == 0 || neigh_region_face_idx.cols() == 0) continue;
            for (int j = 0; j < neigh_region_face_idx.rows(); j++) {
                index_1 = getIndex(knn_face_idx, neigh_region_face_idx(j, 0), 0);
                index_2 = getIndex(knn_face_idx, grow_face_idx(i, 0), 0);
                r = index_1(0, 0); r_1 = index_2(0, 0);
                angle(j, 0) = getVecAngle(normal.row(r), normal.row(r_1)) * 180 / M_PI;
            }
            float min_angle = angle.block(0, 0, neigh_region_face_idx.rows(), 1).minCoeff();
            index = getIndex(knn_face_idx, grow_face_idx(i, 0), 0);
            r = index(0, 0);
            if (min_angle < threshold_1 && dist(r, 0) < threshold_2)
                is_extend(i, 0) = 1;
        }
        int this_extend_face_num = is_extend.block(0, 0, grow_face_num, 1).sum();
        if (this_extend_face_num == 0)
            sign = 0;
        extend_face_idx.block(extend_face_num, 0, this_extend_face_num, 1) = getSubMat_Rows(grow_face_idx, getIndex(is_extend.block(0, 0, grow_face_num, 1), 1, 0));
        extend_face_num += this_extend_face_num;

        loop++;
    }
    if (loop == max_loop) {
        cout << "Reached the maximum number of iterations!" << endl;
    }
    return extend_face_idx.block(0, 0, extend_face_num, 1);
}


MyMatrixXf regionConnect(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& PFneighbor, int max_loop) {
    // Multi-connected supporting region connection algorithm
    // First, region growing until the supporting region is connected
    // Then, build the connectivity graph of the grown region and grow along the shortest path
    int r;
    MyMatrixXf ori_connect_region = findConnectRegion(face, region_face_idx);
    
    int ori_connect_region_count = ori_connect_region.rows();
    int connect_region_count = ori_connect_region_count;
    if (connect_region_count == 1 || connect_region_count == 0)
        return region_face_idx;

    // Region growing
    int loop = 0;
    MyMatrixXf enlarge_face_idx = region_face_idx;
    MyMatrixXf region_point, region_edge;
    MyMatrixXf enlarge_region_face, connect_region;
    while (connect_region_count != 1 && loop < max_loop) {
        enlarge_region_face = getSubMat_Rows(face, enlarge_face_idx);
        getRegionPE(enlarge_region_face, region_point, region_edge);
        enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
        connect_region = findConnectRegion(face, enlarge_face_idx);
        connect_region_count = connect_region.rows();
        loop++;
    }
    if (loop == max_loop) {
        cout << "Reached the maximum number of iterations!" << endl;
        return region_face_idx;
    }
    MyMatrixXf new_connect_region = MyMatrixXf::Zero(ori_connect_region_count, enlarge_face_idx.rows() + 1);
    new_connect_region.block(0, 0, ori_connect_region.rows(), ori_connect_region.cols()) = ori_connect_region;
    // Shortest path growing
    // Build the adjacency matrix and undirected graph of region_point, using edge length as weight
    int pts_num = point.rows();
    enlarge_region_face = getSubMat_Rows(face, enlarge_face_idx);
    getRegionPE(enlarge_region_face, region_point, region_edge);
    int region_edge_num = region_edge.rows();
    Graph G;
    for (int i = 0; i < region_edge_num; i++) {
        float weight = point2pointDist(point.row(region_edge(i, 0)), point.row(region_edge(i, 1)));
        add_edge(region_edge(i, 0), region_edge(i, 1), weight, G);
    }

    // Shortest path search
    // 1. Search for the shortest path between each pair of connected components
    // 2. In each iteration, connect the two closest connected components and update the component table
    loop = 0;
    while (ori_connect_region_count != 1 && loop < max_loop) {
        // Extract the boundary points of each connected component
        MyMatrixXf outline_all = MyMatrixXf::Zero(ori_connect_region_count, 1000);
        for (int i = 0; i < ori_connect_region_count; i++) {
            MyMatrixXf this_region_face_idx = new_connect_region.block(i, 1, 1, new_connect_region(i, 0));
            MyMatrixXf this_region_point, this_region_edge;
            getRegionPE(getSubMat_Rows(face, this_region_face_idx), this_region_point, this_region_edge);
            MyMatrixXf outline = getSubMat_Rows(this_region_point.col(0), getIndex(this_region_point.col(1), 0, 0));
            outline_all(i, 0) = outline.rows();
            outline_all.block(i, 1, 1, outline_all(i, 0)) = outline.transpose();
        }
        MyMatrixXf shortest_path = MyMatrixXf::Zero(ori_connect_region_count * ori_connect_region_count, 1 + pts_num);
        MyMatrixXf shortest_path_connect_region = MyMatrixXf::Zero(ori_connect_region_count * ori_connect_region_count, 2);
        // Shortest path from region to region
        for (int i = 0; i < ori_connect_region_count; i++) {
            for (int j = 0; j < ori_connect_region_count; j++) {
                if (i >= j)
                    continue;
                // Shortest path between boundary points of two regions
                MyMatrixXf this_shortest_path = MyMatrixXf::Zero(outline_all(i, 0) * outline_all(j, 0), min(1 + pts_num, 99999));
                MyMatrixXf d_all = MyMatrixXf::Zero(outline_all(i, 0) * outline_all(j, 0), 1);
                MyMatrixXf p_size = MyMatrixXf::Zero(outline_all(i, 0) * outline_all(j, 0), 1);
                for (int m = 0; m < outline_all(i, 0); m++) {
                    vector<vector<size_t>>path;
                    MyMatrixXf dist = shortestPath(G, path, outline_all(i, m + 1), outline_all.block(j, 1, 1, outline_all(j, 0)).transpose());

                    for (int n = 0; n < outline_all(j, 0); n++) {
                        MyMatrixXf p = vec2mat(path[n]);
                        r = outline_all(j, 0);
                        d_all(m * r + n, 0) = dist(n, 0);
                        this_shortest_path.block(m * r + n, 0, 1, p.rows()) = p.transpose();
                        p_size(m * r + n, 0) = p.rows();
                    }
                }
                MyMatrixXf min_id;
                getMin(d_all, min_id);
                r = min_id(0, 0);
                shortest_path(i * ori_connect_region_count + j, 0) = p_size(r, 0);
                shortest_path.block(i * ori_connect_region_count + j, 1, 1, p_size(r, 0)) = this_shortest_path.block(r, 0, 1, p_size(r, 0));
                shortest_path_connect_region(i * ori_connect_region_count + j, 0) = i;
                shortest_path_connect_region(i * ori_connect_region_count + j, 1) = j;
            }
        }
        for (int i = 0; i < ori_connect_region_count * ori_connect_region_count; i++) {
            if (shortest_path(i, 0) == 0)
                shortest_path(i, 0) = std::numeric_limits<float>::max();
        }
        MyMatrixXf min_id;
        float min_path_ptsnum = getMin(shortest_path.col(0), min_id);
        int r = min_id(0, 0);
        MyMatrixXf neigh_face_idx = getPointNeighFace(shortest_path.block(r, 1, 1, min_path_ptsnum).transpose(), PFneighbor);
        int start_region = shortest_path_connect_region(r, 0);
        int end_region = shortest_path_connect_region(r, 1);
        MyMatrixXf temp_union = MyUnion(new_connect_region.block(start_region, 1, 1, new_connect_region(start_region, 0)), new_connect_region.block(end_region, 1, 1, new_connect_region(end_region, 0)));
        MyMatrixXf new_face_idx = MyUnion(temp_union,neigh_face_idx);
        new_connect_region(start_region, 0) = new_face_idx.rows();
        new_connect_region.block(start_region, 1, 1, new_face_idx.rows()) = new_face_idx.transpose();
        new_connect_region.block(end_region, 0, new_connect_region.rows() - end_region - 1, new_connect_region.cols()) = new_connect_region.block(end_region + 1, 0, new_connect_region.rows() - end_region - 1, new_connect_region.cols());
        ori_connect_region_count -= 1;
    }
    if (loop == max_loop) {
        cout << "max_loop!" << endl;
        return region_face_idx;
    }
    return new_connect_region.block(0, 1, 1, new_connect_region(0, 0)).transpose();
}

MyMatrixXf shortestPath(Graph G, vector<vector<size_t>>& path, size_t start_v,const MyMatrixXf& end_v) {    
    // Get the weight property map
    WeightMap weight_map = get(edge_weight, G);
    // Vectors to hold distances and predecessors
    std::vector<float> distance(num_vertices(G), std::numeric_limits<float>::max());
    std::vector<Vertex> predecessor(num_vertices(G));
    // Run Dijkstra's algorithm from start_v

    dijkstra_shortest_paths(G, start_v,
        predecessor_map(make_iterator_property_map(predecessor.begin(), get(vertex_index, G))).
        distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, G))).weight_map(weight_map));
    // Output the results
    path.resize(end_v.rows());
    MyMatrixXf dist = MyMatrixXf::Zero(end_v.rows(), 1);
    for (std::size_t i = 0; i < end_v.rows(); ++i) {
        dist(i, 0) = distance[end_v(i, 0)];
        path[i].resize(num_vertices(G));
        int j = 0;
        // Trace the path from the source to vertex i
        if (dist(i, 0) != std::numeric_limits<float>::max()) {
            for (Vertex v = end_v(i, 0); v != start_v; v = predecessor[v]) {
                path[i][j++] = v;
            }
        }
        path[i][j++] = start_v;
        path[i].resize(j);
        std::reverse(path[i].begin(), path[i].begin());
    }
    return dist;
}

void removeEdge(Graph& G, size_t start_v, size_t end_v) {
    Edge e;
    bool exists;
    tie(e, exists) = edge(start_v, end_v, G);  // Get the edge between vertices 

    if (exists) {
        // Remove the edge
        remove_edge(e, G);
    }
}
void removeEdge(Graph& G, size_t v) {
    std::vector<Edge> edges_to_remove;
    std::pair<OutEdgeIterator, OutEdgeIterator> edge_pair = out_edges(v, G);
    for (OutEdgeIterator e = edge_pair.first; e != edge_pair.second; ++e) {
        edges_to_remove.push_back(*e);
    }

    // Remove all collected edges
    for (const auto& edge : edges_to_remove) {
        remove_edge(edge, G);
    }
}

MyMatrixXf myMin(const MyMatrixXf& mat, MyMatrixXf& min_id, float error_rate) {
    MyMatrixXf temp_id;
    float min_val = getMin(mat, temp_id);
    MyMatrixXf min_vals = MyMatrixXf::Zero(mat.rows() * mat.cols(), 1);
    MyMatrixXf min_id_ori = MyMatrixXf::Zero(mat.rows() * mat.cols(), 2);
    int count = 0;
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            if (mat(i, j) <= (1 + error_rate) * min_val) {
                min_vals(count, 0) = mat(i, j);
                min_id_ori.row(count++) << i, j;
            }
        }
    }
    min_id = min_id_ori.block(0, 0, count, 2);
    return min_vals.block(0, 0, count, 1);
}


void coMatClustering(MyMatrixXf& coMat) {
    Graph G;
    for (int i = 0; i < coMat.rows(); i++) {
        for (int j = i + 1; j < coMat.cols(); j++) {
            if (coMat(i, j) == 1) {
                add_edge(i, j, 1, G);
            }
        }
    }
    MyMatrixXf sign = MyMatrixXf::Zero(coMat.rows(), 1);
    MyMatrixXf cluster = MyMatrixXf::Zero(coMat.rows(), coMat.rows() + 1);
    int cluster_num = 0;
    for (int i = 0; i < coMat.rows(); i++) {
        if (sign(i, 0) == 1)
            continue;
        if (coMat.row(i).sum() == 0)
            continue;
        // Initialize reachability marker array
        std::vector<bool> reachable(num_vertices(G), false);

        // Perform breadth-first search starting from the start vertex
        ReachableVisitor vis(reachable);
        breadth_first_search(G,
            vertex(i, G),
            visitor(vis));

        // Output all vertices reachable from the start point
        int count = 0;
        for (Vertex v = 0; v < num_vertices(G); ++v) {
            if (reachable[v]) {
                cluster(cluster_num, (count++) + 1) = v;
                cluster(cluster_num, 0) += 1;
                sign(v, 0) = 1;
            }
        }
        cluster_num += 1;
    }
    for (int i = 0; i < cluster_num; i++) {
        for (int j = 1; j < cluster(i, 0) + 1; j++) {
            valueSubMat(coMat, getSubMat(cluster, i * MyMatrixXf::Ones(1, 1), j * MyMatrixXf::Ones(1, 1)), cluster.block(i, 1, 1, cluster(1, 0)), 1);
        }
    }
}

bool segmentRegionGrowth(const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& extend_region_constrained_edge, MyMatrixXf& region_segment) {
    bool error = false;
    int r;
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    int region_face_num = region_face_idx.rows();
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);
    int region_edge_num = region_edge.rows();
    MyMatrixXf outline_edge = getSubMat_Rows(region_edge.block(0, 0, region_edge_num, 2), getIndex(region_edge.col(2), 0, 0));
    int outline_edge_num = outline_edge.rows();

    // Initialize
    MyMatrixXf used_face = MyMatrixXf::Zero(region_face_num, 1);
    int extend_region_constrained_edge_num = extend_region_constrained_edge.rows();

    int loop = 0;
    MyMatrixXf region_segment_ori = MyMatrixXf::Zero(region_face_num, region_face_num + 1);
    int region_segment_count = 0;
    MyMatrixXf this_segment_region_face_idx = MyMatrixXf::Zero(region_face_num, 1);
    MyMatrixXf used_edge = MyMatrixXf::Zero(region_edge_num, 1);
    MyMatrixXf exist_line_id;
    MyMatrixXf candidate_face_idx, seed_edge;
    int seed_face_idx, count, loop_start, loop_end;
    MyMatrixXf this_face = MyMatrixXf::Zero(1, 3);
    MyMatrixXf this_face_edge = MyMatrixXf::Zero(3, 2);
    MyMatrixXf growth_edge = MyMatrixXf::Zero(3, 2);
    MyMatrixXf temp;
    MyMatrixXf idx_temp = MyMatrixXf::Zero(region_face_num, 1);
    for (int k = 0; k < region_face_num; k++)
        idx_temp(k, 0) = k;
    while (used_face.sum() != region_face_num && loop < 10000) {
        // Update the usage of edges in each growth iteration
        used_edge = MyMatrixXf::Zero(region_edge_num, 1);
        for (int i = 0; i < outline_edge_num; i++) {
            // Boundary edges cannot grow
            isExistLine(region_edge.block(0, 0, region_edge_num, 2), outline_edge.row(i), exist_line_id);
            for (int j = 0; j < exist_line_id.rows(); j++) {
                r = exist_line_id(j, 0);
                used_edge(r, 0) += 1;
            }
        }
        for (int i = 0; i < extend_region_constrained_edge_num; i++) {
            // Extended constrained edges cannot grow
            isExistLine(region_edge.block(0, 0, region_edge_num, 2), extend_region_constrained_edge.row(i), exist_line_id);
            for (int j = 0; j < exist_line_id.rows(); j++) {
                r = exist_line_id(j, 0);
                used_edge(r, 0) += 1;
            }
                
        }

        // Seed triangle face
        candidate_face_idx = getIndex(used_face, 0, 0);
        seed_face_idx = candidate_face_idx(0, 0);
        seed_edge = getFaceEdge(region_face.row(seed_face_idx));
        for (int i = 0; i < 3; i++) {
            isExistLine(region_edge.block(0, 0, region_edge_num, 2), seed_edge.row(i), exist_line_id);
            for (int j = 0; j < exist_line_id.rows(); j++) {
                r = exist_line_id(j, 0);
                used_edge(r, 0) += 1;
            }
        }

        this_segment_region_face_idx(0, 0) = seed_face_idx;
        count = 1;
        loop_start = count;
        loop_end = count;
        // Growth
        while (loop_start <= loop_end) {
            for (int i = loop_start; i < loop_end + 1; i++) {
                r = this_segment_region_face_idx(i - 1, 0);
                this_face = region_face.row(r);
                this_face_edge = getFaceEdge(this_face);
                for (int j = 0; j < 3; j++) {
                    isExistLine(region_edge.block(0, 0, region_edge_num, 2), this_face_edge.row(j), exist_line_id);
                    // If an edge has been used more than twice, do not grow further (to avoid duplication)
                    r = exist_line_id(0, 0);
                    if (used_edge(r, 0) >= 2)
                        continue;
                    // Growth
                    float third_face_id = findCertainNeighFace(region_face, idx_temp, this_face_edge.row(j), this_segment_region_face_idx(i - 1, 0));
                    if (!anyIsmember(this_segment_region_face_idx.block(0, 0, count, 1), third_face_id) && third_face_id != -1) {
                        this_segment_region_face_idx(count++, 0) = third_face_id;
                        growth_edge = getFaceEdge(region_face.row(third_face_id));
                        for (int k = 0; k < 3; k++) {
                            isExistLine(region_edge.block(0, 0, region_edge_num, 2), growth_edge.row(k), exist_line_id);
                            for (int m = 0; m < exist_line_id.rows(); m++) {
                                r = exist_line_id(m, 0);
                                used_edge(r, 0) += 1;
                            }
                        }
                    }
                }
                r = this_segment_region_face_idx(i - 1, 0);
                used_face(r, 0) = 1;
            }
            loop_start = loop_end + 1;
            loop_end = count;
            loop++;
        }
        temp = getSubMat_Rows(region_face_idx, this_segment_region_face_idx.block(0, 0, count, 1));
        region_segment_ori.block(region_segment_count++, 0, 1, count + 1) << count, temp.transpose();
    }
    region_segment = region_segment_ori.block(0, 0, region_segment_count, region_segment_ori.cols());
    if (loop == 10000)
        error = true;
    return error;
}

MyMatrixXf extendSegmentRegion(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& base_region_face_idx, const MyMatrixXf& segment_region_face_idx, const MyMatrixXf& line, const MyMatrixXf& constrained_edge, float interval) {
    MyMatrixXf extend_segment_region_face_idx = MyMatrixXf::Zero(region_face_idx.rows(), 1);
    extend_segment_region_face_idx.block(0, 0, segment_region_face_idx.rows(), 1) << segment_region_face_idx;
    int extend_segment_face_num = segment_region_face_idx.rows();
    float dist = line2meshMeanDist(point, face, segment_region_face_idx, line, interval);
    int sign = 1;
    int r;
    MyMatrixXf region_face, outline_point_sort, outline_sort;
    MyMatrixXf segment_region_constrained_edge_ori = MyMatrixXf::Zero(constrained_edge.rows(), 2);
    MyMatrixXf segment_region_constrained_edge;
    MyMatrixXf neigh_face_idx, temp_intersect, this_candidate_face_idx, add_face_id;
    MyMatrixXf candidate_face = MyMatrixXf::Zero(1, 3);
    MyMatrixXf candidate_face_edge = MyMatrixXf::Zero(3, 2);
    MyMatrixXf new_segment_region_face_idx = MyMatrixXf::Zero(region_face_idx.rows(), 1);
    bool is_exist_line_1, is_exist_line_2, is_exist_line_3;
    MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
    MyMatrixXf opt_id;
    while (sign == 1) {
        sign = 0;
        region_face = getSubMat_Rows(face, extend_segment_region_face_idx.block(0, 0, extend_segment_face_num, 1));
        getOutlineSort(region_face, outline_point_sort, outline_sort);
        //Extract constrained line segments on the boundary
        int count = 0;
        for (int i = 0; i < constrained_edge.rows(); i++) {
            if (isExistLine(outline_sort, constrained_edge.row(i)))
                segment_region_constrained_edge_ori.row(count++) = constrained_edge.row(i);
        }
        segment_region_constrained_edge = segment_region_constrained_edge_ori.block(0, 0, count, 2);

        int outline_num = outline_sort.rows();
        // Extract candidate faces
        MyMatrixXf candidate_face_idx = MyMatrixXf::Zero(outline_num, 1);
        MyMatrixXf this_new_dist = MyMatrixXf::Zero(outline_num, 1);
        int candidate_count = 0;
        for (int i = 0; i < outline_num; i++) {
            if (isExistLine(constrained_edge, outline_sort.row(i)))
                continue;// Cannot extend from a constrained edge
            // If there is no adjacent face, cannot extend
            
            neigh_face_idx = findEdgeNeighFace(face, outline_sort.row(i));
            temp_intersect = MyIntersection(neigh_face_idx, region_face_idx);
            this_candidate_face_idx = MyDifference(temp_intersect, extend_segment_region_face_idx.block(0, 0, extend_segment_face_num, 1));
            if (this_candidate_face_idx.rows() == 0)
                continue;
            // After extension, cannot include the constrained edge again
            int is_extend = 1;
            r = this_candidate_face_idx(0, 0);
            candidate_face = face.row(r);
            candidate_face_edge = getFaceEdge(candidate_face);
            for (int j = 0; j < 3; j++) {
                if (isExistLine(segment_region_constrained_edge, candidate_face_edge.row(j))) {
                    is_extend = 0;
                    break;
                }
            }
            if (is_extend == 0)
                continue;
            // After extension, if the distance does not change, holes are not allowed; if the distance decreases, holes are allowed
            if (anyIsmember(base_region_face_idx, this_candidate_face_idx)) {
                new_segment_region_face_idx.block(0, 0, extend_segment_face_num + 1, 1) << extend_segment_region_face_idx.block(0, 0, extend_segment_face_num, 1),
                    this_candidate_face_idx;
                float temp_new_dist = line2meshMeanDist(point, face, new_segment_region_face_idx.block(0, 0, extend_segment_face_num + 1, 1), line, interval);
                if (temp_new_dist > dist)
                    is_extend = 0;
                if (temp_new_dist == dist) {
                    int count = 0;
                    for (int j = 0; j < 3; j++) {
                        if (anyIsmember(outline_point_sort, candidate_face(0, j)))
                            count++;
                    }
                    if (count == 3) {
                        // If all three vertices of the new face are on the boundary, and the two new edges are not on the original boundary, a hole will be created. Reject this case.
                        is_exist_line_1 = isExistLine(outline_sort, candidate_face_edge.row(0));
                        is_exist_line_2 = isExistLine(outline_sort, candidate_face_edge.row(1));
                        is_exist_line_3 = isExistLine(outline_sort, candidate_face_edge.row(2));
                        if (is_exist_line_1 + is_exist_line_2 + is_exist_line_3 == 1)
                            is_extend = 0;
                    }
                }
                if (is_extend == 1) {
                    this_new_dist(candidate_count, 0) = temp_new_dist;
                    candidate_face_idx(candidate_count++, 0) = this_candidate_face_idx(0, 0);
                }
            }
        }
        if (candidate_count != 0) {
            sign = 1;
            MyMatrixXf candidate_dist = MyMatrixXf::Zero(candidate_count, 1);
            // Only extend one face at a time to prevent two faces from simultaneously including the constrained edge
            for (int i = 0; i < candidate_count; i++) {
                tri = getSubMat_Rows(point, face.row(candidate_face_idx(i, 0)));
                candidate_dist(i, 0) = 0.5 * line2faceDist(tri, line, interval) + 0.5 * line2faceDist2(tri, line);
            }
            getMin(candidate_dist, opt_id);
            r = opt_id(0, 0);
            dist = this_new_dist(r, 0);
            extend_segment_region_face_idx(extend_segment_face_num++, 0) = candidate_face_idx(r, 0);
        }
    }
    return extend_segment_region_face_idx.block(0, 0, extend_segment_face_num, 1);
}


MyMatrixXf regionSegment(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx, const MyMatrixXf& base_region_face_idx, const MyMatrixXf& line, const MyMatrixXf& constrained_edge, float interval) {
    int r;
    int pts_num = point.rows();
    MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
    MyMatrixXf region_point, region_edge;
    getRegionPE(region_face, region_point, region_edge);
    int region_point_num = region_point.rows();
    int region_edge_num = region_edge.rows();
    int constrained_edge_num = constrained_edge.rows();

    if (constrained_edge_num == 0)
        return region_face_idx;

    // Extract constrained edges inside the supporting region (excluding the boundary) and connect them
    MyMatrixXf region_constrained_edge_ori = MyMatrixXf::Zero(constrained_edge_num, 2);
    int count = 0;
    MyMatrixXf in_region_edge = getSubMat_Rows(region_edge.block(0, 0, region_edge.rows(), 2), getIndex(region_edge.col(2), 1, 0));
    for (int i = 0; i < constrained_edge_num; i++) {
        if (isExistLine(in_region_edge, constrained_edge.row(i))) {
            region_constrained_edge_ori.row(count++) = constrained_edge.row(i);
        }
    }
    MyMatrixXf region_constrained_edge = region_constrained_edge_ori.block(0, 0, count, 2);
    int region_constrained_edge_num = region_constrained_edge.rows();

    if (region_constrained_edge_num == 0)
        return region_face_idx;

    MyMatrixXf connect_ori = findConnectEdge(region_constrained_edge);
    int connect_constrained_edge_ori_num = connect_ori.rows();
    MyMatrixXf connect_region_constrained_edge_ori = MyMatrixXf::Zero(100, 1000);
    connect_region_constrained_edge_ori.block(0, 0, connect_ori.rows(), connect_ori.cols()) = connect_ori;

    // Extract points that are forbidden to extend (points on constrained_edge)
    MyMatrixXf ban_points = MyMatrixXf::Zero(9999, 1);
    int ban_point_num = 0;
    for (int i = 0; i < connect_ori.rows(); i++) {
        ban_points.block(ban_point_num, 0, connect_ori(i, 0), 1) = connect_ori.block(i, 1, 1, connect_ori(i, 0)).transpose();
        ban_point_num += connect_ori(i, 0);
    }
    // Remove loops with heads
    MyMatrixXf idx;
    for (int i = 0; i < connect_constrained_edge_ori_num; i++) {
        for (int j = 1; j < connect_region_constrained_edge_ori(i, 0) + 1; j++) {
            if (anyIsmember(connect_region_constrained_edge_ori.block(i, 1, 1, j - 1), connect_region_constrained_edge_ori(i, j))) {
                idx = getIndex(connect_region_constrained_edge_ori.block(i, 1, 1, connect_region_constrained_edge_ori(i, 0)), connect_region_constrained_edge_ori(i, j), 0);
                if (idx(0, 1) == connect_region_constrained_edge_ori(i, 0) - 1 && idx(0, 0) != 0) {
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0) = idx(0, 1) - idx(0, 0) + 1;
                    connect_region_constrained_edge_ori.block(connect_constrained_edge_ori_num, 1, 1, connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0)) =
                        connect_region_constrained_edge_ori.block(i, idx(0, 0) + 1, 1, connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0));
                    connect_region_constrained_edge_ori(i, 0) = connect_region_constrained_edge_ori(i, 0) - connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0) + 1;
                    connect_region_constrained_edge_ori.block(i, idx(0, 0) + 2, 1, idx(0, 1) - idx(0, 0)) = MyMatrixXf::Zero(1, idx(0, 1) - idx(0, 0));
                    connect_constrained_edge_ori_num += 1;
                }
                if (idx(0, 1) != connect_region_constrained_edge_ori(i, 0) - 1 && idx(0, 0) == 0) {
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0) = connect_region_constrained_edge_ori(i, 0) - (idx(0, 1) - idx(0, 0));
                    connect_region_constrained_edge_ori.block(connect_constrained_edge_ori_num, 1, 1, connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0)) =
                        connect_region_constrained_edge_ori.block(i, idx(0, 1) + 1, 1, connect_region_constrained_edge_ori(connect_constrained_edge_ori_num, 0));
                    connect_region_constrained_edge_ori.block(i, idx(0, 1) + 2, 1, connect_region_constrained_edge_ori(i, 0) - 1 - idx(0, 1)) = MyMatrixXf::Zero(1, connect_region_constrained_edge_ori(i, 0) - 1 - idx(0, 1));
                    connect_region_constrained_edge_ori(i, 0) = idx(0, 1) - idx(0, 0) + 1;
                    connect_constrained_edge_ori_num += 1;
                }
                break;
            }
        }
    }
    MyMatrixXf connect_region_constrained_edge = MyMatrixXf::Zero(100, connect_region_constrained_edge_ori.cols());
    connect_region_constrained_edge.block(0, 0, connect_constrained_edge_ori_num, connect_region_constrained_edge_ori.cols())
        << connect_region_constrained_edge_ori.block(0, 0, connect_constrained_edge_ori_num, connect_region_constrained_edge_ori.cols());
    int old_connect_constrained_edge_num = connect_constrained_edge_ori_num;
    int connect_constrained_edge_num = old_connect_constrained_edge_num;
    MyMatrixXf coMat_ori = MyMatrixXf::Zero(100, 100);// Association matrix to ensure shortest path search does not self-connect
    // Process closed connected components (divide closed connected constraint edges into two segments, only extend one of them)
    MyMatrixXf new_connect_region_constrained_edge = MyMatrixXf::Zero(1, connect_region_constrained_edge.cols());
    for (int i = 0; i < old_connect_constrained_edge_num; i++) {
        int start_point = connect_region_constrained_edge(i, 1);
        r = connect_region_constrained_edge(i, 0);
        int end_point = connect_region_constrained_edge(i, r);
        if (start_point != end_point) {
            // If the connected component is not closed, find its connection with other components
            for (int j = 0; j < old_connect_constrained_edge_num; j++) {
                if (i == j)
                    continue;
                if (anyIsmember(connect_region_constrained_edge.block(j, 1, 1, connect_region_constrained_edge(j, 0)), start_point) ||
                    anyIsmember(connect_region_constrained_edge.block(j, 1, 1, connect_region_constrained_edge(j, 0)), end_point)) {
                    coMat_ori(i, j) = 1;
                    coMat_ori(j, i) = 1;
                }
            }
            continue;
        }
        int length = round((connect_region_constrained_edge(i, 0) - 1) / 2);
        new_connect_region_constrained_edge(0, 0) = connect_region_constrained_edge(i, 0) - length;
        new_connect_region_constrained_edge.block(0, 1, 1, new_connect_region_constrained_edge(0, 0)) = connect_region_constrained_edge.block(i, length + 1, 1, connect_region_constrained_edge(i, 0) - length);
        connect_region_constrained_edge(i, 0) = length + 1;
        connect_region_constrained_edge.block(i, 2 + length, 1, connect_region_constrained_edge.cols() - (length + 2)) = MyMatrixXf::Zero(1, connect_region_constrained_edge.cols() - (length + 2));// Note: the truncation point must be saved twice to ensure constrained edges are mutually connected
        connect_region_constrained_edge.row(connect_constrained_edge_num) = new_connect_region_constrained_edge;
        coMat_ori(i, connect_constrained_edge_num) = 1;
        coMat_ori(connect_constrained_edge_num++, i) = 1;
    }
    // Output
    
    MyMatrixXf coMat = coMat_ori.block(0, 0, connect_constrained_edge_num, connect_constrained_edge_num);
    // Update coMat using a graph traversal algorithm
    coMatClustering(coMat);
    for (int i = 0; i < connect_constrained_edge_num; i++) {
        coMat(i, i) = 1;// Constrained edge is associated with itself
    }
    // Shortest path search: determine the shortest path from each endpoint to other constrained edges or boundaries, then extend
    // Notes: 1) Extend the shortest path each time, then update; 2) Update the adjacency matrix each time

    // Build the adjacency matrix and undirected graph for region_point
    Graph G;
    for (int i = 0; i < region_edge_num; i++) {
        float weight = point2pointDist(point.row(region_edge(i, 0)), point.row(region_edge(i, 1)));
        add_edge(region_edge(i, 0), region_edge(i, 1), weight, G);
    }

    MyMatrixXf edge_shortest_path = MyMatrixXf::Zero(2 * old_connect_constrained_edge_num, region_point_num + 1);
    MyMatrixXf edge_p_size = MyMatrixXf::Zero(2 * old_connect_constrained_edge_num, 1);

    MyMatrixXf outline_point = getSubMat_Rows(region_point.col(0), getIndex(region_point.col(1), 0, 0));
    
    for (int i = 0; i < 2 * old_connect_constrained_edge_num; i++) {
        //Build new connection graph
        Graph G_temp = G;
        MyMatrixXf min_id;
        MyMatrixXf temp_union;
        int this_end_point;
        if (i % 2 == 0)
            this_end_point = connect_region_constrained_edge(i / 2, 1);
        else {
            r = connect_region_constrained_edge(i / 2, 0);
            this_end_point = connect_region_constrained_edge(i / 2, r);
        }
        if (anyIsmember(outline_point, this_end_point)) {
            edge_shortest_path(i, 0) = 0;
            edge_shortest_path(i, 1) = this_end_point;
            edge_p_size(i, 0) = 1;
            continue;
        }
        for (int j = 0; j < connect_constrained_edge_num; j++) {
            if (coMat(i / 2, j) == 1) {
                for (int k = 1; k < connect_region_constrained_edge(j, 0) + 1; k++) {
                    if (connect_region_constrained_edge(j, k) != this_end_point)
                        removeEdge(G_temp, connect_region_constrained_edge(j, k));
                }
            }
        }
        // Calculate the shortest path for each endpoint
        MyMatrixXf point_set = MyMatrixXf::Zero(region_point_num, 1);
        point_set.block(0, 0, outline_point.rows(), 1) << outline_point;
        int point_set_num = region_point.rows() - region_point.col(1).sum();
        for (int j = 0; j < connect_constrained_edge_num; j++) {
            temp_union = MyUnion(point_set.block(0, 0, point_set_num, 1), connect_region_constrained_edge.block(j, 1, 1, connect_region_constrained_edge(j, 0)));
            point_set.block(0, 0, temp_union.rows(), 1) = temp_union;
            point_set_num = temp_union.rows();
        }
        vector<vector<size_t>> path;
        MyMatrixXf dist = shortestPath(G_temp, path, this_end_point, point_set);
        for (int j = 0; j < dist.rows(); j++) {
            if (dist(j, 0) == 0)
                dist(j, 0) = std::numeric_limits<float>::max();
        }
        // Select the shortest path among all endpoints for extension and update accordingly
        getMin(dist, min_id);
        r = min_id(0, 0);
        edge_shortest_path(i, 0) = dist(r, 0);
        edge_shortest_path.block(i, 1, 1, path[r].size()) = vec2mat(path[r]).transpose();
        edge_shortest_path.block(i, 1, 1, path[r].size()) = matFlip(edge_shortest_path.block(i, 1, 1, path[r].size()));
        edge_p_size(i, 0) = path[r].size();
    }

    // Extend constrained edges in order of distance; stop if another constrained edge or extended edge is encountered
    vector<float> dist_vec = mat2vec(edge_shortest_path.col(0));
    vector<size_t> sort_path_index = sort_indexes(dist_vec, 1);
    MyMatrixXf extend_edge = MyMatrixXf::Zero(1, region_point_num);
    for (int i = 0; i < 2 * old_connect_constrained_edge_num; i++) {
        r = sort_path_index[i];
        if (edge_p_size(r, 0) != 0) {
            // If a path intersects with another constrained edge or path, truncate the path
            for (int j = 1; j < edge_p_size(r, 0) - 1; j++) {
                if (anyIsmember(ban_points.block(0, 0, ban_point_num, 1), edge_shortest_path(r, j + 1))) {
                    edge_p_size(r, 0) = j + 1;
                    break;
                } 
            }
            if (r % 2 == 1) {
                extend_edge.block(0, 0, 1, connect_region_constrained_edge(r / 2, 0) + edge_p_size(r, 0) - 1) << connect_region_constrained_edge.block(r / 2, 1, 1, connect_region_constrained_edge(r / 2, 0)), edge_shortest_path.block(r, 2, 1, edge_p_size(r, 0) - 1);
                connect_region_constrained_edge(r / 2, 0) += (edge_p_size(r, 0) - 1);
                connect_region_constrained_edge.block(r / 2, 1, 1, connect_region_constrained_edge(r / 2, 0)) = extend_edge.block(0, 0, 1, connect_region_constrained_edge(r / 2, 0));
            }
            else {
                extend_edge.block(0, 0, 1, connect_region_constrained_edge(r / 2, 0) + edge_p_size(r, 0) - 1) << matFlip(edge_shortest_path.block(r, 2, 1, edge_p_size(r, 0) - 1)), connect_region_constrained_edge.block(r / 2, 1, 1, connect_region_constrained_edge(r / 2, 0));
                connect_region_constrained_edge(r / 2, 0) += (edge_p_size(r, 0) - 1);
                connect_region_constrained_edge.block(r / 2, 1, 1, connect_region_constrained_edge(r / 2, 0)) = extend_edge.block(0, 0, 1, connect_region_constrained_edge(r / 2, 0));
            }
            ban_points.block(ban_point_num, 0, edge_p_size(r, 0) - 1, 1) = edge_shortest_path.block(r, 1, 1, edge_p_size(r, 0) - 1).transpose();
            ban_point_num += edge_p_size(r, 0) - 1;
        }
    }

    //Region segmentation
    MyMatrixXf extend_region_constrained_edge = MyMatrixXf::Zero(region_edge_num, 2);
    count = 0;
    MyMatrixXf this_constrained_edge = MyMatrixXf::Zero(1, 2);
    for (int i = 0; i < connect_constrained_edge_num; i++) {
        for (int j = 1; j < connect_region_constrained_edge(i, 0); j++) {
            this_constrained_edge << connect_region_constrained_edge(i, j), connect_region_constrained_edge(i, j + 1);
                extend_region_constrained_edge.row(count++) << this_constrained_edge;
        }
    }
    MyMatrixXf region_segment;
    
    bool exist_error = segmentRegionGrowth(face, region_face_idx, extend_region_constrained_edge.block(0, 0, count, 2), region_segment);
    if (exist_error)
        return region_face_idx;
    
    // Select the region with the shortest average node-to-constrained-edge distance after segmentation as the new supporting region
    int region_num = region_segment.rows();
    MyMatrixXf region_dist = MyMatrixXf::Zero(region_num, 1);
    for (int i = 0; i < region_num; i++) {
        region_dist(i, 0) = line2meshMedianDist(point, face, region_segment.block(i, 1, 1, region_segment(i, 0)).transpose(), line, interval);
    }
    MyMatrixXf min_id;
    myMin(region_dist, min_id, 0.005);
    MyMatrixXf new_region_face_idx;

    if (min_id.rows() == 1) {
        r = min_id(0, 0);
        new_region_face_idx = region_segment.block(min_id(0, 0), 1, 1, region_segment(r, 0)).transpose();
    }
    else {
        MyMatrixXf area = MyMatrixXf::Zero(min_id.rows(), 1);
        for (int i = 0; i < min_id.rows(); i++) {
            r = min_id(i, 0);
            area(i, 0) = getRegionArea(point, face, region_segment.block(r, 1, 1, region_segment(r, 0)).transpose());
        }
        MyMatrixXf min_area_id;
        getMin(area, min_area_id);
        r = min_area_id(0, 0); r = min_id(r, 0);
        new_region_face_idx = region_segment.block(r, 1, 1, region_segment(r, 0)).transpose();
    }
    //Support domain expansion to avoid incomplete support domain after segmentation
    MyMatrixXf temp = extendSegmentRegion(point, face, region_face_idx, base_region_face_idx, new_region_face_idx, line, region_constrained_edge, interval);
    return temp;
}


