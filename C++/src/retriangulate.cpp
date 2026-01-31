#include "retriangulate.h"
#include "myMath.h"
#include "meshOper.h"
#include "dist.h"
#include "myMath.h"
#include "time.h"
#include "initialNormal.h"
#include "omp.h"


float getRegionLength(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx);

void postProcess(const MyMatrixXf& input_point, const MyMatrixXf& node_idx, const MyMatrixXf& edge_point_sort, MyMatrixXf& add_constrained_line, MyMatrixXf& map, MyMatrixXf& sign);

MyMatrixXf faceChange(MyMatrixXf old_face, MyMatrixXf map);

MyMatrixXf findConstrainedOutline(const MyMatrixXf& point, const MyMatrixXf& outline_point_sort, const MyMatrixXf& line);

MyMatrixXf lineDirectionChange(const MyMatrixXf& node);

float punishFun(float angle, int sign, float weight);

int getStoragePos(int i, int s, int n);

float getThreshold(vector<float>& distribution, float confidence, int data_num);

int partition(std::vector<float>& nums, int left, int right);

float quickSelect(std::vector<float>& nums, int left, int right, int k);

void update(const MyMatrixXf& point, const MyMatrixXf& node, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, const MyMatrixXf& line,
	const int left, const int right, MyMatrixXf& opt_face, float& opt_val, float weight_1, float weight_2, bool left_sign, bool right_sign);


bool regionRemeshing(const MyMatrixXf& point, MyMatrixXf& face, MyMatrixXf& sign, const MyMatrixXf& region_face_idx, const MyMatrixXf& node_idx, const MyMatrixXf& PFneighbor, MyMatrixXf& constrained_line, float w_1, float w_2,int line_sign) {
	int input_face_num = face.rows();
	int input_region_face_num = region_face_idx.rows();
	MyMatrixXf node = getSubMat_Rows(point, node_idx);
	MyMatrixXf new_face_ori = MyMatrixXf::Zero(input_face_num + 10 * input_region_face_num, 3);
	new_face_ori.block(0, 0, input_face_num, 3) = face;
	bool error = false;

	MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
	MyMatrixXf region_point, region_edge, edge_point_sort, edge_sort;
	getRegionPE(region_face, region_point, region_edge);
	getOutlineSort(region_face, edge_point_sort, edge_sort);
	int outline_num = edge_point_sort.rows();
	MyMatrixXf line = MyMatrixXf::Zero(2, 3);
	line << node.row(0),
		node.row(node.rows() - 1);
		
	MyMatrixXf constrained_outline = findConstrainedOutline(point, edge_point_sort, line);
	if (constrained_outline.rows() != 0) {
		MyMatrixXf new_constrained_line = MyMatrixXf::Zero(constrained_line.rows() + constrained_outline.rows(), 2);
		new_constrained_line << constrained_line,
			constrained_outline;
		constrained_line = new_constrained_line;
		return error;
	}

	MyMatrixXf normal_angle = regionFaceNormalAngle(point, face, region_face_idx, line_sign);
	float weight_1 = w_1 / getRegionLength(point, face, region_face_idx);
	float weight_2 = w_2 / (normal_angle.sum() + 1);

	MyMatrixXf enlarge_face_idx = getPointNeighFace(region_point.col(0), PFneighbor);
	MyMatrixXf neigh_face_normal = MyMatrixXf::Zero(outline_num, 3);
	MyMatrixXf is_constrained = MyMatrixXf::Zero(outline_num, 1);
	MyMatrixXf edge_neigh_face_id, edge_neigh_face_idx, exist_line_id;
	MyMatrixXf neigh_face_idx = MyMatrixXf::Zero(1, 1);
	MyMatrixXf neigh_face = MyMatrixXf::Zero(1, 3);
	for (int i = 0; i < outline_num; i++) {
		edge_neigh_face_id = findEdgeNeighFace(getSubMat_Rows(face, enlarge_face_idx), edge_sort.row(i));
		edge_neigh_face_idx = getSubMat_Rows(enlarge_face_idx, edge_neigh_face_id);
		if (edge_neigh_face_idx.rows() == 1)
			continue;
		neigh_face_idx = MyDifference(edge_neigh_face_idx, region_face_idx);
		if (neigh_face_idx.rows() == 0)
			continue;
		neigh_face = face.row(neigh_face_idx(0, 0));
		is_constrained(i, 0) = isExistLine(constrained_line, edge_sort.row(i), exist_line_id);
		neigh_face_normal.row(i) = getFacesNormal(point, neigh_face);
	}

	float opt_val = std::numeric_limits<float>::max();
	MyMatrixXf opt_face = DP_remeshing(getSubMat_Rows(point, edge_point_sort), node, neigh_face_normal, is_constrained, weight_1, weight_2, opt_val);

	if (opt_val == std::numeric_limits<float>::max()) {
		error = true;
		return error;
	}

	MyMatrixXf add_constrained_line, map;
	postProcess(point, node_idx, edge_point_sort, add_constrained_line, map, sign);
	MyMatrixXf new_constrained_line = MyMatrixXf::Zero(constrained_line.rows() + add_constrained_line.rows(), 2);
	if (constrained_line.rows() != 0 && add_constrained_line.rows() != 0) {
		new_constrained_line << constrained_line,
			add_constrained_line;
	}
	else {
		if (constrained_line.rows() == 0)
			new_constrained_line = add_constrained_line;
		if (add_constrained_line.rows() == 0)
			new_constrained_line = constrained_line;
	}
	constrained_line = new_constrained_line;
	MyMatrixXf opt_face_c = faceChange(opt_face, map);
	int reconstruct_face_num = opt_face_c.rows();

	valueSubMat(new_face_ori, region_face_idx, MyMatrixXf::Zero(1, 1), -1);
	int new_face_num = input_face_num - input_region_face_num + reconstruct_face_num;
	MyMatrixXf new_face = MyMatrixXf::Zero(new_face_num, 3);
	int count = 0;
	for (int i = 0; i < input_face_num; i++) {
		if (new_face_ori(i, 0) != -1) {
			new_face.row(count++) = new_face_ori.row(i);
		}
	}
	new_face.block(count, 0, reconstruct_face_num, 3) = opt_face_c;
	face = new_face;

	MyMatrixXf index = MyMatrixXf::Zero(reconstruct_face_num, 1);
	for (int i = 0; i < reconstruct_face_num; i++) {
		index(i, 0) = new_face_num - reconstruct_face_num + i;
	}

	valueSubMat(sign, getSubMat_Rows(region_point.col(0), getIndex(region_point.col(1), 1, 0)), MyMatrixXf::Zero(1, 1), 0);
	return error;
}

MyMatrixXf faceChange(MyMatrixXf old_face, MyMatrixXf map) {
	int face_num = old_face.rows();
	MyMatrixXf new_face = MyMatrixXf::Zero(face_num, 3);
	int r;
	for (int i = 0; i < face_num; i++) {
		for (int j = 0; j < 3; j++) {
			r = old_face(i, j);
			new_face(i, j) = map(r, 0);
		}
	}
	return new_face;
}


float getRegionLength(const MyMatrixXf& point, const MyMatrixXf& face, const MyMatrixXf& region_face_idx) {
	MyMatrixXf region_face = getSubMat_Rows(face, region_face_idx);
	MyMatrixXf length = MyMatrixXf::Zero(region_face_idx.rows(), 1);
	for (int i = 0; i < region_face_idx.rows(); i++) {
		length(i, 0) = point2pointDist(point.row(region_face(i, 0)), point.row(region_face(i, 1))) +
			point2pointDist(point.row(region_face(i, 1)), point.row(region_face(i, 2))) +
			point2pointDist(point.row(region_face(i, 2)), point.row(region_face(i, 0)));
	}
	return length.sum();
}

void postProcess(const MyMatrixXf& input_point, const MyMatrixXf& node_idx, const MyMatrixXf& edge_point_sort, MyMatrixXf& add_constrained_line, MyMatrixXf& map, MyMatrixXf& sign) {
	int r, c;
	MyMatrixXf exist_point_id_1, exist_point_id_2;
	int node_num = node_idx.rows();
	int input_point_num = node_idx(0, 0);
	MyMatrixXf node = getSubMat_Rows(input_point, node_idx);
	bool is_exist_point_1 = isExistPoint(getSubMat_Rows(input_point, edge_point_sort), node.row(0), exist_point_id_1);
	bool is_exist_point_2 = isExistPoint(getSubMat_Rows(input_point, edge_point_sort), node.row(node_num - 1), exist_point_id_2);

	add_constrained_line = MyMatrixXf::Zero(node_num - 1, 2);

	if (is_exist_point_1 == false) {
		if (is_exist_point_2 == false) {
			map = MyMatrixXf::Zero(edge_point_sort.rows() + node_num, 1);
			map.block(0, 0, edge_point_sort.rows(), 1) = edge_point_sort;
			for (int i = 0; i < node_num; i++) {
				map(edge_point_sort.rows() + i, 0) = i + input_point_num;
			}
			valueSubMat(sign, node_idx, MyMatrixXf::Zero(1, 1), 1);
			for (int i = 0; i < node_num - 1; i++) {
				add_constrained_line.row(i) << i + input_point_num, i + 1 + input_point_num;
			}
		}
		else {
			map = MyMatrixXf::Zero(edge_point_sort.rows() + node_num - 1, 1);
			map.block(0, 0, edge_point_sort.rows(), 1) = edge_point_sort;
			for (int i = 0; i < node_num - 1; i++) {
				map(edge_point_sort.rows() + i, 0) = i + input_point_num;
			}
			valueSubMat(sign, node_idx.block(0, 0, node_num - 1, 1), MyMatrixXf::Zero(1, 1), 1);
			for (int i = 0; i < node_num - 2; i++) {
				add_constrained_line.row(i) << i + input_point_num, i + 1 + input_point_num;
			}
			c = exist_point_id_2(0, 0);
			add_constrained_line.row(node_num - 2) << node_num - 2 + input_point_num, edge_point_sort(c, 0);
		}
	}
	else {
		if (is_exist_point_2 == false) {
			map = MyMatrixXf::Zero(edge_point_sort.rows() + node_num - 1, 1);
			map.block(0, 0, edge_point_sort.rows(), 1) = edge_point_sort;
			for (int i = 1; i < node_num; i++) {
				map(edge_point_sort.rows() + i - 1, 0) = i + input_point_num;
			}
			valueSubMat(sign, node_idx.block(1, 0, node_num - 1, 1), MyMatrixXf::Zero(1, 1), 1);
			r = exist_point_id_1(0, 0);
			add_constrained_line.row(0) << edge_point_sort(r, 0), 1 + input_point_num;
			for (int i = 1; i < node_num - 1; i++) {
				add_constrained_line.row(i) << i + input_point_num, i + 1 + input_point_num;
			}
		}
		else {
			map = MyMatrixXf::Zero(edge_point_sort.rows() + node_num - 2, 1);
			map.block(0, 0, edge_point_sort.rows(), 1) = edge_point_sort;
			for (int i = 1; i < node_num - 1; i++) {
				map(edge_point_sort.rows() + i - 1, 0) = i + input_point_num;
			}
			valueSubMat(sign, node_idx.block(1, 0, node_num - 2, 1), MyMatrixXf::Zero(1, 1), 1);
			r = exist_point_id_1(0, 0); c = exist_point_id_2(0, 0);
			if (node_num > 2) {
				add_constrained_line.row(0) << edge_point_sort(r, 0), 1 + input_point_num;
				for (int i = 1; i < node_num - 2; i++) {
					add_constrained_line.row(i) << i + input_point_num, i + 1 + input_point_num;
				}
				add_constrained_line.row(node_num - 2) << node_num - 2 + input_point_num, edge_point_sort(c, 0);
			}
			else {
				r = exist_point_id_1(0, 0); c = exist_point_id_2(0, 0);
				add_constrained_line.row(0) << edge_point_sort(r, 0), edge_point_sort(c, 0);
			}
		}
	}
}

MyMatrixXf findConstrainedOutline(const MyMatrixXf& point, const MyMatrixXf& outline_point_sort, const MyMatrixXf& line) {
	int outline_point_num = outline_point_sort.rows();
	int left = -1;
	int right = -1;
	for (int i = 0; i < outline_point_num; i++) {
		if (point2pointDist(point.row(outline_point_sort(i, 0)), line.row(0)) < 1e-3) {
			left = i;
		}
		if (point2pointDist(point.row(outline_point_sort(i, 0)), line.row(1)) < 1e-3) {
			right = i;
		}
	}
	MyMatrixXf overlapping = MyMatrixXf::Ones(2, 1);
	if (left != right && left != -1 && right != -1) {
		for (int i = 0; i < outline_point_num; i++) {
			if (i >= min(left, right) && i <= max(left, right)) {
				if (point2segmentDist(point.row(outline_point_sort(i, 0)), line) > 1e-3)
					overlapping(0, 0) = 0;
			}
			else {
				if (point2segmentDist(point.row(outline_point_sort(i, 0)), line) > 1e-3)
					overlapping(1, 0) = 0;
			}
		}
	}
	if (overlapping(0, 0) == 1 && overlapping(1, 0) == 0) {
		MyMatrixXf constrained_outline = MyMatrixXf::Zero(max(left, right) - min(left, right), 2);
		for (int i = min(left, right); i < max(left, right); i++) {
			constrained_outline.row(i - min(left, right)) << outline_point_sort(i, 0), outline_point_sort(i + 1, 0);
		}
		return constrained_outline;
	}
	if (overlapping(1, 0) == 1 && overlapping(0, 0) == 0) {
		MyMatrixXf constrained_outline = MyMatrixXf::Zero(outline_point_num - (max(left, right) - min(left, right)), 2);
		for (int i = 0; i < min(left, right); i++) {
			constrained_outline.row(i) << outline_point_sort(i, 0), outline_point_sort(i + 1, 0);
		}
		for (int i = max(left, right); i < outline_point_num - 1; i++) {
			constrained_outline.row(i - max(left, right) + min(left, right)) << outline_point_sort(i, 0), outline_point_sort(i + 1, 0);
		}
		constrained_outline.row(outline_point_num - 1 - (max(left, right) - min(left, right))) << outline_point_sort(outline_point_num, 0), outline_point_sort(0, 0);
		return constrained_outline;
	}
	return MyMatrixXf::Zero(0, 2);
}

MyMatrixXf DP_remeshing(const MyMatrixXf& point, const MyMatrixXf& node, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, const float weight_1, const float weight_2, float& opt_val) {
	int n = point.rows();
	int node_num = node.rows();
	opt_val = std::numeric_limits<float>::max();
	float punish_weight = 1000;
	MyMatrixXf line = MyMatrixXf::Zero(2, 3);
	line << node.row(0),
		node.row(node_num - 1);

	MyMatrixXf exist_point_id_1, exist_point_id_2;
	bool is_exist_point_1 = isExistPoint(point, node.row(0), exist_point_id_1);
	bool is_exist_point_2 = isExistPoint(point, node.row(node_num - 1), exist_point_id_2);

	int left, right;
	MyMatrixXf opt_face;
	if (is_exist_point_1 == true && is_exist_point_2 == true) {
		left = exist_point_id_1(0, 0);
		right = exist_point_id_2(0, 0);
		update(point, node.block(1, 0, node_num - 2, 3), neigh_face_normal, is_constrained, line, left, right, opt_face, opt_val, weight_1, weight_2, 1, 1);
		return opt_face;
	}

	pcl::KdTreeFLANN<pcl::PointXYZ> tree = createKdtree(point);
	int neigh = min(3, n);
	if (is_exist_point_1 == true && is_exist_point_2 == false) {
		left = exist_point_id_1(0, 0);
		vector<int> index = getKneighbor(tree, node.row(node_num - 1), neigh);
		for (int i = 0; i < neigh; i++) {
			right = index[i];
			if (right == left)
				continue;
			update(point, node.block(1, 0, node_num - 1, 3), neigh_face_normal, is_constrained, line, left, right, opt_face, opt_val, weight_1, weight_2, 1, 0);
			if (opt_val < punish_weight)
				break;
		}
		return opt_face;
	}
	if (is_exist_point_2 == true && is_exist_point_1 == false) {
		right = exist_point_id_2(0, 0);
		vector<int> index = getKneighbor(tree, node.row(0), neigh);
		for (int i = 0; i < neigh; i++) {
			left = index[i];
			if (right == left)
				continue;
			update(point, node.block(0, 0, node_num - 1, 3), neigh_face_normal, is_constrained, line, left, right, opt_face, opt_val, weight_1, weight_2, 0, 1);
			if (opt_val < punish_weight)
				break;
		}
		return opt_face;
	}

	vector<int> index_1 = getKneighbor(tree, node.row(0), neigh);
	vector<int> index_2 = getKneighbor(tree, node.row(node_num - 1), neigh);
	for (int i = 0; i < neigh; i++) {
		for (int j = 0; j < neigh; j++) {
			left = index_1[i];
			right = index_2[j];
			if (right == left)
				continue;
			update(point, node, neigh_face_normal, is_constrained, line, left, right, opt_face, opt_val, weight_1, weight_2, 0, 0);
			if (opt_val < punish_weight) {
				return opt_face;
			}
		}
	}

	return opt_face;
}

void update(const MyMatrixXf& point, const MyMatrixXf& node, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, const MyMatrixXf& line,
	const int left, const int right, MyMatrixXf& opt_face, float& opt_val, float weight_1, float weight_2, bool left_sign, bool right_sign) {
	int pts_num = point.rows();
	int node_num = node.rows();
	int n = pts_num;

	MyMatrixXf neigh_face_normal_1, neigh_face_normal_2, is_constrained_1, is_constrained_2, point_1, point_2, map_1, map_2, distribution_1, distribution_2;
	if (left < right) {
		point_1 = MyMatrixXf::Zero(right - left + 1 + node_num, 3);
		point_1 << point.block(left, 0, right - left + 1, 3),
			lineDirectionChange(node);
		neigh_face_normal_1 = MyMatrixXf::Zero(point_1.rows(), 3);
		neigh_face_normal_1.block(0, 0, right - left, 3) = neigh_face_normal.block(left, 0, right - left, 3);
		is_constrained_1 = MyMatrixXf::Zero(point_1.rows(), 1);
		is_constrained_1.block(0, 0, right - left, 1) = is_constrained.block(left, 0, right - left, 1);
		point_2 = MyMatrixXf::Zero(pts_num - right + left + 1 + node_num, 3);
		point_2 << point.block(right, 0, pts_num - right, 3),
			point.block(0, 0, left + 1, 3),
			node;
		neigh_face_normal_2 = MyMatrixXf::Zero(point_2.rows(), 3);
		neigh_face_normal_2.block(0, 0, pts_num - right + left, 3) << neigh_face_normal.block(right, 0, pts_num - right, 3),
			neigh_face_normal.block(0, 0, left, 3);
		is_constrained_2 = MyMatrixXf::Zero(point_2.rows(), 1);
		is_constrained_2.block(0, 0, pts_num - right + left, 1) << is_constrained.block(right, 0, pts_num - right, 1),
			is_constrained.block(0, 0, left, 1);
		map_1 = MyMatrixXf::Zero(point_1.rows(), 1);
		map_2 = MyMatrixXf::Zero(point_2.rows(), 1);
		for (int i = 0; i < right - left + 1; i++) {
			map_1(i, 0) = i + left;
		}
		for (int i = 0; i < node_num; i++) {
			map_1(i + right - left + 1, 0) = n + node_num - 1 - i;
		}
		for (int i = 0; i < pts_num - right; i++) {
			map_2(i, 0) = i + right;
		}
		for (int i = 0; i < left + 1; i++) {
			map_2(i + pts_num - right, 0) = i;
		}
		for (int i = 0; i < node_num; i++) {
			map_2(i + pts_num - right + left + 1, 0) = n + i;
		}
		distribution_1 = ransacGetAngleDistribution(point_1, lineDirectionChange(line));
		distribution_2 = ransacGetAngleDistribution(point_2, line);
	}
	else {
		point_1 = MyMatrixXf::Zero(left - right + 1 + node_num, 3);
		point_1 << point.block(right, 0, left - right + 1, 3),
			node;
		neigh_face_normal_1 = MyMatrixXf::Zero(point_1.rows(), 3);
		neigh_face_normal_1.block(0, 0, left - right, 3) = neigh_face_normal.block(right, 0, left - right, 3);
		is_constrained_1 = MyMatrixXf::Zero(point_1.rows(), 1);
		is_constrained_1.block(0, 0, left - right, 1) = is_constrained.block(right, 0, left - right, 1);
		point_2 = MyMatrixXf::Zero(pts_num - left + right + 1 + node_num, 3);
		point_2 << point.block(left, 0, pts_num - left, 3),
			point.block(0, 0, right + 1, 3),
			lineDirectionChange(node);
		neigh_face_normal_2 = MyMatrixXf::Zero(point_2.rows(), 3);
		neigh_face_normal_2.block(0, 0, pts_num - left + right, 3) << neigh_face_normal.block(left, 0, pts_num - left, 3),
			neigh_face_normal.block(0, 0, right, 3);
		is_constrained_2 = MyMatrixXf::Zero(point_2.rows(), 1);
		is_constrained_2.block(0, 0, pts_num - left + right, 1) << is_constrained.block(left, 0, pts_num - left, 1),
			is_constrained.block(0, 0, right, 1);
		map_1 = MyMatrixXf::Zero(point_1.rows(), 1);
		map_2 = MyMatrixXf::Zero(point_2.rows(), 1);
		for (int i = 0; i < left - right + 1; i++) {
			map_1(i, 0) = i + right;
		}
		for (int i = 0; i < node_num; i++) {
			map_1(i + left - right + 1, 0) = n + i;
		}
		for (int i = 0; i < pts_num - left; i++) {
			map_2(i, 0) = i + left;
		}
		for (int i = 0; i < right + 1; i++) {
			map_2(i + pts_num - left, 0) = i;
		}
		for (int i = 0; i < node_num; i++) {
			map_2(i + pts_num - left + right + 1, 0) = n + node_num - 1 - i;
		}
		distribution_1 = ransacGetAngleDistribution(point_1, line);
		distribution_2 = ransacGetAngleDistribution(point_2, lineDirectionChange(line));
	}

	MyMatrixXf new_neigh_face_normal_1 = neigh_face_normal_1;
	MyMatrixXf new_is_constrained_1 = is_constrained_1;
	MyMatrixXf new_neigh_face_normal_2 = neigh_face_normal_2;
	MyMatrixXf new_is_constrained_2 = is_constrained_2;
	MyMatrixXf all_points = MyMatrixXf::Zero(pts_num + node_num, 3);
	all_points << point,
		node;
	float opt_val_1 = std::numeric_limits<float>::max();
	float opt_val_2 = std::numeric_limits<float>::max();

	MyMatrixXf opt_face_1 = DP_polygonRemeshing(point_1, neigh_face_normal_1, is_constrained_1, weight_1, weight_2, opt_val_1, distribution_1);
	MyMatrixXf opt_face_1_c = faceChange(opt_face_1, map_1);

	MyMatrixXf this_edge = MyMatrixXf::Zero(1, 2);
	MyMatrixXf edge_neigh_face_id;
	if (left < right) {
		for (int i = pts_num - right + left; i < point_2.rows(); i++) {
			if (i != point_2.rows() - 1) {
				this_edge << map_2(i, 0), map_2(i + 1, 0);
				edge_neigh_face_id = findEdgeNeighFace(opt_face_1_c, this_edge);
				if (i != pts_num - right + left)
					new_is_constrained_2(i, 0) = 1;
				else {
					if (left_sign == true)
						new_is_constrained_2(i, 0) = 1;
				}
			}
			else {
				this_edge << map_2(i, 0), map_2(0, 0);
				edge_neigh_face_id = findEdgeNeighFace(opt_face_1_c, this_edge);
				if (right_sign == true)
					new_is_constrained_2(i, 0) = 1;
			}
			if (edge_neigh_face_id.rows() == 0)
				continue;
			new_neigh_face_normal_2.row(i) = getFacesNormal(all_points, opt_face_1_c.row(edge_neigh_face_id(0, 0)));
		}
	}
	else {
		for (int i = pts_num - left + right; i < point_2.rows(); i++) {
			if (i != point_2.rows() - 1) {
				this_edge << map_2(i, 0), map_2(i + 1, 0);
				edge_neigh_face_id = findEdgeNeighFace(opt_face_1_c, this_edge);
				if (i != pts_num - left + right)
					new_is_constrained_2(i, 0) = 1;
				else {
					if (right_sign == true)
						new_is_constrained_2(i, 0) = 1;
				}
			}
			else {
				this_edge << map_2(i, 0), map_2(0, 0);
				edge_neigh_face_id = findEdgeNeighFace(opt_face_1_c, this_edge);
				if (left_sign == true)
					new_is_constrained_2(i, 0) = 1;
			}
			if (edge_neigh_face_id.rows() == 0)
				continue;
			new_neigh_face_normal_2.row(i) = getFacesNormal(all_points, opt_face_1_c.row(edge_neigh_face_id(0, 0)));
		}
	}
	MyMatrixXf opt_face_2 = DP_polygonRemeshing(point_2, new_neigh_face_normal_2, new_is_constrained_2, weight_1, weight_2, opt_val_2, distribution_2);
	MyMatrixXf opt_face_2_c = faceChange(opt_face_2, map_2);

	float this_opt_val = opt_val_1 + opt_val_2;
	MyMatrixXf opt_face_c = MyMatrixXf::Zero(opt_face_1_c.rows() + opt_face_2_c.rows(), 3);
	opt_face_c << opt_face_1_c,
		opt_face_2_c;
	MyMatrixXf idx = MyMatrixXf::Zero(opt_face_c.rows(), 1);
	for (int i = 0; i < opt_face_c.rows(); i++) {
		idx(i, 0) = i;
	}
	MyMatrixXf normal_angle_1 = regionFaceNormalAngle(all_points, opt_face_c, idx);
	if (normal_angle_1.maxCoeff() > 75) {
		opt_face_2 = DP_polygonRemeshing(point_2, neigh_face_normal_2, is_constrained_2, weight_1, weight_2, opt_val_2, distribution_2);
		opt_face_2_c = faceChange(opt_face_2, map_2);
		if (left < right) {
			for (int i = right - left; i < point_1.rows(); i++) {
				if (i != point_1.rows() - 1) {
					this_edge << map_1(i, 0), map_1(i + 1, 0);
					edge_neigh_face_id = findEdgeNeighFace(opt_face_2_c, this_edge);
					if (i != right - left)
						new_is_constrained_1(i, 0) = 1;
					else {
						if (right_sign == true)
							new_is_constrained_1(i, 0) = 1;
					}
				}
				else {
					this_edge << map_1(i, 0), map_1(0, 0);
					edge_neigh_face_id = findEdgeNeighFace(opt_face_2_c, this_edge);
					if (left_sign == true)
						new_is_constrained_1(i, 0) = 1;
				}
				if (edge_neigh_face_id.rows() == 0)
					continue;
				new_neigh_face_normal_1.row(i) = getFacesNormal(all_points, opt_face_2_c.row(edge_neigh_face_id(0, 0)));
			}
		}
		else {
			for (int i = left - right; i < point_1.rows(); i++) {
				if (i != point_1.rows() - 1) {
					this_edge << map_1(i, 0), map_1(i + 1, 0);
					edge_neigh_face_id = findEdgeNeighFace(opt_face_2_c, this_edge);
					if (i != left - right)
						new_is_constrained_1(i, 0) = 1;
					else {
						if (left_sign == true)
							new_is_constrained_1(i, 0) = 1;
					}
				}
				else {
					this_edge << map_1(i, 0), map_1(0, 0);
					edge_neigh_face_id = findEdgeNeighFace(opt_face_2_c, this_edge);
					if (right_sign == true)
						new_is_constrained_1(i, 0) = 1;
				}
				if (edge_neigh_face_id.rows() == 0)
					continue;
				new_neigh_face_normal_1.row(i) = getFacesNormal(all_points, opt_face_2_c.row(edge_neigh_face_id(0, 0)));
			}
		}
		opt_face_1 = DP_polygonRemeshing(point_1, new_neigh_face_normal_1, new_is_constrained_1, weight_1, weight_2, opt_val_1, distribution_1);
		opt_face_1_c = faceChange(opt_face_1, map_1);

		float new_this_opt_val = opt_val_1 + opt_val_2;
		if (new_this_opt_val < this_opt_val) {
			opt_face_c = MyMatrixXf::Zero(opt_face_1_c.rows() + opt_face_2_c.rows(), 3);
			opt_face_c << opt_face_1_c,
				opt_face_2_c;
			this_opt_val = new_this_opt_val;
			idx = MyMatrixXf::Zero(opt_face_c.rows(), 1);
			for (int i = 0; i < opt_face_c.rows(); i++) {
				idx(i, 0) = i;
			}
		}
	}
	
	if (this_opt_val < opt_val) {
		opt_val = this_opt_val;
		opt_face = opt_face_c;
		
		return;
	}
}


MyMatrixXf lineDirectionChange(const MyMatrixXf& node) {
	int node_num = node.rows();
	MyMatrixXf node_opposite = MyMatrixXf::Zero(node_num, 3);
	for (int i = 0; i < node_num; i++) {
		node_opposite.row(i) = node.row(node_num - 1 - i);
	}
	return node_opposite;
}

void getIJK(int dim_1, int dim_2, int dim_3, int pos, int& i, int& j, int& k) {
	i = floor(pos / (dim_2 * dim_3));
	j = floor(pos % (dim_2 * dim_3) / dim_3);
	k = (pos % (dim_2 * dim_3)) % dim_3;
}


MyMatrixXf DP_polygonRemeshing(const MyMatrixXf& point, const MyMatrixXf& neigh_face_normal, const MyMatrixXf& is_constrained, float weight_1, float weight_2, float& opt_val, MyMatrixXf& distribution) {
	int pts_num = point.rows();
	int n = pts_num;
	float confidence = 0.05;
	vector<faceSet> face((n - 1) * (n - 2) / 2);
	clock_t start = clock();

	weight_2 = weight_2 / weight_1;
	weight_1 = 1;

	MyMatrixXf dist = MyMatrixXf::Zero(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			dist(i, j) = point2pointDist(point.row(i), point.row(j));
			dist(j, i) = dist(i, j);
		}
	}
	MyMatrixXf normal = MyMatrixXf::Zero(n * n * n, 3);
	MyMatrixXf is_tri = MyMatrixXf::Zero(n * n * n, 1);
#pragma omp parallel for
	for (int pos = 0; pos < n * n * n; pos++) {
		MyMatrixXf tri = MyMatrixXf::Zero(3, 3);
		MyMatrixXf index = MyMatrixXf::Zero(1, 3);
		MyMatrixXf angle = MyMatrixXf::Zero(3, 1);
		index << 0, 1, 2;
		int i, j, k;
		getIJK(n, n, n, pos, i, j, k);
		if (i < j && j < k) {
			tri << point.row(i),
				point.row(j),
				point.row(k);
			angle(0, 0) = getLineAngle(point.row(j) - point.row(i), point.row(i) - point.row(k));
			angle(1, 0) = getLineAngle(point.row(k) - point.row(j), point.row(j) - point.row(i));
			angle(2, 0) = getLineAngle(point.row(i) - point.row(k), point.row(k) - point.row(j));
			if (angle.minCoeff() > 0.1 * M_PI / 180)
				is_tri(i * n * n + j * n + k, 0) = 1;
			if (is_tri(i * n * n + j * n + k, 0) == 1) {
				normal.row(i * n * n + j * n + k) = getFacesNormal(tri, index);
				if (norm(normal.row(i * n * n + j * n + k)) == 0) {
					is_tri(i * n * n + j * n + k, 0) = 0;
				}
			}
		}
	}

	MyMatrixXf opt_func = MyMatrixXf::Zero(n * n, n - 2);
	MyMatrixXf neigh_normal = MyMatrixXf::Zero(n * n, 3 * (n - 2));
	MyMatrixXf opt_idx = MyMatrixXf::Zero(n * n, 1);
	MyMatrixXf angle_record = MyMatrixXf::Zero(n * n, 2 * (n - 2));
	MyMatrixXf count_mat = MyMatrixXf::Zero(n * n, 1);
	MyMatrixXf opt_length = MyMatrixXf::Zero(n * n, n - 2);
	MyMatrixXf opt_Dangle = MyMatrixXf::Zero(n * n, n - 2);

	float neigh_angle_1 = 0;
	float neigh_angle_2 = 0;
	for (int i = 0; i < n - 2; i++) {
		int idx = getStoragePos(i, 3, n);
		face[idx].face_set = MyMatrixXf::Zero(1, 3);
		face[idx].face_set.block(0, 0, 1, 3) << i, i + 1, i + 2;

		if (is_tri(i * n * n + (i + 1) * n + (i + 2), 0) == 0) {
			opt_func(i * n + 2, 0) = std::numeric_limits<float>::max();
			continue;
		}

		neigh_angle_1 = 0;
		neigh_angle_2 = 0;
		if (!(neigh_face_normal(i, 0) == 0 && neigh_face_normal(i, 1) == 0 && neigh_face_normal(i, 2) == 0)) {
			neigh_angle_1 = getVecAngle(normal.row(i * n * n + (i + 1) * n + (i + 2)), neigh_face_normal.row(i)) * 180 / M_PI;
			neigh_angle_1 = punishFun(neigh_angle_1, is_constrained(i, 0), weight_2);
		}
		if (!(neigh_face_normal(i + 1, 0) == 0 && neigh_face_normal(i + 1, 1) == 0 && neigh_face_normal(i + 1, 2) == 0)) {
			neigh_angle_2 = getVecAngle(normal.row(i * n * n + (i + 1) * n + (i + 2)), neigh_face_normal.row(i + 1)) * 180 / M_PI;
			neigh_angle_2 = punishFun(neigh_angle_2, is_constrained(i + 1, 0), weight_2);
		}
		opt_length(i * n + 2, 0) = dist(i, i + 1) + dist(i + 1, i + 2) + dist(i, i + 2);
		opt_func(i * n + 2, 0) = weight_1 * opt_length(i * n + 2, 0) + neigh_angle_1 + neigh_angle_2;

		neigh_normal.block(i * n + 2, 0, 1, 3) = normal.row(i * n * n + (i + 1) * n + i + 2);
	}
	start = clock();
	for (int s = 4; s < n + 1; s++) {
		for (int i = 0; i < n - s + 1; i++) {
			int idx = getStoragePos(i, s, n);
			face[idx].face_set = MyMatrixXf::Zero(s - 2, 3 * (s - 2));
		}
	}

	int count = distribution.rows();
	int s_count = count;
	vector<float> new_distribution(9999999);
	vector<float> initial_distribution = mat2vec(distribution);
	std::copy(initial_distribution.begin(), initial_distribution.end(), new_distribution.begin());

	MyMatrixXf perimeter = MyMatrixXf::Zero(n * (n - 2), 1);
	MyMatrixXf Dangle = MyMatrixXf::Zero(n * (n - 2), 1);
	int sta_count = 0;
	for (int s = 4; s < n + 1; s++) {
		perimeter.block(0, 0, (n - s + 1) * (s - 2), 1) = MyMatrixXf::Zero((n - s + 1) * (s - 2), 1);
		Dangle.block(0, 0, (n - s + 1) * (s - 2), 1) = MyMatrixXf::Zero((n - s + 1) * (s - 2), 1);
		sta_count = 0;
		float threshold = getThreshold(new_distribution, confidence, s_count);
#pragma omp parallel for
		for (int pos = 0; pos < (n - 3) * (n - 1); pos++) {
			MyMatrixXf opt_id;
			int i = floor(pos / (n - 1));
			int k = pos % (n - 1);
			if (i >= n - s + 1 || k == 0 || k >= s - 1)
				continue;

			if (is_tri(i * n * n + (i + k) * n + (i + s - 1), 0) == 0) {
				opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
				count_mat(i * n + s - 1, 0) += 1;
				if (count_mat(i * n + s - 1, 0) == s - 2) {
					getMin(opt_func.block(i * n + s - 1, 0, 1, s - 2), opt_id);
					opt_idx(i * n + s - 1, 0) = opt_id(0, 1);
				}
				continue;
			}
			int c = opt_idx(i * n + k, 0);
			if (opt_func(i * n + k, c) == std::numeric_limits<float>::max()) {
				opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
				count_mat(i * n + s - 1, 0) += 1;
				if (count_mat(i * n + s - 1, 0) == s - 2) {
					getMin(opt_func.block(i * n + s - 1, 0, 1, s - 2), opt_id);
					opt_idx(i * n + s - 1, 0) = opt_id(0, 1);
				}
				continue;
			}

			c = opt_idx((i + k) * n + s - k - 1, 0);
			if (opt_func((i + k) * n + s - k - 1, c) == std::numeric_limits<float>::max()) {
				opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
				count_mat(i * n + s - 1, 0) += 1;
				if (count_mat(i * n + s - 1, 0) == s - 2) {
					getMin(opt_func.block(i * n + s - 1, 0, 1, s - 2), opt_id);
					opt_idx(i * n + s - 1, 0) = opt_id(0, 1);
				}
				continue;
			}

			int idx, idx_left, idx_right;
			if (k != 1 && k != s - 2) {
				int r, c;
				float angle_left = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block(i * n + k, 3 * opt_idx(i * n + k, 0), 1, 3)) * 180 / M_PI;
				float angle_right = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block((i + k) * n + s - k - 1, 3 * opt_idx((i + k) * n + s - k - 1, 0), 1, 3)) * 180 / M_PI;
				if (angle_left >= threshold && angle_right >= threshold) {
					MyMatrixXf temp_1 = MyMatrixXf::Zero(k + 1 - 2, 1);
					MyMatrixXf temp_2 = MyMatrixXf::Zero(1, s - (k + 1) - 1);
					MyMatrixXf left = MyMatrixXf::Zero(1, 3);
					MyMatrixXf right = MyMatrixXf::Zero(1, 3);
					for (int p = 1; p < (k + 1) - 1; p++) {
						left << neigh_normal.block(i * n + k, 3 * (p - 1), 1, 3);
						if (left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) {
							temp_1(p - 1, 0) = opt_func(i * n + k, p - 1) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2);
						}
						else {
							temp_1(p - 1, 0) = std::numeric_limits<float>::max();
						}
					}
					float opt_val_left = getMin(temp_1, opt_id);
					r = opt_id(0, 0);

					for (int q = 1; q < s - (k + 1); q++) {
						right << neigh_normal.block((i + k) * n + s - k - 1, 3 * (q - 1), 1, 3);
						if (right(0, 0) != 0 || right(0, 1) != 0 || right(0, 2) != 0) {
							temp_2(0, q - 1) = opt_func((i + k) * n + s - (k + 1), q - 1) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2);
						}
						else {
							temp_2(0, q - 1) = std::numeric_limits<float>::max();
						}
					}
					float opt_val_right = getMin(temp_2, opt_id);
					c = opt_id(0, 1);
					opt_func(i * n + s - 1, k - 1) = opt_val_left + opt_val_right + weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1));

				}
				else {
					if (angle_left < threshold && angle_right >= threshold) {
						MyMatrixXf temp = MyMatrixXf::Zero(1, s - (k + 1) - 1);
						MyMatrixXf right = MyMatrixXf::Zero(1, 3);
						for (int q = 1; q < s - (k + 1); q++) {
							right << neigh_normal.block((i + k) * n + s - k - 1, 3 * (q - 1), 1, 3);
							if (right(0, 0) != 0 || right(0, 1) != 0 || right(0, 2) != 0) {
								temp(0, q - 1) = opt_func((i + k) * n + s - (k + 1), q - 1) +
									weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
									punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2);
							}
							else {
								temp(0, q - 1) = std::numeric_limits<float>::max();
							}
						}
						opt_func(i * n + s - 1, k - 1) = getMin(temp, opt_id);
						r = opt_idx(i * n + k, 0);
						c = opt_id(0, 1);
						MyMatrixXf left = neigh_normal.block(i * n + k, 3 * r, 1, 3);
						if (left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) {
							opt_func(i * n + s - 1, k - 1) += (opt_func(i * n + k, r) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2));
						}
						else {
							opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						}
					}
					if (angle_left >= threshold && angle_right < threshold) {
						MyMatrixXf temp = MyMatrixXf::Zero(k + 1 - 2, 1);
						MyMatrixXf left = MyMatrixXf::Zero(1, 3);
						for (int p = 1; p < (k + 1) - 1; p++) {
							left << neigh_normal.block(i * n + k, 3 * (p - 1), 1, 3);
							if (left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) {
								temp(p - 1, 0) = opt_func(i * n + k, p - 1) +
									weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
									punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2);
							}
							else {
								temp(p - 1, 0) = std::numeric_limits<float>::max();
							}
						}
						opt_func(i * n + s - 1, k - 1) = getMin(temp, opt_id);
						r = opt_id(0, 0);
						c = opt_idx((i + k) * n + s - k - 1, 0);
						MyMatrixXf right = neigh_normal.block((i + k) * n + s - k - 1, 3 * c, 1, 3);
						if (right(0, 0) != 0 || right(0, 1) != 0 || right(0, 2) != 0) {
							opt_func(i * n + s - 1, k - 1) += (opt_func((i + k) * n + s - (k + 1), c) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2));
						}
						else {
							opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						}
					}
					if (angle_left < threshold && angle_right < threshold) {
						r = opt_idx(i * n + k, 0);
						c = opt_idx((i + k) * n + s - k - 1, 0);
						MyMatrixXf left = neigh_normal.block(i * n + k, 3 * r, 1, 3);
						MyMatrixXf right = neigh_normal.block((i + k) * n + s - k - 1, 3 * c, 1, 3);
						if ((left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) && (right(0, 0) != 0 || right(0, 1) != 0 || right(0, 2) != 0)) {
							opt_func(i * n + s - 1, k - 1) = opt_func(i * n + k, r) + opt_func((i + k) * n + s - (k + 1), c) +
								weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2);
						}
						else {
							opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						}
					}
				}
				if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
					opt_length(i * n + s - 1, k - 1) = opt_length(i * n + k, r) + opt_length((i + k) * n + s - (k + 1), c) + dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1);
					opt_Dangle(i * n + s - 1, k - 1) = opt_Dangle(i * n + k, r) + opt_Dangle((i + k) * n + s - (k + 1), c) +
						getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block(i * n + k, 3 * r, 1, 3)) * 180 / M_PI +
						getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block((i + k) * n + s - k - 1, 3 * c, 1, 3)) * 180 / M_PI;
					perimeter(sta_count, 0) = opt_length(i * n + s - 1, k - 1);
					Dangle(sta_count++, 0) = opt_Dangle(i * n + s - 1, k - 1);
				}
				else {
					opt_length(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
					opt_Dangle(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
				}

				idx = getStoragePos(i, s, n);
				idx_left = getStoragePos(i, k + 1, n);
				idx_right = getStoragePos(i + k, s - k, n);
				if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
					face[idx].face_set.block(k - 1, 0, 1, 3 * (k + 1 - 2)) = face[idx_left].face_set.block(r, 0, 1, 3 * (k + 1 - 2));
					face[idx].face_set.block(k - 1, 3 * (k + 1 - 2), 1, 3 * (s - k - 2)) = face[idx_right].face_set.block(c, 0, 1, 3 * (s - k - 2));
					face[idx].face_set.block(k - 1, 3 * (s - 3), 1, 3) << i, i + k, i + s - 1;
				}
				if (opt_func(i * n + s - 1, k - 1) < 500) {
					angle_record(i * n + s - 1, 2 * (k - 1)) = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), normal.row(i * n * n + (i + r + 1) * n + (i + k))) * 180 / M_PI;
					angle_record(i * n + s - 1, 2 * (k - 1) + 1) = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), normal.row((i + k) * n * n + (i + k + c + 1) * n + (i + s - 1))) * 180 / M_PI;
				}
				for (int neigh_count = 0; neigh_count < 2; neigh_count++) {
					if (angle_record(i * n + k, 2 * r + neigh_count) != 0) {
						new_distribution[count++] = angle_record(i * n + k, 2 * r + neigh_count);
					}
					if (angle_record((i + k) * n + s - k - 1, 2 * c + neigh_count) != 0) {
						new_distribution[count++] = angle_record((i + k) * n + s - k - 1, 2 * c + neigh_count);
					}
				}
			}
			else {
				if (k == 1) {
					int c;
					float angle_right = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block((i + k) * n + s - k - 1, 3 * opt_idx((i + k) * n + s - k - 1, 0), 1, 3)) * 180 / M_PI;
					if (angle_right >= threshold) {
						MyMatrixXf temp = MyMatrixXf::Zero(1, s - (k + 1) - 1);
						MyMatrixXf right = MyMatrixXf::Zero(1, 3);
						for (int q = 1; q < s - (k + 1); q++) {
							right << neigh_normal.block((i + k) * n + s - k - 1, 3 * (q - 1), 1, 3);
							if (right(0, 0) != 0 || right(0, 1) != 0 || right(0, 2) != 0) {
								temp(0, q - 1) = opt_func((i + k) * n + s - (k + 1), q - 1) +
									weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
									punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2);
							}
							else {
								temp(0, q - 1) = std::numeric_limits<float>::max();
								continue;
							}
							if (neigh_face_normal(i, 0) != 0 || neigh_face_normal(i, 1) != 0 || neigh_face_normal(i, 2) != 0)
								temp(0, q - 1) += punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_face_normal.row(i)) * 180 / M_PI, is_constrained(i, 0), weight_2);
						}

						opt_func(i * n + s - 1, k - 1) = getMin(temp, opt_id);
						c = opt_id(0, 1);
					}
					else {
						if (neigh_face_normal(i, 0) != 0 || neigh_face_normal(i, 1) != 0 || neigh_face_normal(i, 2) != 0)
							opt_func(i * n + s - 1, k - 1) = punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_face_normal.row(i)) * 180 / M_PI, is_constrained(i, 0), weight_2);
						c = opt_idx((i + k) * n + s - k - 1, 0);
						MyMatrixXf right = neigh_normal.block((i + k) * n + s - k - 1, 3 * c, 1, 3);
						if (right(0, 0) != 0 || right(0, 0) != 0 || right(0, 0) != 0) {
							opt_func(i * n + s - 1, k - 1) += (opt_func((i + k) * n + s - (k + 1), c) +
								weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), right) * 180 / M_PI, -1, weight_2));
						}
						else {
							opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						}
					}
					if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
						opt_length(i * n + s - 1, k - 1) = opt_length((i + k) * n + s - (k + 1), c) + dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1);
						opt_Dangle(i * n + s - 1, k - 1) = opt_Dangle((i + k) * n + s - (k + 1), c) +
							getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block((i + k) * n + s - k - 1, 3 * c, 1, 3)) * 180 / M_PI;
						perimeter(sta_count, 0) = opt_length(i * n + s - 1, k - 1);
						Dangle(sta_count++, 0) = opt_Dangle(i * n + s - 1, k - 1);
					}
					else {
						opt_length(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						opt_Dangle(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
					}

					idx = getStoragePos(i, s, n);
					idx_right = getStoragePos(i + k, s - k, n);
					if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
						face[idx].face_set.block(k - 1, 0, 1, 3 * (s - k - 2)) = face[idx_right].face_set.block(c, 0, 1, 3 * (s - k - 2));
						face[idx].face_set.block(k - 1, 3 * (s - k - 2), 1, 3) << i, i + k, i + s - 1;
					}

					if (opt_func(i * n + s - 1, k - 1) < 500) {
						angle_record(i * n + s - 1, 2 * (k - 1) + 1) = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), normal.row((i + k) * n * n + (i + k + c + 1) * n + (i + s - 1))) * 180 / M_PI;
					}
					for (int neigh_count = 0; neigh_count < 2; neigh_count++) {
						if (angle_record((i + k) * n + s - k - 1, 2 * c + neigh_count) != 0) {
							new_distribution[count++] = angle_record((i + k) * n + s - k - 1, 2 * c + neigh_count);
						}
					}
				}
				if (k == s - 2) {
					int r;
					float angle_left = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block(i * n + k, 3 * opt_idx(i * n + k, 0), 1, 3)) * 180 / M_PI;
					if (angle_left >= threshold) {
						MyMatrixXf temp = MyMatrixXf::Zero(k + 1 - 2, 1);
						MyMatrixXf left = MyMatrixXf::Zero(1, 3);
						for (int p = 1; p < (k + 1) - 1; p++) {
							left << neigh_normal.block(i * n + k, 3 * (p - 1), 1, 3);
							if (left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) {
								temp(p - 1, 0) = opt_func(i * n + k, p - 1) +
									weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
									punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2);
							}
							else {
								temp(p - 1, 0) = std::numeric_limits<float>::max();
								continue;
							}
							if (neigh_face_normal(i + s - 2, 0) != 0 || neigh_face_normal(i + s - 2, 1) != 0 || neigh_face_normal(i + s - 2, 2) != 0)
								temp(p - 1, 0) += punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_face_normal.row(i + s - 2)) * 180 / M_PI, is_constrained(i + s - 2, 0), weight_2);
						}
						opt_func(i * n + s - 1, k - 1) = getMin(temp, opt_id);
						r = opt_id(0, 0);
					}
					else {
						if (neigh_face_normal(i + s - 2, 0) != 0 || neigh_face_normal(i + s - 2, 1) != 0 || neigh_face_normal(i + s - 2, 2) != 0)
							opt_func(i * n + s - 1, k - 1) = punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_face_normal.row(i + s - 2)) * 180 / M_PI, is_constrained(i + s - 2, 0), weight_2);
						r = opt_idx(i * n + k, 0);
						MyMatrixXf left = neigh_normal.block(i * n + k, 3 * r, 1, 3);
						if (left(0, 0) != 0 || left(0, 1) != 0 || left(0, 2) != 0) {
							opt_func(i * n + s - 1, k - 1) += (opt_func(i * n + k, r) +
								weight_1 * (dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1)) +
								punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), left) * 180 / M_PI, -1, weight_2));
						}
						else {
							opt_func(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						}
					}
					if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
						opt_length(i * n + s - 1, k - 1) = opt_length(i * n + k, r) + dist(i, i + k) + dist(i + k, i + s - 1) + dist(i, i + s - 1);
						opt_Dangle(i * n + s - 1, k - 1) = opt_Dangle(i * n + k, r) +
							getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_normal.block(i * n + k, 3 * r, 1, 3)) * 180 / M_PI;
						perimeter(sta_count, 0) = opt_length(i * n + s - 1, k - 1);
						Dangle(sta_count++, 0) = opt_Dangle(i * n + s - 1, k - 1);
					}
					else {
						opt_length(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
						opt_Dangle(i * n + s - 1, k - 1) = std::numeric_limits<float>::max();
					}

					idx = getStoragePos(i, s, n);
					idx_left = getStoragePos(i, k + 1, n);
					if (opt_func(i * n + s - 1, k - 1) != std::numeric_limits<float>::max()) {
						face[idx].face_set.block(k - 1, 0, 1, 3 * (k + 1 - 2)) = face[idx_left].face_set.block(r, 0, 1, 3 * (k + 1 - 2));
						face[idx].face_set.block(k - 1, 3 * (k + 1 - 2), 1, 3) << i, i + k, i + s - 1;
					}
					if (opt_func(i * n + s - 1, k - 1) < 500) {
						angle_record(i * n + s - 1, 2 * (k - 1)) = getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), normal.row(i * n * n + (i + r + 1) * n + (i + k))) * 180 / M_PI;
					}
					for (int neigh_count = 0; neigh_count < 2; neigh_count++) {
						if (angle_record(i * n + k, 2 * r + neigh_count) != 0) {
							new_distribution[count++] = angle_record(i * n + k, 2 * r + neigh_count);
						}
					}
				}
			}
			neigh_normal.block(i * n + s - 1, 3 * (k - 1), 1, 3) = normal.row(i * n * n + (i + k) * n + (i + s - 1));
			if (s == n) {
				if (neigh_face_normal(n - 1, 0) != 0 || neigh_face_normal(n - 1, 1) != 0 || neigh_face_normal(n - 1, 2) != 0)
					opt_func(i * n + s - 1, k - 1) += punishFun(getVecAngle(normal.row(i * n * n + (i + k) * n + (i + s - 1)), neigh_face_normal.row(n - 1)) * 180 / M_PI, is_constrained(n - 1, 0), weight_2);
			}

			count_mat(i * n + s - 1, 0) += 1;
			if (count_mat(i * n + s - 1, 0) == s - 2) {
				getMin(opt_func.block(i * n + s - 1, 0, 1, s - 2), opt_id);
				opt_idx(i * n + s - 1, 0) = opt_id(0, 1);
			}
		}
		weight_2 = 1 / ((Dangle.block(0, 0, sta_count, 1).sum() + sta_count) * (s - 2) / (perimeter.block(0, 0, sta_count, 1).sum() * (s - 3)));

		s_count = count;

	}
	if (dp_count < 4) {
		dp_count++;
		myTime += clock() - start;
	}
	MyMatrixXf opt_id;
	opt_val = getMin(opt_func.block(n - 1, 0, 1, n - 2), opt_id);
	int idx = getStoragePos(0, n, n);
	MyMatrixXf opt_face = MyMatrixXf::Zero(n - 2, 3);
	for (int i = 0; i < n - 2; i++) {
		opt_face.row(i) = face[idx].face_set.block(opt_id(0, 1), 3 * i, 1, 3);
	}
	return opt_face;
}

float punishFun(float angle, int sign, float weight) {
	float lambda = 1000 * M_PI / 180;
	if (sign == 0) {
		if (angle > 105)
			return lambda * angle;
		else
			return 0;
	}
	if (sign == 1) {
		if (angle > 135)
			return lambda * angle;
		else
			return 0;
	}
	else {
		if (angle > 135)
			return lambda * angle;
		else {
			if (angle > 105)
				return 3 * weight * angle;
			else
				return weight * angle;
		}
	}
}

int getStoragePos(int i, int s, int n) {
	return (s - 3) * ((n - 2) + (n - (s - 2))) / 2 + i;
}

float getThreshold(vector<float>& vec, float confidence, int data_num) {
	if (data_num == 0)
		return 0.0;
	sort(vec.begin(), vec.begin() + data_num);
	int max_num = ceil((1 - confidence) * data_num);

	if (max_num != 0)
		max_num -= 1;
	return vec[max_num];
}

int partition(std::vector<float>& nums, int left, int right) {
	float pivot = nums[right];
	int i = left - 1;

	for (int j = left; j < right; ++j) {
		if (nums[j] <= pivot) {
			++i;
			std::swap(nums[i], nums[j]);
		}
	}
	std::swap(nums[i + 1], nums[right]);
	return i + 1;
}

float quickSelect(std::vector<float>& nums, int left, int right, int k) {
	if (left == right) {
		return nums[left];
	}

	int pivotIndex = partition(nums, left, right);

	if (pivotIndex == k) {
		return nums[k];
	}
	else if (pivotIndex < k) {
		return quickSelect(nums, pivotIndex + 1, right, k);
	}
	else {
		return quickSelect(nums, left, pivotIndex - 1, k);
	}
}