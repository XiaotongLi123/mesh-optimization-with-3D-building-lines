#include "lineProcess.h"
#include "myMath.h"
#include "dist.h"

MyMatrixXf getNode(const MyMatrixXf& line, float interval) {
	float length = point2pointDist(line.block(0, 0, 1, 3), line.block(1, 0, 1, 3));
	int segment_num = round(length / interval);
	if (segment_num == 0) {
		segment_num = 1;
	}

	interval = point2pointDist(line.block(0, 0, 1, 3), line.block(1, 0, 1, 3)) / segment_num;
	int node_num = segment_num + 1;

	MyMatrixXf para = lineEquation(line.block(0, 0, 1, 3), line.block(1, 0, 1, 3));

	MyMatrixXf node = MyMatrixXf::Zero(node_num, 3);
	float t;
	for (int i = 0; i < node_num; i++) {
		if (i != node_num) {
			t = i * interval / sqrt(norm(para));
			node(i, 0) = line(0, 0) + para(0, 0) * t;
			node(i, 1) = line(0, 1) + para(0, 1) * t;
			node(i, 2) = line(0, 2) + para(0, 2) * t;
		}
		else {
			node.block(i, 0, 1, 3) = line.block(1, 0, 1, 3);
		}
	}
	return node;
}

void lineSort(MyMatrixXf& line) {
	int line_num = line.rows() / 2;
	vector<float> length(line_num);
	for (int i = 0; i < line_num; i++) {
		length[i] = point2pointDist(line.block(2 * i, 0, 1, 3), line.block(2 * i + 1, 0, 1, 3));
	}
	vector<size_t> sort = sort_indexes(length, 0);
	MyMatrixXf line_ori = line;
	for (int i = 0; i < line_num; i++) {
		line.block(2 * i, 0, 2, 3) = line_ori.block(2 * sort[i], 0, 2, 3);
	}
}
