#include "initialNormal.h"
#include "myMath.h"
#include "dist.h"

MyMatrixXf point2lineVec(const MyMatrixXf& point, const MyMatrixXf& line);

float incrementGetMean(float mean, int n, float new_data);

MyMatrixXf ransacGetAngleDistribution(const MyMatrixXf& point, const MyMatrixXf& line) {
	int pts_num = point.rows();
	MyMatrixXf vector = MyMatrixXf::Zero(pts_num, 3);
	for (int i = 0; i < pts_num; i++) {
		vector.row(i) = point2lineVec(point.row(i), line);
	}
	MyMatrixXf sign = MyMatrixXf::Zero(pts_num, 1);
	for (int i = 0; i < pts_num; i++) {
		if (norm(vector.row(i)) != 0)
			sign(i, 0) = 1;
	}
	MyMatrixXf non_zero_vector = getSubMat_Rows(vector, getIndex(sign, 1, 0));
	MyMatrixXf center = ransacFindCenter(non_zero_vector, 7.5);
	int initial_num = 0;
	MyMatrixXf angle_d = MyMatrixXf::Zero(non_zero_vector.rows(), 1);
	for (int i = 0; i < non_zero_vector.rows(); i++) {
		float temp = getVecAngle(center, non_zero_vector.row(i)) * 180 / M_PI;
		if (temp <= 60) {
			angle_d(initial_num++, 0) = temp;
		}
	}
	return angle_d.block(0, 0, initial_num, 1);
}

MyMatrixXf ransacFindCenter(const MyMatrixXf& non_zero_data, float sigma) {
	MyMatrixXf center = MyMatrixXf::Zero(1, 3);
	int data_num = non_zero_data.rows();
	for (int i = 0; i < data_num; i++) {
		center += non_zero_data.row(i);
	}
	center = center / data_num;

	float pretotal = round(data_num * 0.6) - 0.1;
	int count = 0;
	for (int i = 0; i < data_num; i++) {
		float total = 0;
		for (int j = 0; j < data_num; j++) {
			if (getVecAngle(non_zero_data.row(i), non_zero_data.row(j)) * 180 / M_PI <= sigma)
				total += 1;
		}
		if (total > pretotal) {
			pretotal = total;
			center = non_zero_data.row(i);
			count += 1;
		}
		if (total == pretotal) {
			for (int m = 0; m < 3; m++) {
				center(0, m) = incrementGetMean(center(0, m), count, non_zero_data(i, m));
				count += 1;
			}
		}
	}
	return center;
}

MyMatrixXf point2lineVec(const MyMatrixXf& point, const MyMatrixXf& line) {
	MyMatrixXf vector = MyMatrixXf::Zero(1, 3);
	if (point2pointDist(point, line.row(0)) < 1e-6 || point2pointDist(point, line.row(1)) < 1e-6) {
		return vector;
	}
	MyMatrixXf para = lineEquation(line.row(0), line.row(1));

	float d = -(para(0, 0) * point(0, 0) + para(0, 1) * point(0, 1) + para(0, 2) * point(0, 2));

	float t = -(para(0, 0) * line(0, 0) + para(0, 1) * line(0, 1) + para(0, 2) * line(0, 2) + d) / (norm(para) * norm(para));
	MyMatrixXf O = line.row(0) + t * para;
	vector = point - O;
	if (norm(vector) != 0) {
		vector = vector / norm(vector);
	}
	return vector;
}

float incrementGetMean(float mean, int n, float new_data) {
	return (n * mean + new_data) / (n + 1);
}