#pragma once
#include "class.h"

// Get equally spaced nodes along the line segment
MyMatrixXf getNode(const MyMatrixXf& line, float interval);

// Sort line segments by length
void lineSort(MyMatrixXf& line);


