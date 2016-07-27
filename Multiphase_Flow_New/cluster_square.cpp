/*
 * cluster_square.cpp
 *
 *  Created on: 27.06.2016
 *      Author: pascal
 */

#include "cluster_square.h"

Cluster_Square::Cluster_Square() {
	pos_x_min = -1;
	pos_x_max = -1;
	pos_y_min = -1;
	pos_y_max = -1;
	parent = nullptr;
	grid_fine = nullptr;
}

Cluster_Square::~Cluster_Square() {
	if (grid_fine != nullptr)
		delete grid_fine;
}

