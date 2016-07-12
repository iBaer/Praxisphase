/*
 * cluster_square.h
 *
 *  Created on: 27.06.2016
 *      Author: pascal
 */

#include <set>
#include "grid.h"

#ifndef CLUSTER_SQUARE_H_
#define CLUSTER_SQUARE_H_

class Cluster_Square {
public:
	Cluster_Square();
	virtual ~Cluster_Square();

	int pos_x_min;
	int pos_x_max;
	int pos_y_min;
	int pos_y_max;
	Grid* parent;
	Grid* grid_fine;
	std::set<int> adjacencies;

};

#endif /* CLUSTER_SQUARE_H_ */
