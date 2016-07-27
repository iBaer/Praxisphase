/*
 * cluster_builder.h
 *
 *  Created on: 25.07.2016
 *      Author: pascal
 */

#ifndef CLUSTER_BUILDER_H_
#define CLUSTER_BUILDER_H_

#include "cluster_square.h"
#include <set>

class Cluster_Builder {
public:
	Cluster_Builder();
	virtual ~Cluster_Builder();

	void build_clusters(Grid* grid, int* &marked_cells, std::vector<Cluster_Square>& clusters);
	int binary_clustering(int* &marked_cells, int x_max, int y_max);
	void clusterize(int* &clustered_cells, int* grid_size, int nth_cluster, int x, int y);
	void square_clustering(std::vector<Cluster_Square> &clusters, int* &marked_cells, Grid* grid, int cluster_amount);
	int cluster_adjacency_check(std::vector<Cluster_Square> &clusters, int* &marked_cells, int cluster_amount);
	void square_cluster_merge(std::vector<Cluster_Square> &clusters, int* &marked_cells, Grid* grid, int &cluster_amount);

	int* marked_cells;
	Grid* grid_main;

};

#endif /* CLUSTER_BUILDER_H_ */
