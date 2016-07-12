/*
 * adaptive_mesh.h
 *
 *  Created on: 20.06.2016
 *      Author: pascal
 */

#ifndef ADAPTIVE_MESH_H_
#define ADAPTIVE_MESH_H_

#include "numerical_method.h"
#include "grid.h"
#include "constants.h"
#include "cluster_square.h"
#include <set>

class Adaptive_Mesh {
public:
	Adaptive_Mesh(Solver* solver, Grid *grid, Constants *constants, Time_Step_Calculation* time_calculation);
	virtual ~Adaptive_Mesh();
	void amr();
	void grid_marker(int* &marked_cells, Grid* grid_one, Grid* grid_two, double tolerance);
	int binary_clustering(int* &marked_cells, int x_max, int y_max);
	void clusterize(int* &clustered_cells, int* grid_size, int nth_cluster, int x, int y);
	void square_clustering(vector<Cluster_Square> &clusters, int* &marked_cells, int x_max, int y_max, int cluster_amount);
	int cluster_adjacency_check(vector<Cluster_Square> &clusters, int* &marked_cells, int x_max, int y_max, int cluster_amount);
	void square_cluster_merge(vector<Cluster_Square> &clusters, int* &marked_cells, int x_max, int y_max, int &cluster_amount);
	void create_fine_grid(Cluster_Square& cluster, Grid* parent_grid);


	Solver* solver;
	Grid* grid_main;
	Grid* grid_twostep;
	Grid* grid_doublestep;
	Constants * constants;
	Time_Step_Calculation* time_calculation;

	int* marked_cells;

};

#endif /* ADAPTIVE_MESH_H_ */
