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
#include "cluster_builder.h"
#include <set>

class Adaptive_Mesh {
public:
	Adaptive_Mesh(Solver* solver, Grid *grid, Constants *constants, Time_Step_Calculation* time_calculation);
	virtual ~Adaptive_Mesh();
	void amr();
	int grid_marker(int* &marked_cells, Grid* grid_one, Grid* grid_two, double tolerance);
	void create_fine_grid(Cluster_Square& cluster, Grid* parent_grid);
	void grid_to_file(std::string filename, Grid* grid);
	void fine_to_coarse(Cluster_Square& cluster, double dt);
	void coarse_correction(Cluster_Square& cluster, double dt);
	void calculate_grids(Grid* refinable_grid, int grid_level, double& dt);

	Solver* solver;
	Grid* grid_main;
	Grid* grid_twostep;
	Grid* grid_doublestep;
	Constants * constants;
	Time_Step_Calculation* time_calculation;
	Cluster_Builder* cluster_builder;

	int* marked_cells;
	//TODO: dx, dt sollte nur an einer Stelle im Programm bestimmt werden! -> Grid!?
	double dx;
	double dy;
};

#endif /* ADAPTIVE_MESH_H_ */
