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

class Adaptive_Mesh {
public:
	Adaptive_Mesh(Solver* solver, Grid *grid, Constants *constants, Time_Step_Calculation* time_calculation);
	virtual ~Adaptive_Mesh();
	void amr();

	Solver* solver;
	Grid* grid_main;
	Grid* grid;
	Constants * constants;
	Time_Step_Calculation* time_calculation;

};

#endif /* ADAPTIVE_MESH_H_ */
