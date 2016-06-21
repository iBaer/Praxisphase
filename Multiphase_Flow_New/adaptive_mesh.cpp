/*
 * adaptive_mesh.cpp
 *
 *  Created on: 20.06.2016
 *      Author: pascal
 */

#include "adaptive_mesh.h"

#include "grid.h"
#include "solver.h"
#include "constants.h"

Adaptive_Mesh::Adaptive_Mesh(Solver * solver, Grid *grid, Constants *constants, Time_Step_Calculation* time_calculation) {
	// TODO Auto-generated constructor stub
	this->solver = solver;
	this->grid_main = grid;
	this->constants = constants;
	this->grid = new Grid(grid_main->grid_size_total[0],grid_main->grid_size_total[1], constants);
	this->time_calculation = time_calculation;

}

Adaptive_Mesh::~Adaptive_Mesh() {
	// TODO Auto-generated destructor stub
	delete grid;
}

void Adaptive_Mesh::amr() {
	// 1. Gröbstes Netz berechnen
	double dt = time_calculation->dt;
	solver->calc_method_flux(dt);

	// 2. Welche Zellen verfeinern

	// 3. Feineres Gitter erzeugen

	// 4. Rekursiv feinere Gitter berechnen

	// 5. Gröberes Gitter mit feinerem Gitter verbessern

	// 6. Next



}
