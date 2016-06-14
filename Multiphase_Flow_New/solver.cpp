#include "solver.h"

Solver::Solver(string solver_name, Constants* constants, Computation *computation, Grid* grid)
{

	name = solver_name;
	this->constants = constants;
	this->computation = computation;
	this->grid = grid;
    dimension = constants->dimension;

    size_total = new int[dimension];
    size_m1 = new int[dimension];

    for (int i = 0; i < dimension; i++){
    	size_total[i] = grid->grid_size_total[i];
    	size_m1[i] = size_total[i] - 1;
    }

    //TODO: dz; oder "delta array"
    dx = (constants->pos_x_max - constants->pos_x_min)/(double)grid->grid_size[0];
    if (dimension ==2){
        dy = (constants->pos_y_max - constants->pos_y_min)/(double)grid->grid_size[1];
    }

    time_calculation = NULL;
}

Solver::~Solver(){
	delete[] size_total;
	delete[] size_m1;

}

