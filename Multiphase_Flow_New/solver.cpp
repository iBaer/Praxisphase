#include "solver.h"

Solver::Solver(string methodname, Constants* constants, Computation *computation, Grid* grid)
{
	name = methodname;
	konstanten = constants;
	gs = computation;
	this->grid = grid;
    ordnung = 1;

    dimension = konstanten->dimension;
    CELLS = new int(dimension);

    mor = konstanten->mor;
    mol = konstanten->mol;
    mur = konstanten->mur; //2D
    mul = konstanten->mul; //2D
    timeou = konstanten->timeou;
    steps = 0;
    maxnt = konstanten->maxnt;
    teilerend = konstanten->teilerend;
    teiler = konstanten->teiler;
    variante = (int) konstanten->variante;


    CELLS[0] = konstanten->CELLSX;
    CELLS[1] = konstanten->CELLSY;
    dx = (mor-mol)/(double)CELLS[0];
    dy = (mur-mul)/(double)CELLS[1];

}
