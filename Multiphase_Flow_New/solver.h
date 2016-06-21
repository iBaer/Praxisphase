#ifndef SOLVER_H_
#define SOLVER_H_

#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>

#include "computation.h"
#include "constants.h"
#include "time_step_calculation.h"

using namespace std;

/*!
 * @class Solver
 * Abstrakte Klasse für sämtliche Methoden zur Lösung von Differential Gleichungen
 */

class Solver {
public:
	/**
	 * Konstruktor.
	 * @param dim Setzt Dimension für die Berechnung.
	 * @param ordn Ordnung der numerischen Methode.
	 * @param cells Anzahl der Zellen.
	 * @param method Name der Methode für unterscheidung bei Output.
	 */
	Solver(string solver_name, Constants* constants, Computation *computation,
			Grid* grid);

	virtual ~Solver();



	/**
	 * Abstrakte Methode zur berechnung des Flusses der jeweiligen numerischen Methode.
	 * @return Matrix der Flüsse (1D)
	 */
	virtual void calc_method_flux(double dt) = 0;

	/**
	 * Konstanten Objekt welches für die berechnungen benötigt wird.
	 * @see Konstanten
	 */
	Constants *constants;
	/**
	 * Raster in den gerechnet wird.
	 * @see Raster
	 */
	Grid *grid;
	/**
	 * Gleichungssystem Objekt.
	 * @see Gleichungssystem
	 */
	Computation *computation;

	string name;

	/**
	 * Dimension in der gerechnet wird.
	 */
	int dimension;
	/**
	 * Delta x.
	 */
	double dx;
	/**
	 * Delta y für die 2. Dimension
	 */
	double dy;

	/**
	 * Größe des Rasters
	 */
	int *size_total;
	/**
	 * Größe des Rasters für die Flussberechnung
	 */
	int *size_m1;

	int split_method;

	Time_Step_Calculation* time_calculation;

};

#endif /* SOLVER_H_ */
