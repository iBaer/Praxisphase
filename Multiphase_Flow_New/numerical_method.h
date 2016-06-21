#ifndef NUMERICAL_METHOD_H
#define NUMERICAL_METHOD_H
#include <iostream>

#include "computation.h"
#include "grid.h"
#include "solver.h"

/*!
 * @class Numerical_Method
 * Die abstrakte Klasse numerische_methode gibt den Rahmen
 * für alle erbenden numerischen Methoden vor.
 */

class Numerical_Method {
public:
	/**
	 * Konstruktor.
	 * @param solver Pointer auf das Objekt des jeweiligen numerischen Lösers für die PDG.
	 * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
	 * @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird
	 * @param grid Pointer auf Raster-Objekt.
	 */
	Numerical_Method(Solver* solver, Constants* constants, Computation *computation, Grid *grid);

	/**
	 * Destruktor
	 */
	~Numerical_Method();

	/**
	 * Startet die numerische Berechnung der Mehrphasenströmung.
	 */
	void start_method();

	/**
	 * Name der Methode.
	 * Benötigt zur Trennung der Outfiles der jeweiligen Löser
	 */
	//TODO: Remove
	std::string solver_name;

	/**
	 * Array, dass entsprecht der Rasterdimension die totale Größe der jeweiligen Richtung des Rasters enthält.
	 */
	int* grid_size;

	/**
	 * Gamma Konstante
	 */
	double gamma;

	/**
	 * Grenze für dt
	 */
	double time_limit;

	/**
	 * Dimension des Rasters
	 */
	int dimension;

	/**
	 * Ab 2D, Integrationsschema splitting oder unsplitting
	 */
	int with_splitting;

	/**
	 * Ordnung des Verfahrens.
	 */
	int order;

	/**
	 * Anzahl an Schritten die gemacht wurden.
	 */
	int steps;

	/**
	 * Mit Output für jeden Berechnungsschritt.
	 */
	int output_per_step;

	/**
	 * Maximale Anzahl an Schritten der numerischen Methode
	 * (falls es einen Fehler gibt oder andere Umstände).
	 */
	int step_limit;

	/**
	 * Delta t.
	 */
	double dt;

	/**
	 * Delta t_fixed.
	 */
	//double dt_fixed; //TODO: implementieren?

	/**
	 * Delta x.
	 */
	double dx;

	/**
	 * Faktor für die ersten Delta t Schritte
	 */
	double divider;

	/**
	 * Ende der Multiplikation der Zeitschritte mit divider
	 */
	double divider_last;

	/**
	 * Delta y für die 2. Dimension
	 */
	//TODO: für alle Dimensionen verallgemeinern
	double dy;

	/**
	 * Linke Grenze des Rasters (Koordinate)
	 */
	//TODO: für alle Dimensionen verallgemeinern
	int pos_x_min;

	/**
	 * Rechte Grenze des Rasters
	 */
	int pos_x_max;

	/**
	 * Obere Grenze des Rasters
	 */
	int pos_y_min;

	/**
	 * Untere Grenze des Rasters
	 */
	int pos_y_max;

	/**
	 * Variante der EOS (Equation of States/Zustandsgleichung)
	 */
	int cfl_option;

	/**
	 * Die konstante Größe CT aus dem Fortranprogramm ist äquivalent zu K_g
	 */
	double ct;

	/**
	 * Pointer auf das Objekt, welches die Konstanten enthält.
	 * @see Constants
	 */
	Constants* constants;

	/**
	 * Pointer auf Raster-Objekt.
	 * @see Grid
	 */
	Grid* grid_main;

	/**
	 * Pointer auf das Objekt, das für die Berechnungen benötigt wird.
	 * @see Computation
	 */
	Computation* computation;

	/**
	 * Pointer auf das Objekt des jeweiligen numerischen Lösers für die PDG.
	 * @see Solver
	 */
	Solver* solver;

	Time_Step_Calculation* time_calculation;


	/**
	 * Schreibt Ergebnisse in Dateien für d,p und u,ur für die jeweilige Dimension
	 */
	void write();

};

#endif // NUMERICAL_METHOD_H
