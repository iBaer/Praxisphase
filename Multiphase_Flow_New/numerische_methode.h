#ifndef NUMERISCHE_METHODE_H
#define NUMERISCHE_METHODE_H
#include <iostream>

#include "computation.h"
#include "grid.h"
#include "solver.h"

/*!
 * @class numerische_methode
 * Die abstrakte Klasse numerische_methode gibt den Rahmen
 * für alle erbenden numerischen Methoden vor.
 */

class numerische_methode {
public:
	/**
	 * Konstruktor.
	 * @param dim Setzt Dimension für die Berechnung.
	 * @param ordn Ordnung der numerischen Methode.
	 * @param cells Anzahl der Zellen.
	 * @param method Name der Methode für unterscheidung bei Output.
	 */
	numerische_methode(Solver* solver, Constants* constants,
			Computation *computation, Grid *grid);
	~numerische_methode();
	/**
	 * Startet die Berechnung.
	 */
	void start_method();
	/**
	 * Name der Methode.
	 * Wird hauptsächlich dazu benötigt um die
	 * die Outputs der verschiedenen Methoden ausseinander
	 * halten zu können.
	 */
	std::string name;
	/**
	 * Array welches die Anzahl der Zellen in der entsprechenden Dimension zeigt.
	 */
	int* grid_size;
	/**
	 * Gamma Konstante
	 */
	double gamma;
	/**
	 * Zeit Output
	 */
	double time_output;
	/**
	 * Dimension in der gerechnet wird.
	 */
	int dimension;
	/**
	 * Bei dimensions==2, Integrationsschema splitting oder unsplitting
	 */
	int with_splitting;
	/**
	 * ordnung des Verfahrens.
	 */
	int order;
	/**
	 * Anzahl an Schritten die gemacht wurden.
	 */
	int steps;
	/**
	 * Mit Output für jeden Berechnungsschritt.
	 */
	int step_output;
	/**
	 * Gesetztes Maximum, damit die Methode nicht unendlich läuft
	 * (falls es einen fehler gibt oder andere umstände).
	 */
	int maxnt;
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
	 * teiler, Faktor für die ersten Delta t Schritte
	 */
	double divider;
	/**
	 * teilerend, Ende der Multiplikation der Zeitschritte mit teiler
	 */
	double divider_end;
	/**
	 * Delta y für die 2. Dimension
	 */
	double dy;
	/**
	 * Linke Grenze.
	 */
	int pos_x_min;
	/**
	 * Rechte Grenze.
	 */
	int pos_x_max;
	/**
	 * Obere Grenze.
	 */
	int pos_y_min;
	/**
	 * Untere Grenze.
	 */
	int pos_y_max;
	/**
	 * Variante der EOS.
	 */
	int clf_option;

	double ct;
	/**
	 * Konstanten Objekt welches für die berechnungen benötigt wird.
	 * @see Konstanten
	 */
	Constants* constants;
	/**
	 * Raster in den gerechnet wird.
	 * @see Raster
	 */
	Grid* grid_main;
	/**
	 * Gleichungssystem Objekt.
	 * @see Gleichungssystem
	 */
	Computation* computation;
	/**
	 * Solver Objekt.
	 * @see Solver
	 */
	Solver* solver;
	/**
	 * CFL Bedingung anwenden und neue Zeit berechnen.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 * @return neue Zeit.
	 */
	double cfl_condition(int n, double time);
	/**
	 * Abstrakte methode zur berechnung des Flusses der jeweiligen numerischen Methode.
	 * @return Matrix der Flüsse (1D)
	 */
	//virtual std::vector< std::vector< std::vector< std::vector <double> > > >calc_method_flux(int dir) =0;
	/**
	 * Aktualisiert alle zelle mithilfe des berechneten Flusses.
	 */
	void update(
			double* fi,
			int dir);
	/**
	 * Schreibt Ergebnisse in Dateien u,d,ur,p
	 */
	void write();
	void splitting();
	void unsplitting();
	/**
	 * bestimmt die Elemente der Jacobi-Matrix in 3 Varianten, abhängig von der Wahl der EOS
	 */
	void matrix_1d(double * values, int n, double uone, double utwo,
			double uthree, double p, double done, double dtwo, double ccl,
			double g, double ct, double cref, int variante);

	/**
	 * bestimmt die Elemente der Jacobi-Matrix in 3 Varianten, abhängig von der Wahl der EOS
	 */
	void matrix_2d(double * values_x, double * values_y, int n, double uone,
			double utwo, double uthree, double ufour, double ufive, double p,
			double done, double dtwo, double ccl, double g, double ct,
			double cref);

	//void cfl_1d_eigenvalues(int n, double &time);
	//void cfl_1d_approx(int n, double &time);
	//void cfl_2d_eigenvalues(int n, double &time);
	//void cfl_1d_eigenvalues(int n, double &time);
	//void cfl_1d_eigenvalues(int n, double &time);

};

#endif // NUMERISCHE_METHODE_H
