/*
 * time_step_calculation.h
 *
 *  Created on: 08.06.2016
 *      Author: pascal
 */

#ifndef TIME_STEP_CALCULATION_H_
#define TIME_STEP_CALCULATION_H_

#include "constants.h"
#include "grid.h"

class Time_Step_Calculation {
public:
	Time_Step_Calculation(int neqs, Grid * grid);
	virtual ~Time_Step_Calculation();
	Constants* constants;
	Grid* grid_main;
	int n_eqns;
	double gamma;
	double ct;
	int n;
	double time;
	double time_old;

	int with_splitting;
	double* v_max_new;
	double* v_max_old;
	int compare_step;

	/**
	 * Grenze für dt
	 */
	double time_limit;

	/**
	 * Dimension des Rasters
	 */
	int dimension;

	/**
	 * Delta t.
	 */
	double dt;
	double dt_old;

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
	 * Variante der EOS (Equation of States/Zustandsgleichung)
	 */
	int cfl_option;

	/**
	 * CFL Bedingung anwenden und neue Zeit berechnen.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 * @return neue Zeit.
	 */
	double cfl_condition();

	double halve_dt();

	int compare_dt();


	/**
	 * Unterfunktion der Methode cfl_condition().
	 * Berechnung der CFL Bedingung in 1D mit Eigenwerten.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 */
	double cfl_1d_eigenvalues();

	/**
	 * Unterfunktion der Methode cfl_condition().
	 * Berechnung der CFL Bedingung in 1D über Annäherung.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 */
	void cfl_1d_approx();

	/**
	 * Berechnet die CFL Bedingung in 1D
	 * @param v_max größer Eigenwert oder Approximation.
	 */
	void cfl_1d(double v_max);

	/**
	 * Berechnet die CFL Bedingung in 2D
	 * @param v_max größer Eigenwerte oder Approximation.
	 */
	void cfl_2d(double* v_max);

	/**
	 * Unterfunktion der Methode cfl_condition().
	 * Berechnung der CFL Bedingung in 2D mit Eigenwerten.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 */
	double* cfl_2d_eigenvalues();

	/**
	 * Unterfunktion der Methode cfl_condition().
	 * Berechnung der CFL Bedingung in 2D über Annäherung.
	 * @param n aktueller Zeitschritt.
	 * @param time aktuelle Zeit.
	 */
	void cfl_2d_approx();

	/**
	 * Bestimmt die Elemente der Jacobi-Matrix für eine Dimension
	 */
	void matrix_1d(double * values, double * u, double p, double dtwo, int variante);

	/**
	 * Bestimmt die Elemente der Jacobi-Matrix für zwei Dimensionen
	 */
	void matrix_2d(double * values_x, double * values_y,  double * u, double p, double dtwo);

};

#endif /* TIME_STEP_CALCULATION_H_ */
