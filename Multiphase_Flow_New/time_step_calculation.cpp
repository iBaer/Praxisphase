/*
 * eigenvalues.cpp
 *
 *  Created on: 08.06.2016
 *      Author: pascal
 */

#include <lapackpp/gmd.h>
#include <lapackpp/lavd.h>
#include <lapackpp/laslv.h>
#include "time_step_calculation.h"
#include "cfl_1d.h"
#include "cfl_2d.h"
#include <iomanip>

using namespace std;

Time_Step_Calculation::Time_Step_Calculation(int neqs, Grid * grid) {
	// TODO Auto-generated constructor stub
	this->constants = &Constants::instance();
	this->grid = grid;

	this->n_eqns = neqs;
	gamma = constants->gamma;
	ct = constants->ct;
	dimension = constants->dimension;
	time_limit = constants->timeou;
	//steps = 0;
	//step_limit = constants->maxnt;
	divider_last = constants->teilerend;
	divider = constants->teiler;
	double pos_x_max = constants->pos_x_max;
	double pos_x_min = constants->pos_x_min;
	double pos_y_max = constants->pos_y_max; //2D
	double pos_y_min = constants->pos_y_min; //2D
	dx = (pos_x_max - pos_x_min) / (double) constants->grid_size_x;
	dy = (pos_y_max - pos_y_min) / (double) constants->grid_size_y;

	dt = 0.0;
	dt_old = 0.0;

	cfl_option = (int) constants->variante;

	n=0;
	time = 0;
	time_old = 0;
	with_splitting = 1;

	v_max_new = new double[2];
	v_max_old = new double[2];

	compare_step = 0;

}

Time_Step_Calculation::~Time_Step_Calculation() {
	delete[] v_max_new;
	delete[] v_max_old;
}

/**
 *****************************************************************************************
 * CFL Bedingung anwenden und neue Zeit berechnen.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 * @return neue Zeit.
 *****************************************************************************************/
double Time_Step_Calculation::cfl_condition(Grid * grid) {

	this->grid = grid;

	switch (dimension) {
	case (1): {
		//1D
		if ((int) constants->cfl_with_eig == 1)
		//1 = true, v_max wird über Eigenwerte bestimmt
		{
			cfl_1d(cfl_1d_eigenvalues());
		}

		else {
			cfl_1d_approx();
		}

		break;

	}
	case (2): {
		//2 Dimensionen

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, v_max wird über Eigenwerte bestimmt

		if ((int) constants->cfl_with_eig == 0) {
			cfl_2d_approx();
		}

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, v_max wird über Eigenwerte bestimmt
		else {
			//Eigenvalues* eig = new Eigenvalues(computation->neqs,grid_main);
			cfl_2d(cfl_2d_eigenvalues());
		}

		break;
	}
	}
	return dt;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D mit Eigenwerten.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Time_Step_Calculation::cfl_1d(double v_max) {
	double cfl = constants->cfl;


	dt = cfl * dx / v_max;

	if (n <= divider_last)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D mit Eigenwerten.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Time_Step_Calculation::cfl_2d(double* v_max) {
	double cfl = constants->cfl;

	if (compare_step == 0) {
		//time_old = time;
	}
	else{

	}


	if (with_splitting == 1) {
		dt = cfl / (v_max[0] / dx + v_max[1] / dy);
	} else {
		dt = cfl / max(v_max[0] / dx, v_max[1] / dy);
	}

	if (n <= divider_last)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;

	if(with_splitting>=2){
		if(compare_step==0){
			//time = time + dt;
		}
	}
	else{
		//time = time + dt;
	}

	if (with_splitting >= 1) { // 2

		if (compare_step == 0) {
			v_max_old[0] = v_max[0];
			v_max_old[1] = v_max[1];
			dt_old = dt;
			compare_step = 1;
		} else {
			v_max_new[0] = v_max[0];
			v_max_new[1] = v_max[1];
			compare_step = 0;
		}
	}
	cout << "Größte Eigenwerte: " << v_max[0] << " und " << v_max[1] << endl;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;

	delete[] v_max;
}

double Time_Step_Calculation::halve_dt(double dt) {
	cout << "=== HALVED DT ==="<<endl;
	dt = dt * 0.5	;

	return dt;
}

double Time_Step_Calculation::set_new_time(double dt) {
	time = time + dt;

	return time;
}

int Time_Step_Calculation::compare_dt() {

	double dt_comp = dt - dt_old;
	cout << "Delta ts" << endl;
	cout << "dt:  " << dt << " | dtold: " << dt_old << endl;
	cout << "dtc: " << dt_comp << endl;


	if (dt < dt_old) {
		return -1;
	}
	else {
		//dt <= dt_new;
		return 0;
	}
}

void Time_Step_Calculation::reset_step() {

		//time = time_old;
		dt = dt_old;
		v_max_new[0] = v_max_old[0];
		v_max_new[1] = v_max_old[1];

		compare_step = 1;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D mit Eigenwerten.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
double Time_Step_Calculation::cfl_1d_eigenvalues() {
	double cref = constants->cref;
	//double cfl = constants->cfl;
	double gamma_inv = 1.0 / gamma;

	double p = 0.0, dtwo = 0.0;

	//ehemals smax für "maximum wave speed"
	double v_max = 0.0;

	//int n_eqns = computation->neqs;

	double* u_eqns = new double[n_eqns];

	double values[9];
	LaVectorDouble real(3);
	LaVectorDouble img(3);
	LaGenMatDouble vr(3, 3);

	//Schritt 1: Maxima finden
	for (int x = 0; x < grid->grid_size_total[0]; x++) {

		// ACHTUNG, NUR VARIANTE 1 IST IN DEN GLEICHUNGEN IMPLEMENTIERT!
		switch (cfl_option) {
		case (1):

			p = ct * pow(grid->cellsgrid[x][0], gamma);

			grid->cellsgrid[x][1] = p;
			dtwo = pow((p / cref), gamma_inv);
			break;
			/*case(2):
			 dtwo = ccl/((1/raster->cellsgrid[x][0])-((1-ccl)/done));
			 p = cref * pow(dtwo,g);
			 raster->cellsgrid[x][1] = p;;
			 break;
			 case(3):
			 dtwo = ccl/((1/raster->cellsgrid[x][0])-((1-ccl)/done));
			 p = ct* pow(raster->cellsgrid[x][0],g);
			 raster->cellsgrid[x][1] = p;;
			 break;*/
		}

		u_eqns[0] = grid->cellsgrid[x][0];
		u_eqns[1] = grid->cellsgrid[x][2] * u_eqns[0];
		u_eqns[2] = grid->cellsgrid[x][3];

		//Schritt 2: Einsetzen in die Jacobi-Matrix
		matrix_1d(values, u_eqns, p, dtwo, cfl_option);

		//Schritt 3: Berechnen der Eigenwerte
		LaGenMatDouble A(values, 3, 3, true);
		LaEigSolve(A, real, img, vr);

		//Schritt 4: Höchsten Eigenwert suchen
		for (int n = 0; n < 3; n++) {
			if (v_max < fabs(real(n)))
				v_max = fabs(real(n));
		}
	}

	printf("using eig, computed in all cells: %10.4e\n", v_max);
/*
	dt = cfl * dx / v_max;

	if (n <= divider_last)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;
*/
	delete[] u_eqns;

	return v_max;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D über Annäherung.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Time_Step_Calculation::cfl_1d_approx() {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double ccl = constants->ccl;
	double rho_one = constants->rho_one;
	double gamma_inv = 1.0 / gamma;

	double rho = 0.0, p = 0.0, uxr = 0.0, ux = 0.0, rho_two = 0.0;

	//TODO: max_u, max_u_step
	double v_max = 0.0, v_max_step = 0.0;

	//dt über Näherung berechnen
	for (int i = 0; i < grid->grid_size_total[0]; i++) {

		rho = grid->cellsgrid[i][0];
		uxr = grid->cellsgrid[i][3];
		ux = grid->cellsgrid[i][2];

		switch (cfl_option) {
		case (1): {
			p = ct * pow(rho, gamma);
			grid->cellsgrid[i][1] = p;

			rho_two = pow((p / cref), gamma_inv);
			break;
		}
		case (2): {
			rho_two = ccl / ((1 / rho) - ((1 - ccl) / rho_one));
			p = cref * pow(rho_two, gamma);
			grid->cellsgrid[i][1] = p;

			break;
		}
		case (3): {
			rho_two = ccl / ((1 / rho) - ((1 - ccl) / rho_one));
			p = ct * pow(rho, gamma);
			grid->cellsgrid[i][1] = p;

			break;
		}
		}

		//ehemals maxs
		v_max_step = fabs(ux) + AM_1d(gamma, p, rho) + XS_1d(rho, rho_two, rho_one, ccl, uxr);

		if (v_max_step > v_max)
			v_max = v_max_step;
	}
	printf("using lambda exact %10.4e\n", v_max);

	dt = cfl * dx / v_max;

	if (n <= divider_last)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;

}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 2D mit Eigenwerten.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
double* Time_Step_Calculation::cfl_2d_eigenvalues() {
	double cref = constants->cref;
	//double cfl = constants->cfl;
	double gamma_inv = 1.0 / this->gamma;
	double rho = 0.0;
	double u_x = 0.0, ur_x = 0.0, u_y = 0.0, ur_y = 0.0;

	double p = 0.0, rho_two = 0.0;

	double* v_max = new double[2];
	v_max[0] = 0, v_max[1] = 0;

	//int n_eqns = computation->neqs;

	double* u_eqns = new double[n_eqns];

	//ueber Eigenwerte

	int n_eqns_squared = n_eqns * n_eqns;
	double values1[n_eqns_squared];
	double values2[n_eqns_squared];

	/**** Eventuell als Unteraufruf, um dt überprüfen zu können bei Splitting ****/

	//Schritt 1: Daten holen
	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {

			rho = grid->cellsgrid[x + y * grid->grid_size_total[0]][0];
			u_x = grid->cellsgrid[x + y * grid->grid_size_total[0]][2];
			u_y = grid->cellsgrid[x + y * grid->grid_size_total[0]][4];
			ur_x = grid->cellsgrid[x + y * grid->grid_size_total[0]][3];
			ur_y = grid->cellsgrid[x + y * grid->grid_size_total[0]][5];

			p = ct * pow(grid->cellsgrid[x + y * grid->grid_size_total[0]][0],gamma);
			grid->cellsgrid[x + y * grid->grid_size_total[0]][1] = p;
			rho_two = pow((p / cref), gamma_inv);

			u_x = u_x * rho;
			u_y = u_y * rho;

			//Schritt 2: Einsetzen
			u_eqns[0] = rho;
			u_eqns[1] = u_x;
			u_eqns[2] = u_y;
			u_eqns[3] = ur_x;
			u_eqns[4] = ur_y;

			matrix_2d(values1, values2, u_eqns, p, rho_two);

			//Schritt 3: Berechnen der Eigenwerte
			LaGenMatDouble A(values1, n_eqns, n_eqns, true);
			LaGenMatDouble B(values2, n_eqns, n_eqns, true);
			LaVectorDouble real(n_eqns);
			LaVectorDouble real_b(n_eqns);
			LaVectorDouble img(n_eqns);
			LaVectorDouble img_b(n_eqns);
			LaGenMatDouble vr(n_eqns, n_eqns);
			LaGenMatDouble vr_b(n_eqns, n_eqns);

			LaEigSolve(A, real, img, vr);
			LaEigSolve(B, real_b, img_b, vr_b);

			//Schritt 4: Höchsten Eigenwert suchen
			for (int n = 0; n < n_eqns; n++) {
				if (v_max[0] < fabs(real(n)))
					v_max[0] = fabs(real(n));
				if (v_max[1] < fabs(real_b(n)))
					v_max[1] = fabs(real_b(n));
			}
		}
	}

	delete[] u_eqns;

	return v_max;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 2D über Annäherung.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Time_Step_Calculation::cfl_2d_approx() {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double ccl = constants->ccl;
	double rho_one = constants->rho_one;
	double gamma_inv = 1.0 / gamma;

	double rho = 0.0, p = 0.0, uxr = 0.0, ux = 0.0, uy = 0.0, uyr = 0.0, rho_two = 0.0;

	//TODO: max_u, max_u_step
	double v_max = 0.0, v_max_step = 0.0;

	//über Näherung
	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {

			int index = x + y * grid->grid_size_total[0];
			rho = grid->cellsgrid[index][0];
			uxr = grid->cellsgrid[index][3];
			ux = grid->cellsgrid[index][2];
			uy = grid->cellsgrid[index][4];
			uyr = grid->cellsgrid[index][5];
			p = ct * pow(rho, gamma);
			grid->cellsgrid[x + y * grid->grid_size_total[0]][1] = p;
			rho_two = pow((p / cref), gamma_inv);

			v_max_step = fabs(ux) + fabs(uy) + AM_2d(gamma, p, rho) + XS_2d(rho, rho_two, rho_one, ccl, uxr, uxy);

			if (v_max_step > v_max) {
				v_max = v_max_step;
			}
		}
	}

	cout << "Größte geratene Eigenwerte: " << v_max << endl;

	dt = cfl * min(dx, dy) / v_max;

	if (n <= divider_last)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;


}



/**
 *****************************************************************************************
 *  Bestimmt die Elemente der Jacobi-Matrix für eine Dimension
 *****************************************************************************************/

//TODO: Parameter kürzen
//TODO: Matrix hängt von neqs², ist hier aber fix
void Time_Step_Calculation::matrix_1d(double * values, double * u, double p, double dtwo, int variante) {
	double cref = constants->cref;
	double ccl = constants->ccl;
	double rho_one = constants->rho_one;

	// Index für u[index] um 1 niedriger als in den Gleichungen
	switch (variante) {
	case (1):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[2] * u[2] + gamma * ct * pow(u[0], gamma - 1.0);
		values[4] = (2 * u[1]) / u[0];
		values[5] = u[0] * ccl * (1.0 - ccl) * 2.0 * u[2];
		values[6] = -u[1] / (u[0] * u[0]) * u[2] + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / rho_one;
		values[7] = u[2] / u[0];
		values[8] = (u[1] / u[0]) + (1.0 - 2.0 * ccl) * u[2];
		break;
	case (2):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[2] * u[2]
				+ cref * gamma * pow(rho_one, gamma - 1.0) * (ccl * rho_one * rho_one) / pow(rho_one + (ccl - 1.0) * u[0], 2.0);
		values[4] = (2 * u[1]) / u[0];
		values[5] = u[0] * ccl * (1.0 - ccl) * 2.0 * u[2];
		values[6] = -u[1] / (u[0] * u[0]) * u[2]
				+ cref * (gamma * pow(dtwo, gamma - 2.0) - (gamma * pow(dtwo, gamma - 1.0)) / rho_one) * (ccl * rho_one * rho_one) / pow(rho_one + (ccl - 1.0) * u[0], 2.0);
		values[7] = u[2] / u[0];
		values[8] = (u[1] / u[0]) + (1.0 - 2.0 * ccl) * u[2];
		break;
	case (3):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[2] * u[2] + gamma * ct * pow(u[0], gamma - 1.0);
		values[4] = (2.0 * u[1]) / u[0];
		values[5] = u[0] * ccl * (1.0 - ccl) * 2.0 * u[2];
		values[6] = -u[1] / (u[0] * u[0]) * u[2] + (1.0 / ccl) * ((1.0 / u[0]) - (1.0 / rho_one)) * ct * gamma * pow(u[0], gamma - 1.0);
		values[7] = u[2] / u[0];
		values[8] = (u[1] / u[0]) + (1.0 - 2.0 * ccl) * u[2];
	}
}

/**
 *****************************************************************************************
 *  Bestimmt die Elemente der Jacobi-Matrix für zwei Dimensionen
 *****************************************************************************************/
void Time_Step_Calculation::matrix_2d(double * values_x, double * values_y, double * u, double p, double dtwo) {

	double cref = constants->cref;
	double ccl = constants->ccl;
	double rho_one = constants->rho_one;

	values_x[0] = 0;
	values_x[1] = 1;
	values_x[2] = 0;
	values_x[3] = 0;
	values_x[4] = 0;
	values_x[5] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[3] * u[3] + gamma * ct * pow(u[0], gamma - 1.0);
	values_x[6] = (2 * u[1]) / u[0];
	values_x[7] = 0;
	values_x[8] = ccl * (1 - ccl) * 2 * u[0] * u[3];
	values_x[9] = 0;
	values_x[10] = -(u[1] * u[2]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[3] * u[4];
	values_x[11] = u[2] / u[0];
	values_x[12] = u[1] / u[0];
	values_x[13] = ccl * (1 - ccl) * u[0] * u[4];
	values_x[14] = ccl * (1 - ccl) * u[0] * u[3];
	values_x[15] = -(u[1] * u[3]) / (u[0] * u[0]) + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / (rho_one);
	values_x[16] = u[3] / u[0];
	values_x[17] = 0;
	values_x[18] = u[1] / u[0] + (1 - 2 * ccl) * u[3];
	values_x[19] = 0;
	values_x[20] = -(u[1] * u[4]) / (u[0] * u[0]) - u[2] * u[3] / (u[0] * u[0]);
	values_x[21] = u[4] / u[0];
	values_x[22] = u[3] / u[0];
	values_x[23] = u[2] / u[0] + (1 - 2 * ccl) * u[4];
	values_x[24] = u[1] / u[0] + (1 - 2 * ccl) * u[3];

	values_y[0] = 0;
	values_y[1] = 0;
	values_y[2] = 1;
	values_y[3] = 0;
	values_y[4] = 0;
	values_y[5] = -(u[1] * u[2]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[3] * u[4];
	values_y[6] = u[2] / u[0];
	values_y[7] = u[1] / u[0];
	values_y[8] = ccl * (1 - ccl) * u[0] * u[4];
	values_y[9] = ccl * (1 - ccl) * u[0] * u[3];
	values_y[10] = -(u[2] * u[2]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[4] * u[4] + gamma * ct * pow(u[0], gamma - 1.0);
	values_y[11] = 0;
	values_y[12] = (2 * u[2]) / u[0];
	values_y[13] = 0;
	values_y[14] = ccl * (1 - ccl) * 2 * u[0] * u[4];
	values_y[15] = -(u[1] * u[4]) / (u[0] * u[0]) - u[2] * u[3] / (u[0] * u[0]);
	values_y[16] = u[4] / u[0];
	values_y[17] = u[3] / u[0];
	values_y[18] = u[2] / u[0] + (1 - 2 * ccl) * u[4];
	values_y[19] = u[1] / u[0] + (1 - 2 * ccl) * u[3];
	values_y[20] = -(u[2] * u[4]) / (u[0] * u[0]) + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / (rho_one);
	values_y[21] = 0;
	values_y[22] = u[4] / u[0];
	values_y[23] = 0;
	values_y[24] = u[2] / u[0] + (1 - 2 * ccl) * u[4];
}
