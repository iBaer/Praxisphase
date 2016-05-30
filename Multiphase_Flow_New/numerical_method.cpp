#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <lapackpp/gmd.h>
#include <lapackpp/lavd.h>
#include <lapackpp/laslv.h>

#include "cfl_1d.h"
#include "cfl_2d.h"
#include "constants.h"
#include "numerical_method.h"

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor.
 * @param solver Pointer auf das Objekt des jeweiligen numerischen Lösers für die PDG.
 * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
 * @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird
 * @param grid Pointer auf Raster-Objekt.
 *******************************************1**********************************************/
Numerical_Method::Numerical_Method(Solver* solver, Constants* constants, Computation *computation, Grid *grid) {

	this->constants = constants;
	this->computation = computation;
	grid_main = grid;
	solver_name = solver->name;
	this->solver = solver;
	order = 1;

	dimension = constants->dimension;
	grid_size = new int[dimension];

	pos_x_max = constants->pos_x_max;
	pos_x_min = constants->pos_x_min;
	pos_y_max = constants->pos_y_max; //2D
	pos_y_min = constants->pos_y_min; //2D
	time_limit = constants->timeou;
	steps = 0;
	step_limit = constants->maxnt;
	divider_end = constants->teilerend;
	divider = constants->teiler;
	cfl_option = (int) constants->variante;

	grid_size[0] = constants->grid_size_x;
	grid_size[1] = constants->grid_size_y;
	gamma = constants->gamma;
	dx = (pos_x_max - pos_x_min) / (double) grid_size[0];
	dy = (pos_y_max - pos_y_min) / (double) grid_size[1];

	dt = 0.0;

	ct = constants->ct;

	with_splitting = 1;
	step_output = 0;
	/*if (dimension == 2) {
	 std::cout << "Wahl!" << endl << " 1: unsplitting, 2: splitting:";
	 std::cin >> with_splitting;
	 }*/
}

/**
 *****************************************************************************************
 * Destruktor
 *****************************************************************************************/
Numerical_Method::~Numerical_Method() {
	delete[] grid_size;

}

/**
 *****************************************************************************************
 * Startet die numerische Berechnung der Mehrphasenströmung.
 *****************************************************************************************/
void Numerical_Method::start_method() {
	double time = 0.0, timetol = 0.000001, timedif = 1.0;

	// Falls kein Schleifendurchgang gemacht wird
	if (step_output == 1)
		write();

	for (int n = 1; n <= step_limit && timedif > timetol; n++) {
		cout << n << " : " << step_limit << endl;

		// set boundary conditions

		grid_main->apply_boundary_conditions();
		if (step_output == 1)
			write();

		// compute time step
		time = cfl_condition(n, time);

		// update
		if (dimension == 1)
			solver->calc_method_flux(dt, 1);

		else if (with_splitting == 1)
			solver->calc_method_flux(dt, 0);

		else if (with_splitting == 2)
			solver->calc_method_flux(dt, 1);

		timedif = fabs(time - time_limit);
		steps = n;

	}
	write();
}

/**
 *****************************************************************************************
 * CFL Bedingung anwenden und neue Zeit berechnen.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 * @return neue Zeit.
 *****************************************************************************************/
double Numerical_Method::cfl_condition(int n, double time) {

	switch (dimension) {
	case (1): {
		//1D
		if ((int) constants->calceigv == 1)
		//1 = true, smax wird über Eigenwerte bestimmt
		{
			cfl_1d_eigenvalues(n, time);
		}

		else {
			cfl_1d_approx(n, time);
		}

		break;

	}
	case (2): {
		//2 Dimensionen

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, smax wird über Eigenwerte bestimmt

		if ((int) constants->calceigv == 0) {
			cfl_2d_approx(n, time);
		}

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, smax wird über Eigenwerte bestimmt
		else {
			cfl_2d_eigenvalues(n, time);
		}

		break;
	}
	}
	return time;
}

/**
 *****************************************************************************************
 *  Schreibt Ergebnisse in Dateien für d,p und u,ur für die jeweilige Dimension
 *****************************************************************************************/
void Numerical_Method::write() {
	double xpos = 0.0;
	double ypos = 0.0;
	double p = 0.0;
	string added;

	if (dimension == 1) {
		added = to_string(grid_size[0]) + "_" + solver_name + "_" + to_string(dimension) + "d_" + to_string(cfl_option) + ".variant_" + "div" + to_string(divider)
				+ "till" + to_string((int) divider_end) + "_" + to_string(steps) + "Steps";

	} else {
		added = to_string(grid_size[0]) + "x" + to_string(grid_size[1]) + "_" + solver_name + "_" + to_string(dimension) + "d_split" + to_string(with_splitting)
				+ "_IC" + to_string(grid_main->choice) + "_div" + to_string(int(1. / divider)) + "till" + to_string((int) divider_end) + "_" + to_string(steps)
				+ "Steps";
	}

	switch (grid_main->dimension) {
	case (1): {
		string d_path = "d_" + added;
		string uxr_path = "uxr_" + added;
		string ux_path = "ux_" + added;
		string p_path = "p_" + added;

		ofstream d_out(d_path.c_str());
		ofstream uxr_out(uxr_path.c_str());
		ofstream ux_out(ux_path.c_str());
		ofstream p_out(p_path.c_str());

		for (int i = order; i < grid_main->grid_size_total[0] - grid_main->orderofgrid; i++) {
			xpos = (pos_x_min + (pos_x_max - pos_x_min) * ((double) (i - order) / grid_size[0]));

			p = ct * pow(grid_main->cellsgrid[i][0], gamma);

			d_out << fixed << setprecision(8) << xpos << " \t" << setprecision(10) << grid_main->cellsgrid[i][0] << "\n";
			uxr_out << fixed << setprecision(8) << xpos << " \t" << setprecision(10) << grid_main->cellsgrid[i][3] << "\n";
			ux_out << fixed << setprecision(8) << xpos << " \t" << setprecision(10) << grid_main->cellsgrid[i][2] << "\n";
			p_out << fixed << setprecision(8) << xpos << " \t" << scientific << setprecision(10) << p << "\n";

		}

		d_out.close();
		uxr_out.close();
		ux_out.close();
		p_out.close();
		break;
	}
	case (2): {
		//2D

		// output files fuer das 2-d Feld
		string d_path = "d_" + added;
		string uxr_path = "uxr_" + added;
		string uyr_path = "uyr_" + added;
		string ux_path = "ux_" + added;
		string uy_path = "uy_" + added;
		string p_path = "p_" + added;
		ofstream d_out(d_path.c_str());
		ofstream uxr_out(uxr_path.c_str());
		ofstream ux_out(ux_path.c_str());
		ofstream uyr_out(uyr_path.c_str());
		ofstream uy_out(uy_path.c_str());
		ofstream p_out(p_path.c_str());

		// ACHTUNG, FUER 2. ORDNUNG NOCHMAL ALLE GRENZEN KONTROLLIEREN

		// Output ohne Randwerte
		for (int x = order; x < grid_main->grid_size_total[0] - grid_main->orderofgrid; x++) {
			for (int y = order; y < grid_main->grid_size_total[1] - grid_main->orderofgrid; y++) {
				xpos = pos_x_min + (pos_x_max - pos_x_min) * (((double) x - order) / (double) grid_size[0]);
				ypos = pos_y_min + (pos_y_max - pos_y_min) * (((double) y - order) / (double) grid_size[1]);

				int index = x + y * grid_main->grid_size_total[0];
				p = ct * pow(grid_main->cellsgrid[index][0], gamma);

				d_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << setprecision(10)
						<< grid_main->cellsgrid[index][0] << "\n";
				uxr_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << setprecision(10)
						<< grid_main->cellsgrid[index][3] << "\n";
				ux_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << setprecision(10)
						<< grid_main->cellsgrid[index][2] << "\n";
				uyr_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << setprecision(10)
						<< grid_main->cellsgrid[index][5] << "\n";
				uy_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << setprecision(10)
						<< grid_main->cellsgrid[index][4] << "\n";
				p_out << fixed << setprecision(8) << xpos << " " << fixed << setprecision(8) << ypos << " \t" << scientific << setprecision(10) << p << "\n";
			}
		}

		// Dateien schliessen
		d_out.close();
		uxr_out.close();
		ux_out.close();
		uyr_out.close();
		uy_out.close();
		p_out.close();

		// output files fuer eine 1-d Linie im 2-d Feld
		string d_path_d1 = "d_" + added + "_d1";
		string u_path_d1 = "u_" + added + "_d1";
		string ur_path_d1 = "ur_" + added + "_d1";
		string p_path_d1 = "p_" + added + "_d1";
		ofstream d_out_d1(d_path_d1.c_str());
		ofstream u_out_d1(u_path_d1.c_str());
		ofstream ur_out_d1(ur_path_d1.c_str());
		ofstream p_out_d1(p_path_d1.c_str());

		//d => Abstand 1d betrachtung x-y (punkte einer geraden~)
		double d = 0.0, xout, uout, urout, ux, uy, sign;

		for (int x = order; x < grid_main->grid_size_total[0] - grid_main->orderofgrid; x++) {
			for (int y = order; y < grid_main->grid_size_total[1] - grid_main->orderofgrid; y++) {
				if (grid_main->choice < 3) {
					xpos = pos_x_min + (pos_x_max - pos_x_min) * (((double) x - order) / (double) grid_size[0]);
					ypos = pos_y_min + (pos_y_max - pos_y_min) * (((double) y - order) / (double) grid_size[1]);

					int index = x + y * grid_main->grid_size_total[0];

					p = ct * pow(grid_main->cellsgrid[index][0], gamma);

					if (grid_main->choice == 0) {
						d = fabs(ypos);
					}
					if (grid_main->choice == 1) {
						d = fabs(xpos + ypos);
					}
					if (grid_main->choice == 2) {
						d = fabs(xpos + 2 * ypos);
					}

					if (d < 0.0001) {
						sign = xpos < 0 ? -1.0 : 1.0;
						xout = sign * sqrt(xpos * xpos + ypos * ypos);
						ux = grid_main->cellsgrid[index][2];
						uy = grid_main->cellsgrid[index][4];
						sign = ux < 0 ? -1.0 : 1.0;
						uout = sign * sqrt(ux * ux + uy * uy);
						ux = grid_main->cellsgrid[index][3];
						uy = grid_main->cellsgrid[index][5];
						sign = ux < 0 ? -1.0 : 1.0;
						urout = sign * sqrt(ux * ux + uy * uy);

						d_out_d1 << fixed << setprecision(8) << xout << " \t" << setprecision(10) << grid_main->cellsgrid[index][0] << "\n";
						u_out_d1 << fixed << setprecision(8) << xout << " \t" << setprecision(10) << uout << "\n";
						ur_out_d1 << fixed << setprecision(8) << xout << " \t" << setprecision(10) << urout << "\n";
						p_out_d1 << fixed << setprecision(8) << xout << " \t" << scientific << setprecision(10) << p << "\n";
					}
				}
			}
		}
		// Dateien schliessen
		d_out_d1.close();
		u_out_d1.close();
		ur_out_d1.close();
		p_out_d1.close();
	}
	}
}

/**
 *****************************************************************************************
 *  Bestimmt die Elemente der Jacobi-Matrix für eine Dimension
 *****************************************************************************************/

//TODO: Parameter kürzen
void Numerical_Method::matrix_1d(double * values, int n, double * u, double p, double dtwo, int variante) {
	double cref = constants->cref;
	double ccl = constants->ccl;
	double done = constants->rho_one;

	switch (variante) {
	case (1):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[2] * u[2] + gamma * ct * pow(u[0], gamma - 1.0);
		values[4] = (2 * u[1]) / u[0];
		values[5] = u[0] * ccl * (1.0 - ccl) * 2.0 * u[2];
		values[6] = -u[1] / (u[0] * u[0]) * u[2] + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / done;
		values[7] = u[2] / u[0];
		values[8] = (u[1] / u[0]) + (1.0 - 2.0 * ccl) * u[2];
		break;
	case (2):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(u[1] * u[1]) / (u[0] * u[0]) + ccl * (1 - ccl) * u[2] * u[2]
				+ cref * gamma * pow(done, gamma - 1.0) * (ccl * done * done) / pow(done + (ccl - 1.0) * u[0], 2.0);
		values[4] = (2 * u[1]) / u[0];
		values[5] = u[0] * ccl * (1.0 - ccl) * 2.0 * u[2];
		values[6] = -u[1] / (u[0] * u[0]) * u[2]
				+ cref * (gamma * pow(dtwo, gamma - 2.0) - (gamma * pow(dtwo, gamma - 1.0)) / done) * (ccl * done * done) / pow(done + (ccl - 1.0) * u[0], 2.0);
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
		values[6] = -u[1] / (u[0] * u[0]) * u[2] + (1.0 / ccl) * ((1.0 / u[0]) - (1.0 / done)) * ct * gamma * pow(u[0], gamma - 1.0);
		values[7] = u[2] / u[0];
		values[8] = (u[1] / u[0]) + (1.0 - 2.0 * ccl) * u[2];
	}
}

/**
 *****************************************************************************************
 *  Bestimmt die Elemente der Jacobi-Matrix für zwei Dimensionen
 *****************************************************************************************/
void Numerical_Method::matrix_2d(double * values_x, double * values_y, int n, double * u, double p,
		double dtwo) {

	double cref = constants->cref;
	double ccl = constants->ccl;
	double done = constants->rho_one;

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
	values_x[15] = -(u[1] * u[3]) / (u[0] * u[0]) + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / (done);
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
	values_y[20] = -(u[2] * u[4]) / (u[0] * u[0]) + pow(cref / ct, 1.0 / gamma) * gamma * ct * pow(u[0], gamma - 2.0) - (gamma * ct * pow(u[0], gamma - 1.0)) / (done);
	values_y[21] = 0;
	values_y[22] = u[4] / u[0];
	values_y[23] = 0;
	values_y[24] = u[2] / u[0] + (1 - 2 * ccl) * u[4];
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D mit Eigenwerten.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Numerical_Method::cfl_1d_eigenvalues(int n, double &time) {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double gi = 1.0 / gamma;

	double* max_u = new double[dimension];
	double* max_ur = new double[dimension];

	max_u[0] = 0.0, max_ur[0] = 0.0, max_u[1] = 0.0, max_ur[1] = 0.0;
	double p = 0.0, dtwo = 0.0;

	//TODO: max_u, max_u_step
	double smax = 0.0;

	int n_eqns;  // number of equations, up in a first line in "formeln....in"

	//TODO: Nicht von Dimensionen abhängig
	if (dimension == 1)
		n_eqns = 3;
	else
		n_eqns = 5;

	double* u_eqns = new double[n_eqns];

	string line;

	double values[9];
	LaVectorDouble real(3);
	LaVectorDouble img(3);
	LaGenMatDouble vr(3, 3);

	//Schritt 1: Maxima finden
	for (int x = 0; x < grid_main->grid_size_total[0]; x++) {

		// ACHTUNG, NUR VARIANTE 1 IST IN DEN GLEICHUNGEN IMPLEMENTIERT!
		switch (cfl_option) {
		case (1):

			p = ct * pow(grid_main->cellsgrid[x][0], gamma);

			grid_main->cellsgrid[x][1] = p;
			dtwo = pow((p / cref), gi);
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

		u_eqns[0] = grid_main->cellsgrid[x][0];
		u_eqns[1] = grid_main->cellsgrid[x][2] * u_eqns[0];
		u_eqns[2] = grid_main->cellsgrid[x][3];

		//Schritt 2: Einsetzen in die Jacobi-Matrix
		matrix_1d(values, n, u_eqns, p, dtwo, cfl_option);

		//Schritt 3: Berechnen der Eigenwerte
		LaGenMatDouble A(values, 3, 3, true);
		LaEigSolve(A, real, img, vr);

		//Schritt 4: Höchsten Eigenwert suchen
		for (int n = 0; n < 3; n++) {
			if (smax < fabs(real(n)))
				smax = fabs(real(n));
		}
	}

	printf("using eig, computed in all cells: %10.4e\n", smax);

	dt = cfl * dx / smax;

	//TODO: else if, unteres if prior?
	if (n <= divider_end)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;

	delete[] u_eqns;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 1D über Annäherung.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Numerical_Method::cfl_1d_approx(int n, double &time) {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double ccl = constants->ccl;
	double done = constants->rho_one;
	double gi = 1.0 / gamma;

	double* max_u = new double[dimension];
	double* max_ur = new double[dimension];

	max_u[0] = 0.0, max_ur[0] = 0.0, max_u[1] = 0.0, max_ur[1] = 0.0;
	double d = 0.0, p = 0.0, uxr = 0.0, ux = 0.0, dtwo = 0.0;

	//TODO: max_u, max_u_step
	double smax = 0.0, maxs = 0.0;

	string line;

	//dt über Näherung berechnen
	for (int i = 0; i < grid_main->grid_size_total[0]; i++) {

		d = grid_main->cellsgrid[i][0];
		uxr = grid_main->cellsgrid[i][3];
		ux = grid_main->cellsgrid[i][2];

		switch (cfl_option) {
		case (1): {
			p = ct * pow(d, gamma);
			grid_main->cellsgrid[i][1] = p;

			dtwo = pow((p / cref), gi);
			break;
		}
		case (2): {
			dtwo = ccl / ((1 / d) - ((1 - ccl) / done));
			p = cref * pow(dtwo, gamma);
			grid_main->cellsgrid[i][1] = p;

			break;
		}
		case (3): {
			dtwo = ccl / ((1 / d) - ((1 - ccl) / done));
			p = ct * pow(d, gamma);
			grid_main->cellsgrid[i][1] = p;

			break;
		}
		}

		maxs = fabs(ux) + AM_1d(gamma, p, d) + XS_1d(d, dtwo, done, ccl, uxr);

		if (maxs > smax)
			smax = maxs;
	}
	printf("using lambda exact %10.4e\n", smax);

	dt = cfl * dx / smax;

	if (n <= divider_end)
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
void Numerical_Method::cfl_2d_eigenvalues(int n, double &time) {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double gi = 1.0 / gamma;

	double max_d = 0.0;
	double* max_u = new double[dimension];
	double* max_ur = new double[dimension];

	max_u[0] = 0.0, max_ur[0] = 0.0, max_u[1] = 0.0, max_ur[1] = 0.0;
	double p = 0.0, dtwo = 0.0;

	double smax1 = 0, smax2 = 0;

	int n_eqns;  // number of equations, up in a first line in "formeln....in"

	//TODO: Nicht von Dimensionen abhängig
	if (dimension == 1)
		n_eqns = 3;
	else
		n_eqns = 5;

	double* u_eqns = new double[n_eqns];

	string line;

	//ueber Eigenwerte

	int n_eqns_squared = n_eqns * n_eqns;
	double values1[n_eqns_squared];
	double values2[n_eqns_squared];

	/**** Eventuell als Unteraufruf, um dt überprüfen zu können bei Splitting ****/

	//Schritt 1: Daten holen
	for (int x = 0; x < grid_main->grid_size_total[0]; x++) {
		for (int y = 0; y < grid_main->grid_size_total[1]; y++) {

			max_d = grid_main->cellsgrid[x][0];
			max_u[0] = grid_main->cellsgrid[x][2];
			max_u[1] = grid_main->cellsgrid[x][4];
			max_ur[0] = grid_main->cellsgrid[x][3];
			max_ur[1] = grid_main->cellsgrid[x][5];

			p = ct * pow(grid_main->cellsgrid[x + y * grid_main->grid_size_total[0]][0], gamma);
			grid_main->cellsgrid[x + y * grid_main->grid_size_total[0]][1] = p;
			dtwo = pow((p / cref), gi);

			max_u[0] = max_u[0] * max_d;
			max_u[1] = max_u[1] * max_d;

			//Schritt 2: Einsetzen
			u_eqns[0] = max_d;
			u_eqns[1] = max_u[0];
			u_eqns[2] = max_u[1];
			u_eqns[3] = max_ur[0];
			u_eqns[4] = max_ur[1];
			matrix_2d(values1, values2, n_eqns_squared, u_eqns, p, dtwo);

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
				if (smax1 < fabs(real(n)))
					smax1 = fabs(real(n));
				if (smax2 < fabs(real_b(n)))
					smax2 = fabs(real_b(n));
			}
		}
	}

	/*****/
	if (with_splitting == 1) {
		dt = cfl / (smax1 / dx + smax2 / dy);
	} else {
		dt = cfl / max(smax1 / dx, smax2 / dy);
	}

	if (n <= divider_end)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Größte Eigenwerte: " << smax1 << " und " << smax2 << endl;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;

	delete[] u_eqns;
}

/**
 *****************************************************************************************
 * Unterfunktion der Methode cfl_condition().
 * Berechnung der CFL Bedingung in 2D über Annäherung.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 *****************************************************************************************/
void Numerical_Method::cfl_2d_approx(int n, double &time) {
	double cref = constants->cref;
	double cfl = constants->cfl;
	double ccl = constants->ccl;
	double done = constants->rho_one;
	double gi = 1.0 / gamma;

	double* max_u = new double[dimension];
	double* max_ur = new double[dimension];

	max_u[0] = 0.0, max_ur[0] = 0.0, max_u[1] = 0.0, max_ur[1] = 0.0;
	double d = 0.0, p = 0.0, uxr = 0.0, ux = 0.0, uy = 0.0, uyr = 0.0, dtwo = 0.0;

	//TODO: max_u, max_u_step
	double smax = 0.0, maxs = 0.0;

	string line;

	//über Näherung
	for (int x = 0; x < grid_main->grid_size_total[0]; x++) {
		for (int y = 0; y < grid_main->grid_size_total[1]; y++) {

			int index = x + y * grid_main->grid_size_total[0];
			d = grid_main->cellsgrid[index][0];
			uxr = grid_main->cellsgrid[index][3];
			ux = grid_main->cellsgrid[index][2];
			uy = grid_main->cellsgrid[index][4];
			uyr = grid_main->cellsgrid[index][5];
			p = ct * pow(d, gamma);
			grid_main->cellsgrid[x + y * grid_main->grid_size_total[0]][1] = p;
			dtwo = pow((p / cref), gi);

			maxs = fabs(ux) + fabs(uy) + AM_2d(gamma, p, d) + XS_2d(d, dtwo, done, ccl, uxr, uxy);

			if (maxs > smax) {
				smax = maxs;
			}
		}
	}

	cout << "Größte geratene Eigenwerte: " << smax << endl;

	dt = cfl * min(dx, dy) / smax;

	if (n <= divider_end)
		dt = dt * divider;

	if ((time + dt) > time_limit)
		dt = time_limit - time;
	time = time + dt;
	cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t" << time << endl;

}
