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
#include "adaptive_mesh.h"

#include "time_step_calculation.h"

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
	time_calculation = new Time_Step_Calculation(computation->neqs,grid_main);
	solver->time_calculation=time_calculation;

	dimension = constants->dimension;
	grid_size = new int[dimension];

	pos_x_max = constants->pos_x_max;
	pos_x_min = constants->pos_x_min;
	pos_y_max = constants->pos_y_max; //2D
	pos_y_min = constants->pos_y_min; //2D
	time_limit = constants->timeou;
	steps = 0;
	step_limit = constants->maxnt;
	divider_last = constants->teilerend;
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

	output_per_step = 0;

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
	delete time_calculation;
}

/**
 *****************************************************************************************
 * Startet die numerische Berechnung der Mehrphasenströmung.
 *****************************************************************************************/
void Numerical_Method::start_method() {
	double time = 0.0, timetol = 0.000001, timedif = 1.0;
	time_calculation->with_splitting = with_splitting;

	// Falls kein Schleifendurchgang gemacht wird
	if (output_per_step == 1)
		write();

	for (int n = 1; n <= step_limit && timedif > timetol; n++) {
	//for (int n = 1; n < 500; n++) {
		cout << n << " : " << step_limit << endl;

		// set boundary conditions

		grid_main->apply_boundary_conditions();
		if (output_per_step == 1)
			write();

		if(10){
			Adaptive_Mesh* adaptive_mesh = new Adaptive_Mesh(solver, grid_main, constants, time_calculation);
			adaptive_mesh->amr();
			delete adaptive_mesh;
			exit(0);
		}

		// compute time step
		time_calculation->n = n;
		time_calculation->time = time;

		dt = time_calculation->cfl_condition(grid_main);
		cout << "Current time " << time << endl;

		// update
		if (dimension == 1){
			int spliting_method = 1;
			solver->split_method = spliting_method;
			solver->calc_method_flux(dt,grid_main);
		}
		else if (dimension == 2){
			int spliting_method = with_splitting - 1;
			solver->split_method = spliting_method;
			solver->calc_method_flux(dt,grid_main);
		}
		if(time != time_calculation->time)
			time = time_calculation->time;

		timedif = fabs(time - time_limit);
		steps = n;

	}
	write();
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
				+ "till" + to_string((int) divider_last) + "_" + to_string(steps) + "Steps";

	} else {
		added = to_string(grid_size[0]) + "x" + to_string(grid_size[1]) + "_" + solver_name + "_" + to_string(dimension) + "d_split" + to_string(with_splitting)
				+ "_IC" + to_string(grid_main->choice) + "_div" + to_string(int(1. / divider)) + "till" + to_string((int) divider_last) + "_" + to_string(steps)
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
