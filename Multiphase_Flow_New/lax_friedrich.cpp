#include "lax_friedrich.h"

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor von der Klasse LaxFriedrichMethod.
 * Ruft einfach den Konstrukter von der geerbten Klasse auf.
 *****************************************************************************************/
Lax_Friedrich::Lax_Friedrich(Constants *constants, Computation *computation, Grid *grid) :
		Solver("Lax-Friedrich", constants, computation, grid) {

	//temporär
	size_total[1] = grid->grid_size_total[1];
	size_m1[1] = grid->grid_size_total[1] - 1;

	neqs = computation->neqs;

	uall = new double[neqs * size_total[0] * size_total[1]];
	fall = new double[neqs * size_total[0] * size_total[1]];
	gall = new double[neqs * size_total[0] * size_total[1]];

	cs = new double**[neqs];
	f = new double**[neqs];
	g = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[size_total[0]];
		f[i] = new double*[size_total[0]];
		g[i] = new double*[size_total[0]];
		for (int j = 0; j < size_total[0]; j++) {
			cs[i][j] = uall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f[i][j] = fall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g[i][j] = gall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);

		}
	}

	f_lax = new double*[dimension];
	for (int i=0;i<dimension;i++){
		f_lax[i] = new double[neqs * (size_m1[0]) * (size_m1[1])];
	}



}

Lax_Friedrich::~Lax_Friedrich() {
	delete[] uall;
	delete[] fall;
	delete[] gall;
	for (int i = 0; i < neqs; i++) {
		delete[] cs[i];
		delete[] f[i];
		delete[] g[i];
	}
	delete[] cs;
	delete[] f;
	delete[] g;

	for (int i = 0; i < dimension; i++) {
		delete[] f_lax[i];
	}
	delete[] f_lax;
}

/**
 *****************************************************************************************
 * Berechnung des Lax-Friedrich Flusses.
 * @return Das zurückgelieferte Objekt ist ein Vektor mit 4 Dimensionen
 *         (Formel,grid x-Koordinate,grid y-Koordinate, Flussrichtung)
 *
 * Zur Zeit wird hier noch viel zu viel berechnet, wenn der Fluss nur in eine Richtung
 * benötigt wird!!
 * Weiterhin müssen die Laufvariablen neu sortiert werden !!
 *****************************************************************************************/
double* Lax_Friedrich::calc_method_flux(double dt, int dir) {
	cout << "Lax-Friedrich Fluss berechnen..." << endl;


	switch (dimension) {
	// Eine Dimension
	case (1): {
		computation->compute_u_1d(cs, grid);
		computation->compute_f_1d(f, grid);
		//Berechne Lax-Friedrich-Fluss
		for (int k = 0; k < neqs; k++) {
			for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {
				f_lax[0][0 + (size_m1[1]) * i + (size_m1[1]) * (size_m1[0]) * k] =
						0.5 * (f[k][i][0] + f[k][i + 1][0]) + 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
			}
		}
		break;
	}

		// Zwei Dimensionen
	case (2): {

		if (dir == 0) {

			//Berechne U, F und G
			computation->compute_u_2d(cs, grid);

			computation->compute_f_2d(f, grid);

			computation->compute_g_2d(g, grid);

			// Berechne Lax-Friedrich Flüsse
			// Faktor 0.25 nach Formel für unsplitting,
			// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting
			//int index = 0;
			for (int k = 0; k < computation->neqs; k++) {
				for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
					for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {

						f_lax[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] = 0.5
								* (f[k][x][y] + f[k][x + 1][y])

						//+ 0.25*(dx/dt)*(cs[k][x][y] - cs[k][x+1][y]);
								+ 0.5 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);

						f_lax[1][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] =

						0.5 * (g[k][x][y] + g[k][x][y + 1])
						//+ 0.25*(dy/dt)*(cs[k][x][y] - cs[k][x][y+1]);
								+ 0.5 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);

					}
				}
			}

			break;
		}

		else if (dir == 1) {

			//Berechne U, F und G
			computation->compute_u_2d(cs, grid);
			computation->compute_f_2d(f, grid);

			// Berechne Lax-Friedrich Flüsse
			// Faktor 0.25 nach Formel für unsplitting,
			// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting
			//int index = 0;

			for (int k = 0; k < neqs; k++) {
				for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
					for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
						f_lax[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] =
							0.5 * (f[k][x][y] + f[k][x + 1][y])
								//+ 0.25*(dx/dt)*(cs[k][x][y] - cs[k][x+1][y]);
								+ 0.5 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);
					}
				}
			}

			//break;
			return f_lax[0];

		} else {

			//Berechne U, F und G
			computation->compute_u_2d(cs, grid);
			computation->compute_g_2d(g, grid);

			// Berechne Lax-Friedrich Flüsse
			// Faktor 0.25 nach Formel für unsplitting,
			// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting
			//int index = 0;

			for (int k = 0; k < neqs; k++) {
				for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
					for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {

						f_lax[1][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] =
								0.5 * (g[k][x][y] + g[k][x][y + 1])
									//+ 0.25*(dy/dt)*(cs[k][x][y] - cs[k][x][y+1]);
									+ 0.5 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);
					}
				}
			}
			return f_lax[1];

		}
	}
	}

	return f_lax[0];
}
