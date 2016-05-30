#include "lax_friedrich.h"

using namespace std;

/**
 *****************************************************************************************
* Konstruktor von der Klasse Lax_Friedrich.
* Ruft den Konstruktor der Superklasse Solver auf.
* @param constants Pointer auf das Objekt, welches die Konstanten enthält.
* @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird.
* @param grid Pointer auf Raster-Objekt.
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

/**
 *****************************************************************************************
* Destruktor
 *****************************************************************************************/
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
* Implementierung der virtuellen Methode zur Berechnung des Flußes.
* Wählt lediglich das gewünschte Schema aus und delegiert die Berechnung weiter.
* @param dt Delta t.
* @param dir Unsplitting = 0, Splitting = 1.
 *****************************************************************************************/
void Lax_Friedrich::calc_method_flux(double dt, int dir) {
	cout << "Lax-Friedrich Fluss berechnen..." << endl;

	switch (dimension) {
	// Eine Dimension
	case (1): {
		solve_1d(dt);
		break;
	}

		// Zwei Dimensionen
	case (2): {

		if (dir == 0) {
			solve_2d_unsplit(dt);
			break;
		}

		else{
			solve_2d_split(dt);
		}
	}
	}

}

/**
 *****************************************************************************************
* Untermethode, die die Berechnung des Lax-Friedrich-Flusses
* und das Updaten der Zellen für eine Dimension durchführt.
* @param dt Delta t.
 *****************************************************************************************/
void Lax_Friedrich::solve_1d(double dt){
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double d, ux, uxd, uxr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;

	computation->compute_u_1d(cs, grid);
	computation->compute_f_1d(f, grid);
	//Berechne Lax-Friedrich-Fluss
	for (int k = 0; k < neqs; k++) {
		for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {
			f_lax[0][0 + (size_m1[1]) * i + (size_m1[1]) * (size_m1[0]) * k] =
					0.5 * (f[k][i][0] + f[k][i + 1][0]) + 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
		}
	}

	// Updateschritt
	for (int i = order; i < grid->grid_size_total[0] - grid->orderofgrid; i++) {
		d = grid->cellsgrid[i][0];
		ux = grid->cellsgrid[i][2];
		uxd = d * ux;
		uxr = grid->cellsgrid[i][3];

		// Form: index_x1 + index_x2 * (XPOS) + index_end * (VARIABLE (D=0,UXD=1, ...))
		// TODO: Form könnte für 1D Fall vereinfacht werden
		int index_x1 = 0 + (height_m1) * (i - 1);
		int index_x2 = 0 + (height_m1) * (i);
		int index_x3 = (height_m1) * (width_m1);

		d = d + dtodx * (f_lax[0][index_x1 + index_x3 * (0)] - f_lax[0][index_x2 + index_x3 * (0)]);
		uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x3 * (1)] - f_lax[0][index_x2 + index_x3 * (1)]);
		uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x3 * (2)] - f_lax[0][index_x2 + index_x3 * (2)]);

		grid->cellsgrid[i][0] = d;
		grid->cellsgrid[i][2] = uxd / d;
		grid->cellsgrid[i][3] = uxr;
	}
}

/**
 *****************************************************************************************
* Untermethode, die die Berechnung des Lax-Friedrich-Flusses
* und das Updaten der Zellen für zwei Dimension durchführt.
* @param dt Delta t.
 *****************************************************************************************/
void Lax_Friedrich::solve_2d_unsplit(double dt){
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

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

	cout << "update x mit dtodx=" << dtodx << endl;
	for (int x = order; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
		for (int y = order; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
			pos = x + y * grid->grid_size_total[0];
			d = grid->cellsgrid[pos][0];
			ux = grid->cellsgrid[pos][2];
			uy = grid->cellsgrid[pos][4];
			uxr = grid->cellsgrid[pos][3];
			uyr = grid->cellsgrid[pos][5];

			uxd = ux * d;
			uyd = uy * d;


			int index_x1 = y + (height_m1) * (x - 1);
			int index_x2 = y + (height_m1) * (x);
			int index_x3 = (height_m1) * (width_m1);

			int index_y1 = (y-1) + (height_m1) * (x);
			int index_y2 = (y) + (height_m1) * (x);
			int index_y3 = (height_m1) * (width_m1);

			d = d + dtodx * (f_lax[0][index_x1 + index_x3 * (0)] - f_lax[0][index_x2 + index_x3 * (0)])
					+ dtody * (f_lax[1][index_y1 + index_y3 * (0)] - f_lax[1][index_y2 + index_y3 * (0)]);
			uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x3 * (1)] - f_lax[0][index_x2 + index_x3 * (1)])
					+ dtody * (f_lax[1][index_y1 + index_y3 * (1)] - f_lax[1][index_y2 + index_y3 * (1)]);
			uyd = uyd + dtodx * (f_lax[0][index_x1 + index_x3 * (2)] - f_lax[0][index_x2 + index_x3 * (2)])
					+ dtody * (f_lax[1][index_y1 + index_y3 * (2)] - f_lax[1][index_y2 + index_y3 * (2)]);
			uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x3 * (3)] - f_lax[0][index_x2 + index_x3 * (3)])
					+ dtody * (f_lax[1][index_y1 + index_y3 * (3)] - f_lax[1][index_y2 + index_y3 * (3)]);
			uyr = uyr + dtodx * (f_lax[0][index_x1 + index_x3 * (4)] - f_lax[0][index_x2 + index_x3 * (4)])
					+ dtody * (f_lax[1][index_y1 + index_y3 * (4)] - f_lax[1][index_y2 + index_y3 * (4)]);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;
		}
	}
}

/**
 *****************************************************************************************
* Untermethode, die die Berechnung des Lax-Friedrich-Flusses
* und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
* @param dt Delta t.
 *****************************************************************************************/
void Lax_Friedrich::solve_2d_split(double dt){
	int order = grid->orderofgrid;
		double dtodx = dt / dx;
		double dtody = dt / dy;
		double d, ux, uy, uxd, uyd, uxr, uyr;
		int width_m1 = grid->grid_size_total[0] - 1;
		int height_m1 = grid->grid_size_total[1] - 1;
		int pos;

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

		cout << "update x mit dtodx=" << dtodx << endl;
		for (int x = order; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
			for (int y = order; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
				pos = x + y * grid->grid_size_total[0];
				d = grid->cellsgrid[pos][0];
				ux = grid->cellsgrid[pos][2];
				uy = grid->cellsgrid[pos][4];
				uxr = grid->cellsgrid[pos][3];
				uyr = grid->cellsgrid[pos][5];

				uxd = ux * d;
				uyd = uy * d;

				// Form: index_x1 + index_x2 * (XPOS) + index_end * (VARIABLE (D=0,UXD=1, ...))

				int index_x1 = y + (height_m1) * (x - 1);
				int index_x2 = y + (height_m1) * (x);
				int index_x3 = (height_m1) * (width_m1);
				//cout << "Dir 1 fi["<< index_x1 + index_x3 * (0)<<"]: "<<fi[index_x1 + index_x3 * (0)]<<endl;

				d = d + dtodx * (f_lax[0][index_x1 + index_x3 * (0)] - f_lax[0][index_x2 + index_x3 * (0)]);
				uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x3 * (1)] - f_lax[0][index_x2 + index_x3 * (1)]);
				uyd = uyd + dtodx * (f_lax[0][index_x1 + index_x3 * (2)] - f_lax[0][index_x2 + index_x3 * (2)]);
				uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x3 * (3)] - f_lax[0][index_x2 + index_x3 * (3)]);
				uyr = uyr + dtodx * (f_lax[0][index_x1 + index_x3 * (4)] - f_lax[0][index_x2 + index_x3 * (4)]);

				grid->cellsgrid[pos][0] = d;
				grid->cellsgrid[pos][2] = uxd / d;
				grid->cellsgrid[pos][4] = uyd / d;
				grid->cellsgrid[pos][3] = uxr;
				grid->cellsgrid[pos][5] = uyr;
			}
		}

		grid->apply_boundary_conditions();

		// ACHTUNG: HIER SOLLTE UEBERPRUEFT WERDEN, OB DER
		// ZEITSCHRITT NICHT ZU GROSS IST FUER DIE 2 RICHTUNG MIT
		// DEN NEUEN WERTEN, SONST KANN EINEM DAS SYSTEM DIVERGIEREN!

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

		cout << "update y mit dtody=" << dtody << endl;
		for (int x = order; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
			for (int y = order; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
				pos = x + y * grid->grid_size_total[0];
				d = grid->cellsgrid[pos][0];
				ux = grid->cellsgrid[pos][2];
				uy = grid->cellsgrid[pos][4];
				uxr = grid->cellsgrid[pos][3];
				uyr = grid->cellsgrid[pos][5];

				uxd = ux * d;
				uyd = uy * d;

				// Form: 1 + (dimension * (YPOS)) + index_y1 + index_y2 * (VARIABLE (D=0,UXD=1, ...))
				int index_y1 = (y-1) + (height_m1) * (x);
				int index_y2 = (y) + (height_m1) * (x);
				int index_y3 = (height_m1) * (width_m1);


				d = d + dtody * (f_lax[1][index_y1 + index_y3 * (0)] - f_lax[1][index_y2 + index_y3 * (0)]);
				uxd = uxd + dtody * (f_lax[1][index_y1 + index_y3 * (1)] - f_lax[1][index_y2 + index_y3 * (1)]);
				uyd = uyd + dtody * (f_lax[1][index_y1 + index_y3 * (2)] - f_lax[1][index_y2 + index_y3 * (2)]);
				uxr = uxr + dtody * (f_lax[1][index_y1 + index_y3 * (3)] - f_lax[1][index_y2 + index_y3 * (3)]);
				uyr = uyr + dtody * (f_lax[1][index_y1 + index_y3 * (4)] - f_lax[1][index_y2 + index_y3 * (4)]);

				grid->cellsgrid[pos][0] = d;
				grid->cellsgrid[pos][2] = uxd / d;
				grid->cellsgrid[pos][4] = uyd / d;
				grid->cellsgrid[pos][3] = uxr;
				grid->cellsgrid[pos][5] = uyr;

			}
		}
}
