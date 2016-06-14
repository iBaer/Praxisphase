#include "lax_friedrich.h"
#include <stdlib.h> // exit

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
	for (int i = 0; i < dimension; i++) {
		f_lax[i] = new double[neqs * (size_m1[0]) * (size_m1[1])];
	}

	if(dimension==2){
		split_grid = new Grid*[2];
		split_grid[0] = new Grid(size_total[0] * size_total[1]);
		split_grid[1] = new Grid(size_total[0] * size_total[1]);

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
void Lax_Friedrich::calc_method_flux(double dt, int split_method) {
	cout << "Lax-Friedrich Fluss berechnen..." << endl;

	switch (dimension) {
	// Eine Dimension
	case (1): {
		solve_1d(dt);
		break;
	}

		// Zwei Dimensionen
	case (2): {

		if (split_method == 0) {
			solve_2d_unsplit(dt);
			break;
		}

		else if (split_method == 1) {
			solve_2d_split_xtoy(dt,0,0);
		}
		else if (split_method == 2) {
			solve_2d_split_ytox(dt,0);
		}
		else {
			solve_2d_split_xtoy(dt,1,0);
			solve_2d_split_ytox(dt,1);
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
void Lax_Friedrich::solve_1d(double dt) {
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
					0.5 * (f[k][i][0] + f[k][i + 1][0])
					+ 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
		}
	}

	int index_end = (height_m1) * (width_m1);

	// Updateschritt
	for (int i = order; i < grid->grid_size_total[0] - grid->orderofgrid; i++) {

		// Index-Form: index_x1/2 + index_end * (VARIABLE (D=0,UXD=1, ...))
		int index_x1 = 0 + (height_m1) * (i - 1);
		int index_x2 = 0 + (height_m1) * (i);

		d = grid->cellsgrid[i][0];
		ux = grid->cellsgrid[i][2];
		uxd = d * ux;
		uxr = grid->cellsgrid[i][3];

		d = d + dtodx * (f_lax[0][index_x1 + index_end * (0)] - f_lax[0][index_x2 + index_end * (0)]);
		uxd = uxd + dtodx * (f_lax[0][index_x1 + index_end * (1)] - f_lax[0][index_x2 + index_end * (1)]);
		uxr = uxr + dtodx * (f_lax[0][index_x1 + index_end * (2)] - f_lax[0][index_x2 + index_end * (2)]);

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
void Lax_Friedrich::solve_2d_unsplit(double dt) {
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
	int index = 0;
	for (int k = 0; k < computation->neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {

				//TODO: Index weiter splitten
				index = y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k;

				f_lax[0][index] =
						0.5 * (f[k][x][y] + f[k][x + 1][y])
						+ 0.5 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);	//	0.25~ -> sonst 60 Grad unsplitting falsch

				f_lax[1][index] =
						0.5 * (g[k][x][y] + g[k][x][y + 1])
						+ 0.5 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);	// 0.25~

			}
		}
	}

	cout << "update x mit dtodx=" << dtodx << endl;

	int index_x1 = 0, index_x2 = 0, index_x3 = 0, index_y1 = 0, index_y2 = 0, index_y3 = 0;

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

			//Index könnte noch über Schleife (in X-Teile / Y-Teile) geteilt und optimiert werden
			index_x1 = y + (height_m1) * (x - 1);
			index_x2 = y + (height_m1) * (x);
			index_x3 = (height_m1) * (width_m1);

			index_y1 = (y - 1) + (height_m1) * (x);
			index_y2 = (y) + (height_m1) * (x);
			index_y3 = (height_m1) * (width_m1);

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
void Lax_Friedrich::solve_2d_split_xtoy(double dt, int with_average, int rerun) {
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;
	Grid* set_grid;

	// X-Richtung
	cout << "update x mit dtodx=" << dtodx << endl;

	// Falls xtoy und ytox gemittelt werden sollen, Werte in einem seperaten Grid setzen um anschließend ins Hauptgrid "zu mitteln"
	if(with_average==1){
		set_grid = split_grid[0];
	}
	else {
		set_grid = grid;
	}

	//Berechne U, F
	computation->compute_u_2d(cs, grid);
	computation->compute_f_2d(f, grid);

	// Berechne Lax-Friedrich Flüsse
	// Faktor 0.25 nach Formel für unsplitting,
	// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting

	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
				f_lax[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] =
						0.5 * (f[k][x][y] + f[k][x + 1][y])
						+ 0.25 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]); // 0.25~
			}
		}
	}

	int index_x1 = 0, index_x2 = 0, index_x_end = (height_m1) * (width_m1);

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

			index_x1 = y + (height_m1) * (x - 1);
			index_x2 = y + (height_m1) * (x);

			d = d + dtodx * (f_lax[0][index_x1 + index_x_end * (0)] - f_lax[0][index_x2 + index_x_end * (0)]);
			uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x_end * (1)] - f_lax[0][index_x2 + index_x_end * (1)]);
			uyd = uyd + dtodx * (f_lax[0][index_x1 + index_x_end * (2)] - f_lax[0][index_x2 + index_x_end * (2)]);
			uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x_end * (3)] - f_lax[0][index_x2 + index_x_end * (3)]);
			uyr = uyr + dtodx * (f_lax[0][index_x1 + index_x_end * (4)] - f_lax[0][index_x2 + index_x_end * (4)]);

			set_grid->cellsgrid[pos][0] = d;
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;
		}
	}

	grid->apply_boundary_conditions();

	// ACHTUNG: HIER SOLLTE UEBERPRUEFT WERDEN, OB DER
	// ZEITSCHRITT NICHT ZU GROSS IST FUER DIE 2 RICHTUNG MIT
	// DEN NEUEN WERTEN, SONST KANN EINEM DAS SYSTEM DIVERGIEREN!


	//TODO: Neues Delta T und Eigenwerte überprüfen
	//Neues Delta T darf nicht kleiner sein, als vorher

	if (rerun == 0) {
		time_calculation->cfl_condition();
		int dt_comp = time_calculation->compare_dt();

		if (dt_comp == -1) {
			dt = time_calculation->halve_dt();

			cout << "New delta t is smaller than before!" << endl;
			cout << "Restarting update with dt/2!" << endl;

			solve_2d_split_xtoy(dt, with_average,1);

			cout << "Finished recalculation!"<<endl;
			return;

		} else {
			// go on
		}
	}

	// Y-Richtung
	cout << "update y mit dtody=" << dtody << endl;


	//Berechne U und G
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
						+ 0.25 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]); // 0.25~
			}
		}
	}
	int index_y1 = 0, index_y2 = 0, index_y_end = (height_m1) * (width_m1);

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

			index_y1 = (y - 1) + (height_m1) * (x);
			index_y2 = (y) + (height_m1) * (x);

			d = d + dtody * (f_lax[1][index_y1 + index_y_end * (0)] - f_lax[1][index_y2 + index_y_end * (0)]);
			uxd = uxd + dtody * (f_lax[1][index_y1 + index_y_end * (1)] - f_lax[1][index_y2 + index_y_end * (1)]);
			uyd = uyd + dtody * (f_lax[1][index_y1 + index_y_end * (2)] - f_lax[1][index_y2 + index_y_end * (2)]);
			uxr = uxr + dtody * (f_lax[1][index_y1 + index_y_end * (3)] - f_lax[1][index_y2 + index_y_end * (3)]);
			uyr = uyr + dtody * (f_lax[1][index_y1 + index_y_end * (4)] - f_lax[1][index_y2 + index_y_end * (4)]);

			set_grid->cellsgrid[pos][0] = d;
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;

		}
	}
}

/**
 *****************************************************************************************
 * Untermethode, die die Berechnung des Lax-Friedrich-Flusses
 * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
 * @param dt Delta t.
 *****************************************************************************************/
void Lax_Friedrich::solve_2d_split_ytox(double dt, int with_average) {
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

	// X-Richtung
	cout << "update x mit dtodx=" << dtodx << endl;


	//Berechne U, F
	computation->compute_u_2d(cs, grid);
	computation->compute_f_2d(f, grid);

	// Berechne Lax-Friedrich Flüsse
	// Faktor 0.25 nach Formel für unsplitting,
	// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting

	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
				f_lax[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] =
						0.5 * (f[k][x][y] + f[k][x + 1][y])
						+ 0.5 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]); // 0.25~
			}
		}
	}

	int index_x1 = 0, index_x2 = 0, index_x_end = (height_m1) * (width_m1);

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

			index_x1 = y + (height_m1) * (x - 1);
			index_x2 = y + (height_m1) * (x);

			d = d + dtodx * (f_lax[0][index_x1 + index_x_end * (0)] - f_lax[0][index_x2 + index_x_end * (0)]);
			uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x_end * (1)] - f_lax[0][index_x2 + index_x_end * (1)]);
			uyd = uyd + dtodx * (f_lax[0][index_x1 + index_x_end * (2)] - f_lax[0][index_x2 + index_x_end * (2)]);
			uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x_end * (3)] - f_lax[0][index_x2 + index_x_end * (3)]);
			uyr = uyr + dtodx * (f_lax[0][index_x1 + index_x_end * (4)] - f_lax[0][index_x2 + index_x_end * (4)]);

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


	// Y-Richtung
	cout << "update y mit dtody=" << dtody << endl;


	//Berechne U und G
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
						+ 0.5 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]); // 0.25~
			}
		}
	}
	int index_y1 = 0, index_y2 = 0, index_y_end = (height_m1) * (width_m1);

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

			index_y1 = (y - 1) + (height_m1) * (x);
			index_y2 = (y) + (height_m1) * (x);

			d = d + dtody * (f_lax[1][index_y1 + index_y_end * (0)] - f_lax[1][index_y2 + index_y_end * (0)]);
			uxd = uxd + dtody * (f_lax[1][index_y1 + index_y_end * (1)] - f_lax[1][index_y2 + index_y_end * (1)]);
			uyd = uyd + dtody * (f_lax[1][index_y1 + index_y_end * (2)] - f_lax[1][index_y2 + index_y_end * (2)]);
			uxr = uxr + dtody * (f_lax[1][index_y1 + index_y_end * (3)] - f_lax[1][index_y2 + index_y_end * (3)]);
			uyr = uyr + dtody * (f_lax[1][index_y1 + index_y_end * (4)] - f_lax[1][index_y2 + index_y_end * (4)]);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;

		}
	}
}
