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

	//TODO: Konstruktor auch bei Force für 1D optimieren

	allocate_cache(grid);
	/*size_total[1] = grid->grid_size_total[1];
	size_m1[1] = grid->grid_size_total[1] - 1;

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
		split_grid[0] = new Grid(size_total[0],size_total[1], constants);
		split_grid[1] = new Grid(size_total[0],size_total[1], constants);
	}

	with_halved_dt = 0;

	set_grid = grid;*/
}

/**
 *****************************************************************************************
 * Destruktor
 *****************************************************************************************/
Lax_Friedrich::~Lax_Friedrich() {
	delete_cache();
}

void Lax_Friedrich::delete_cache() {
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
	if(dimension==2){
		delete split_grid[0];
		delete split_grid[1];
		delete[] split_grid;
	}
}

void Lax_Friedrich::allocate_cache(Grid * grid) {
	// 1D Optimierung
	size_total[0] = grid->grid_size_total[0];
	size_m1[0] = grid->grid_size_total[0] - 1;
	size_total[1] = grid->grid_size_total[1];
	size_m1[1] = grid->grid_size_total[1] - 1;

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
		split_grid[0] = new Grid(size_total[0],size_total[1], constants);
		split_grid[1] = new Grid(size_total[0],size_total[1], constants);
	}

	with_halved_dt = 0;

	set_grid = grid;
}


/**
 *****************************************************************************************
 * Implementierung der virtuellen Methode zur Berechnung des Flußes.
 * Wählt lediglich das gewünschte Schema aus und delegiert die Berechnung weiter.
 * @param dt Delta t.
 * @param dir Unsplitting = 0, Splitting = 1.
 *****************************************************************************************/
void Lax_Friedrich::calc_method_flux(double dt, Grid * grid) {
	cout << "Lax-Friedrich Fluss berechnen..." << endl;
	if (this->grid!=grid){
		cout << "New grid! Reallocating cache!"<<endl;
		this->grid = grid;
		delete_cache();
		allocate_cache(grid);
	}
	else{
		cout << "Same grid as before, using old cache!"<<endl;

	}
	cout << "dt="<<dt<<endl;
	/*dt = time_calculation->halve_dt(dt);
	cout << "dt="<<dt<<endl;*/
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
			//time_calculation->set_new_time(dt);

			break;
		}

		else if (split_method == 1) {
			//dt = time_calculation->halve_dt();
			//solve_2d_split_xtoy(dt, 0, 1);
			solve_2d_split_xtoy(dt, 0, 0);
			if (with_halved_dt==1){
				dt = time_calculation->halve_dt(dt);
				with_halved_dt=0;
			}
			//time_calculation->set_new_time(dt);

		}
		else if (split_method == 2) {
			solve_2d_split_ytox(dt, 0, 0);
			if (with_halved_dt==1){
				dt = time_calculation->halve_dt(dt);
				with_halved_dt=0;
			}
			//time_calculation->set_new_time(dt);
		}
		else {
			//int dt_tmp;

			Grid* grid_x_y;
			Grid* grid_y_x;
			grid_x_y = solve_2d_split_xtoy(dt, 1, 0);
			//dt_tmp = dt;

			if (with_halved_dt==1) {
				dt = time_calculation->halve_dt(dt);
				grid_y_x = solve_2d_split_ytox(dt, 1, 1);
				with_halved_dt = 0;
			}
			else{
				time_calculation->reset_step();
				grid_y_x = solve_2d_split_ytox(dt, 1, 0);
				if (with_halved_dt==1) {
					dt = time_calculation->halve_dt(dt);
					grid_x_y = solve_2d_split_xtoy(dt, 1, 1);
					with_halved_dt = 0;
				}
			}
			//time_calculation->set_new_time(dt);

			split_mean(grid_x_y, grid_y_x);

		}
	}
	}
	time_calculation->set_new_time(dt);

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
		grid->cellsgrid[i][1] = constants->ct * pow(grid->cellsgrid[i][0], constants->gamma);
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
						+ 0.25 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);	//	0.25~ -> sonst 60 Grad unsplitting falsch

				f_lax[1][index] =
						0.5 * (g[k][x][y] + g[k][x][y + 1])
						+ 0.25 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);	// 0.25~

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
			grid->cellsgrid[pos][1] = constants->ct * pow(grid->cellsgrid[pos][0], constants->gamma);
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
Grid* Lax_Friedrich::solve_2d_split_xtoy(double dt, int with_average, int rerun) {
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

	// X-Richtung
	cout << "update x mit dtodx=" << dtodx << endl;

	set_grid = split_grid[0];

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

	/*for (int y = 0; y < grid->grid_size_total[1]; y++) {
		pos = y * grid->grid_size_total[0];
		for(int k=0;k<6;k++){
			set_grid->cellsgrid[pos][k] = grid->cellsgrid[pos][k];
		}
	}
	for (int y = 0; y < grid->grid_size_total[1]; y++) {
		pos = (grid->grid_size_total[0] - grid->orderofgrid) + y * grid->grid_size_total[0];
		for(int k=0;k<6;k++){
			set_grid->cellsgrid[pos][k] = grid->cellsgrid[pos][k];
		}
	}
	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		pos = x;
		for(int k=0;k<6;k++){
			set_grid->cellsgrid[pos][k] = grid->cellsgrid[pos][k];
		}
	}
	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		pos = x + (grid->grid_size_total[1] - grid->orderofgrid) * grid->grid_size_total[0];
		for(int k=0;k<6;k++){
			set_grid->cellsgrid[pos][k] = grid->cellsgrid[pos][k];
		}
	}*/

	// TODO: shorten?
	if(grid->copy_to(set_grid) == -1){
		cout << "Couldn't copy Grid";
		exit(EXIT_FAILURE);
	}

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
			set_grid->cellsgrid[pos][1] = constants->ct * pow(d, constants->gamma);
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;

		}
	}

	set_grid->apply_boundary_conditions();


	// ACHTUNG: HIER SOLLTE UEBERPRUEFT WERDEN, OB DER
	// ZEITSCHRITT NICHT ZU GROSS IST FUER DIE 2 RICHTUNG MIT
	// DEN NEUEN WERTEN, SONST KANN EINEM DAS SYSTEM DIVERGIEREN!


	//TODO: Neues Delta T und Eigenwerte überprüfen
	//Neues Delta T darf nicht kleiner sein, als vorher

	if (rerun == 0) {
		time_calculation->cfl_condition(set_grid);
		int dt_comp = time_calculation->compare_dt();

		if (dt_comp == -1) {
			dt = time_calculation->halve_dt(dt);
			with_halved_dt = 1;

			cout << "New delta t is smaller than before!" << endl;
			cout << "Restarting update with dt/2!" << endl;

			set_grid = solve_2d_split_xtoy(dt, with_average,1);

			cout << "Finished recalculation!"<<endl;
			return set_grid;

		} else {
			// go on
		}
	}

	// Y-Richtung
	cout << "update y mit dtody=" << dtody << endl;


	//Berechne U und G
	computation->compute_u_2d(cs, set_grid);
	computation->compute_g_2d(g, set_grid);

	// Berechne Lax-Friedrich Flüsse
	// Faktor 0.25 nach Formel für unsplitting,
	// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting

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

			set_grid = split_grid[0];

			d = set_grid->cellsgrid[pos][0];
			ux = set_grid->cellsgrid[pos][2];
			uy = set_grid->cellsgrid[pos][4];
			uxr = set_grid->cellsgrid[pos][3];
			uyr = set_grid->cellsgrid[pos][5];

			uxd = ux * d;
			uyd = uy * d;

			index_y1 = (y - 1) + (height_m1) * (x);
			index_y2 = (y) + (height_m1) * (x);

			d = d + dtody * (f_lax[1][index_y1 + index_y_end * (0)] - f_lax[1][index_y2 + index_y_end * (0)]);
			uxd = uxd + dtody * (f_lax[1][index_y1 + index_y_end * (1)] - f_lax[1][index_y2 + index_y_end * (1)]);
			uyd = uyd + dtody * (f_lax[1][index_y1 + index_y_end * (2)] - f_lax[1][index_y2 + index_y_end * (2)]);
			uxr = uxr + dtody * (f_lax[1][index_y1 + index_y_end * (3)] - f_lax[1][index_y2 + index_y_end * (3)]);
			uyr = uyr + dtody * (f_lax[1][index_y1 + index_y_end * (4)] - f_lax[1][index_y2 + index_y_end * (4)]);

			// Falls xtoy und ytox gemittelt werden sollen, Werte in einem seperaten Grid setzen um anschließend ins Hauptgrid "zu mitteln"
			if(with_average==0){
				set_grid = grid;
			}

			set_grid->cellsgrid[pos][0] = d;
			set_grid->cellsgrid[pos][1] = constants->ct * pow(d, constants->gamma);
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;

		}
	}
	return set_grid;
}

/**
 *****************************************************************************************
 * Untermethode, die die Berechnung des Lax-Friedrich-Flusses
 * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
 * @param dt Delta t.
 *****************************************************************************************/
Grid* Lax_Friedrich::solve_2d_split_ytox(double dt, int with_average, int rerun) {
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

	// Y-Richtung
	cout << "update y mit dtody=" << dtody << endl;

	set_grid = split_grid[1];

	//Berechne U und G
	computation->compute_u_2d(cs, grid);
	computation->compute_g_2d(g, grid);

	// Berechne Lax-Friedrich Flüsse
	// Faktor 0.25 nach Formel für unsplitting,
	// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting

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

	//TODO: Kürzen?
	if(grid->copy_to(set_grid) == -1){
		cout << "Couldn't copy Grid";
		exit(EXIT_FAILURE);
	}

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
			set_grid->cellsgrid[pos][1] = constants->ct * pow(grid->cellsgrid[pos][0], constants->gamma);
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;

		}
	}

	set_grid->apply_boundary_conditions();

	// ACHTUNG: HIER SOLLTE UEBERPRUEFT WERDEN, OB DER
	// ZEITSCHRITT NICHT ZU GROSS IST FUER DIE 2 RICHTUNG MIT
	// DEN NEUEN WERTEN, SONST KANN EINEM DAS SYSTEM DIVERGIEREN!


	//TODO: Neues Delta T und Eigenwerte überprüfen
	//Neues Delta T darf nicht kleiner sein, als vorher

	if (rerun == 0) {
		time_calculation->cfl_condition(set_grid);
		int dt_comp = time_calculation->compare_dt();

		if (dt_comp == -1) {
			dt = time_calculation->halve_dt(dt);
			with_halved_dt = 1;

			cout << "New delta t is smaller than before!" << endl;
			cout << "Restarting update with dt/2!" << endl;

			set_grid = solve_2d_split_ytox(dt, with_average,1);

			cout << "Finished recalculation!"<<endl;
			return set_grid;

		}
		else {
			// go on
		}
	}

	// X-Richtung
	cout << "update x mit dtodx=" << dtodx << endl;

	//Berechne U, F
	computation->compute_u_2d(cs, set_grid);
	computation->compute_f_2d(f, set_grid);

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

			set_grid = split_grid[1];

			d = set_grid->cellsgrid[pos][0];
			ux = set_grid->cellsgrid[pos][2];
			uy = set_grid->cellsgrid[pos][4];
			uxr = set_grid->cellsgrid[pos][3];
			uyr = set_grid->cellsgrid[pos][5];

			uxd = ux * d;
			uyd = uy * d;

			index_x1 = y + (height_m1) * (x - 1);
			index_x2 = y + (height_m1) * (x);

			d = d + dtodx * (f_lax[0][index_x1 + index_x_end * (0)] - f_lax[0][index_x2 + index_x_end * (0)]);
			uxd = uxd + dtodx * (f_lax[0][index_x1 + index_x_end * (1)] - f_lax[0][index_x2 + index_x_end * (1)]);
			uyd = uyd + dtodx * (f_lax[0][index_x1 + index_x_end * (2)] - f_lax[0][index_x2 + index_x_end * (2)]);
			uxr = uxr + dtodx * (f_lax[0][index_x1 + index_x_end * (3)] - f_lax[0][index_x2 + index_x_end * (3)]);
			uyr = uyr + dtodx * (f_lax[0][index_x1 + index_x_end * (4)] - f_lax[0][index_x2 + index_x_end * (4)]);

			// Falls xtoy und ytox gemittelt werden sollen, Werte in einem seperaten Grid setzen um anschließend ins Hauptgrid "zu mitteln"
			if(with_average==0){
				set_grid = grid;
			}

			set_grid->cellsgrid[pos][0] = d;
			set_grid->cellsgrid[pos][1] = constants->ct * pow(grid->cellsgrid[pos][0], constants->gamma);
			set_grid->cellsgrid[pos][2] = uxd / d;
			set_grid->cellsgrid[pos][4] = uyd / d;
			set_grid->cellsgrid[pos][3] = uxr;
			set_grid->cellsgrid[pos][5] = uyr;

		}
	}
	return set_grid;
}

void Lax_Friedrich::split_mean(Grid* grid_one, Grid* grid_two){
	double d_one, p_one, ux_one, uy_one, uxr_one, uyr_one;
	double d_two, p_two, ux_two, uy_two, uxr_two, uyr_two;
	double d, p, ux, uy, uxr, uyr;

	int pos;
	int order = grid->orderofgrid;

	for (int x = order; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
		for (int y = order; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
			pos = x + y * grid->grid_size_total[0];

			d_one = grid_one->cellsgrid[pos][0];
			p_one = grid_one->cellsgrid[pos][1];
			ux_one = grid_one->cellsgrid[pos][2];
			uy_one = grid_one->cellsgrid[pos][4];
			uxr_one = grid_one->cellsgrid[pos][3];
			uyr_one = grid_one->cellsgrid[pos][5];

			d_two = grid_two->cellsgrid[pos][0];
			p_two = grid_one->cellsgrid[pos][1];
			ux_two = grid_two->cellsgrid[pos][2];
			uy_two = grid_two->cellsgrid[pos][4];
			uxr_two = grid_two->cellsgrid[pos][3];
			uyr_two = grid_two->cellsgrid[pos][5];


			d = 0.5 * (d_one + d_two);
			p = 0.5 * (p_one + p_two);
			ux = 0.5 * (ux_one + ux_two);
			uy = 0.5 * (uy_one + uy_two);
			uxr = 0.5 * (uxr_one + uxr_two);
			uyr = 0.5 * (uyr_one + uyr_two);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][1] = p;
			grid->cellsgrid[pos][2] = ux;
			grid->cellsgrid[pos][4] = uy;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;

		}
	}
}

