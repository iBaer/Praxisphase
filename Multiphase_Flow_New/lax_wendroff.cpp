#include "lax_wendroff.h"

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor von der Klasse Lax_Wendroff.
 * Ruft den Konstruktor der Superklasse Solver auf.
 * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
 * @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird.
 * @param grid Pointer auf Raster-Objekt.
 *****************************************************************************************/
Lax_Wendroff::Lax_Wendroff(Constants *constants, Computation *computation, Grid *grid) :
		Solver("Lax-Wendroff", constants, computation, grid) {

	allocate_cache(grid);
	// temporär
	/*size_total[1] = grid->grid_size_total[1];
	size_m1[1] = grid->grid_size_total[1] - 1;

	uall = new double[neqs * size_total[0] * size_total[1]];
	fall = new double[neqs * size_total[0] * size_total[1]];
	gall = new double[neqs * size_total[0] * size_total[1]];
	f_laxall = new double[neqs * size_total[0] * size_total[1]];
	f_rieall = new double[neqs * size_total[0] * size_total[1]];
	g_laxall = new double[neqs * size_total[0] * size_total[1]];
	g_rieall = new double[neqs * size_total[0] * size_total[1]];

	cs = new double**[neqs];
	fd = new double**[neqs];
	gd = new double**[neqs];
	f_lax_half = new double**[neqs];
	f_rie = new double**[neqs];
	g_lax_half = new double**[neqs];
	g_rie = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[size_total[0]];
		fd[i] = new double*[size_total[0]];
		gd[i] = new double*[size_total[0]];
		f_lax_half[i] = new double*[size_total[0]];
		f_rie[i] = new double*[size_total[0]];
		g_lax_half[i] = new double*[size_total[0]];
		g_rie[i] = new double*[size_total[0]];

		for (int j = 0; j < size_total[0]; j++) {
			cs[i][j] = uall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			fd[i][j] = fall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			gd[i][j] = gall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_lax_half[i][j] = f_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_rie[i][j] = f_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_lax_half[i][j] = g_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_rie[i][j] = g_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
		}
	}

	f_laxwen = new double[neqs * (size_m1[0]) * (size_m1[1])];*/

}

/**
 *****************************************************************************************
 * Destruktor
 *****************************************************************************************/
Lax_Wendroff::~Lax_Wendroff() {
	delete_cache();
}

void Lax_Wendroff::delete_cache(){
	delete[] uall;
	delete[] fall;
	delete[] gall;
	delete[] f_laxall;
	delete[] f_rieall;
	delete[] g_laxall;
	delete[] g_rieall;

	for (int i = 0; i < neqs; i++) {
		delete[] cs[i];
		delete[] fd[i];
		delete[] gd[i];
		delete[] f_lax_half[i];
		delete[] f_rie[i];
		delete[] g_lax_half[i];
		delete[] g_rie[i];
	}
	delete[] cs;
	delete[] fd;
	delete[] gd;
	delete[] f_lax_half;
	delete[] f_rie;
	delete[] g_lax_half;
	delete[] g_rie;

	delete[] f_laxwen;
}

void Lax_Wendroff::allocate_cache(Grid * grid){
	//TODO: Eigene Funktion
	//TODO: previous grid

	size_total[0] = grid->grid_size_total[0];
	size_m1[0] = grid->grid_size_total[0] - 1;
	size_total[1] = grid->grid_size_total[1];
	size_m1[1] = grid->grid_size_total[1] - 1;

	uall = new double[neqs * size_total[0] * size_total[1]];
	fall = new double[neqs * size_total[0] * size_total[1]];
	gall = new double[neqs * size_total[0] * size_total[1]];
	f_laxall = new double[neqs * size_total[0] * size_total[1]];
	f_rieall = new double[neqs * size_total[0] * size_total[1]];
	g_laxall = new double[neqs * size_total[0] * size_total[1]];
	g_rieall = new double[neqs * size_total[0] * size_total[1]];

	cs = new double**[neqs];
	fd = new double**[neqs];
	gd = new double**[neqs];
	f_lax_half = new double**[neqs];
	f_rie = new double**[neqs];
	g_lax_half = new double**[neqs];
	g_rie = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[size_total[0]];
		fd[i] = new double*[size_total[0]];
		gd[i] = new double*[size_total[0]];
		f_lax_half[i] = new double*[size_total[0]];
		f_rie[i] = new double*[size_total[0]];
		g_lax_half[i] = new double*[size_total[0]];
		g_rie[i] = new double*[size_total[0]];

		for (int j = 0; j < size_total[0]; j++) {
			cs[i][j] = uall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			fd[i][j] = fall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			gd[i][j] = gall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_lax_half[i][j] = f_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_rie[i][j] = f_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_lax_half[i][j] = g_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_rie[i][j] = g_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
		}
	}

	f_laxwen = new double[neqs * (size_m1[0]) * (size_m1[1])];
}

/**
 *****************************************************************************************
 * Implementierung der virtuellen Methode zur Berechnung des Flußes.
 * Wählt lediglich das gewünschte Schema aus und delegiert die Berechnung weiter.
 * @param dt Delta t.
 * @param dir Unsplitting = 0, Splitting = 1.
 *****************************************************************************************/
void Lax_Wendroff::calc_method_flux(double dt, Grid * grid) {
	cout << "Lax-Wendroff Fluss berechnen..." << endl;

	if (this->grid!=grid){
		cout << "New grid! Reallocating cache!"<<endl;
		this->grid = grid;
		delete_cache();
		allocate_cache(grid);
	}

	switch (dimension) {
	// Eine Dimension
	case (1): {
		solve_1d(dt);
	}

		// Zwei Dimensionen
	case (2): {

		if (split_method == 0) {
			solve_2d_unsplit(dt);
			break;
		}

		else{
			break;
		}
	}
	}
	time_calculation->set_new_time(dt);

}

/**
 *****************************************************************************************
 * Untermethode, die die Berechnung des Lax-Wendroff-Flusses
 * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
 * @param dt Delta t.
 *****************************************************************************************/
void Lax_Wendroff::solve_1d(double dt){
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double d, ux, uxd, uxr;

	int pos;

	computation->compute_u_1d(cs, grid);
	computation->compute_f_1d(fd, grid);

	Grid* u_lax_half= new Grid(grid->grid_size_total[0]);

	//LaxFriedrich Halbschritt
	for (int k = 0; k < neqs; k++) {
		for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {
			pos = i;
			u_lax_half->cellsgrid[pos][0] = 0.5 * (cs[0][i][0] + cs[0][i+1][0]) + 0.5 * dtodx * (fd[0][i][0] - fd[0][i+1][0]);
			u_lax_half->cellsgrid[pos][2] = (0.5 * (cs[1][i][0] + cs[1][i+1][0]) + 0.5 * dtodx * (fd[1][i][0] - fd[1][i+1][0])) / u_lax_half->cellsgrid[pos][0];
			u_lax_half->cellsgrid[pos][3] = 0.5 * (cs[2][i][0] + cs[2][i+1][0]) + 0.5 * dtodx * (fd[2][i][0] - fd[2][i+1][0]);

			/*u_lax_half->cellsgrid[pos][0] = 0.5 * (cs[0][i+1][0] + cs[0][i][0]) - 0.5 * dtodx * (fd[0][i+1][0] - fd[0][i][0]);
			u_lax_half->cellsgrid[pos][2] = (0.5 * (cs[1][i+1][0] + cs[1][i][0]) - 0.5 * dtodx * (fd[1][i+1][0] - fd[1][i][0])) / u_lax_half->cellsgrid[pos][0];
			u_lax_half->cellsgrid[pos][3] = 0.5 * (cs[2][i+1][0] + cs[2][i][0]) - 0.5 * dtodx * (fd[2][i+1][0] - fd[2][i][0]);*/

			u_lax_half->cellsgrid[pos][1] = constants->ct * pow(u_lax_half->cellsgrid[pos][0], constants->gamma);
		}
	}

	computation->compute_f_1d(fd, u_lax_half);

	//Letzter Halbschritt
	// Updateschritt

	for (int i = order; i < grid->grid_size_total[0] - grid->orderofgrid; i++) {
		d = grid->cellsgrid[i][0];
		ux = grid->cellsgrid[i][2];
		uxd = d * ux;
		uxr = grid->cellsgrid[i][3];

		d = d 		+ dtodx * (fd[0][i-1][0] - fd[0][i][0]);
		uxd = uxd 	+ dtodx * (fd[1][i-1][0] - fd[1][i][0]);
		uxr = uxr	+ dtodx * (fd[2][i-1][0] - fd[2][i][0]);

		/*d = d 		- dtodx * (fd[0][i][0] - fd[0][i-1][0]);
		uxd = uxd 	- dtodx * (fd[1][i][0] - fd[1][i-1][0]);
		uxr = uxr	- dtodx * (fd[2][i][0] - fd[2][i-1][0]);*/

		grid->cellsgrid[i][0] = d;
		grid->cellsgrid[i][1] = constants->ct * pow(grid->cellsgrid[i][0], constants->gamma);
		grid->cellsgrid[i][2] = uxd / d;
		grid->cellsgrid[i][3] = uxr;
	}

	delete u_lax_half;
}

/**
 *****************************************************************************************
 * Untermethode, die die Berechnung des Lax-Wendroff-Flusses
 * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
 * @param dt Delta t.
 *****************************************************************************************/
void Lax_Wendroff::solve_2d_unsplit(double dt) {
	int order = grid->orderofgrid;
	double dtodx = dt / dx ;
	double dtody = dt / dy ;
	double d, ux, uy, uxd, uyd, uxr, uyr;

	computation->compute_u_2d(cs, grid);
	computation->compute_f_2d(fd, grid);
	computation->compute_g_2d(gd, grid);

	Grid* u_lax_xhalf= new Grid(grid->grid_size_total[0], grid->grid_size_total[1], constants);
	Grid* u_lax_yhalf = new Grid(grid->grid_size_total[0], grid->grid_size_total[1], constants);
	Grid* u_lax_allhalf = new Grid(grid->grid_size_total[0], grid->grid_size_total[1], constants);


	int pos = 0;
	//int index = 0;
	//Halber Lax-Friedrich Schritt berechnen in y- und in x-Richtung

	for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid - 0; y++) {
		for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid + 1; x++) {
			pos = x + y * grid->grid_size_total[0];
			x = x - 1;
			u_lax_yhalf->cellsgrid[pos][0] = 	0.5 * (cs[0][x + 1][y + 1] + cs[0][x + 1][y]) + 0.25 * dtody * (gd[0][x + 1][y] - gd[0][x + 1][y + 1]);
			u_lax_yhalf->cellsgrid[pos][2] = 	(0.5 * (cs[1][x + 1][y + 1] + cs[1][x + 1][y]) + 0.25 * dtody * (gd[1][x + 1][y] - gd[1][x + 1][y + 1])) / u_lax_yhalf->cellsgrid[pos][0];
			u_lax_yhalf->cellsgrid[pos][4] = 	(0.5 * (cs[2][x + 1][y + 1] + cs[2][x + 1][y]) + 0.25 * dtody * (gd[2][x + 1][y] - gd[2][x + 1][y + 1])) / u_lax_yhalf->cellsgrid[pos][0];
			u_lax_yhalf->cellsgrid[pos][3] = 	0.5 * (cs[3][x + 1][y + 1] + cs[3][x + 1][y]) + 0.25 * dtody * (gd[3][x + 1][y] - gd[3][x + 1][y + 1]);
			u_lax_yhalf->cellsgrid[pos][5] = 	0.5 * (cs[4][x + 1][y + 1] + cs[4][x + 1][y]) + 0.25 * dtody * (gd[4][x + 1][y] - gd[4][x + 1][y + 1]);

			u_lax_yhalf->cellsgrid[pos][1] = constants->ct * pow(u_lax_yhalf->cellsgrid[pos][0], constants->gamma);

			x = x + 1;

			/*u_lax_yhalf->cellsgrid[pos][0] = 	0.5 * (cs[0][x + 1][y + 1] + cs[0][x + 1][y]) + 0.25 * dtody * (gd[0][x + 1][y + 1] - gd[0][x + 1][y]);
			u_lax_xhalf->cellsgrid[pos][0] = 	0.5 * (cs[0][x + 1][y + 1] + cs[0][x][y + 1]) + 0.25 * dtodx * (fd[0][x + 1][y + 1] - fd[0][x][y + 1]);

			u_lax_yhalf->cellsgrid[pos][2] = 	(0.5 * (cs[1][x + 1][y + 1] + cs[1][x + 1][y]) + 0.25 * dtody * (gd[1][x + 1][y + 1] - gd[1][x + 1][y])) / u_lax_yhalf->cellsgrid[pos][0];
			u_lax_xhalf->cellsgrid[pos][2] = 	(0.5 * (cs[1][x + 1][y + 1] + cs[1][x][y + 1]) + 0.25 * dtodx * (fd[1][x + 1][y + 1] - fd[1][x][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];

			u_lax_yhalf->cellsgrid[pos][4] = 	(0.5 * (cs[2][x + 1][y + 1] + cs[2][x + 1][y]) + 0.25 * dtody * (gd[2][x + 1][y + 1] - gd[2][x + 1][y])) / u_lax_yhalf->cellsgrid[pos][0];
			u_lax_xhalf->cellsgrid[pos][4] = 	(0.5 * (cs[2][x + 1][y + 1] + cs[2][x][y + 1]) + 0.25 * dtodx * (fd[2][x + 1][y + 1] - fd[2][x][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];

			u_lax_yhalf->cellsgrid[pos][3] = 	0.5 * (cs[3][x + 1][y + 1] + cs[3][x + 1][y]) + 0.25 * dtody * (gd[3][x + 1][y + 1] - gd[3][x + 1][y]);
			u_lax_xhalf->cellsgrid[pos][3] = 	0.5 * (cs[3][x + 1][y + 1] + cs[3][x][y + 1]) + 0.25 * dtodx * (fd[3][x + 1][y + 1] - fd[3][x][y + 1]);

			u_lax_yhalf->cellsgrid[pos][5] = 	0.5 * (cs[4][x + 1][y + 1] + cs[4][x + 1][y]) + 0.25 * dtody * (gd[4][x + 1][y + 1] - gd[4][x + 1][y]);
			u_lax_xhalf->cellsgrid[pos][5] = 	0.5 * (cs[4][x + 1][y + 1] + cs[4][x][y + 1]) + 0.25 * dtodx * (fd[4][x + 1][y + 1] - fd[4][x][y + 1]);*/

			// Randwerte in X werden mit Nachbarpunkten gleichgesetzt
			/*if(x == grid->grid_size_total[0] - grid->orderofgrid - 1){
				for(int i=0;i<6;i++){
					u_lax_yhalf->cellsgrid[pos+1][i] = u_lax_yhalf->cellsgrid[pos][i];
					u_lax_xhalf->cellsgrid[pos+1][i] = u_lax_xhalf->cellsgrid[pos][i];
				}
			}
			// Randwerte in Y werden mit Nachbarpunkten gleichgesetzt
			if(y == grid->grid_size_total[1] - grid->orderofgrid - 1){
				for(int i=0;i<6;i++){
					u_lax_yhalf->cellsgrid[pos + grid->grid_size_total[0]][i] = u_lax_yhalf->cellsgrid[pos][i];
					u_lax_xhalf->cellsgrid[pos + grid->grid_size_total[0]][i] = u_lax_xhalf->cellsgrid[pos][i];

				}
			}*/

			/*for(int i=0;i<6;i++){
				if(isnan(u_lax_xhalf->cellsgrid[pos][i])){
					cout << i << "half x: "<< u_lax_xhalf->cellsgrid[pos][i] << " | y "<< u_lax_yhalf->cellsgrid[pos][i] <<endl;
					cout << "ux cs "<< cs[1][x + 1][y + 1] + cs[1][x][y + 1] << " | fd "<< fd[1][x][y + 1] - fd[1][x + 1][y + 1] <<endl;
				}
			}*/

		}
	}

	for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid + 1; y++) {
		for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid - 0; x++) {
			pos = x + y * grid->grid_size_total[0];

			y=y-1;
			u_lax_xhalf->cellsgrid[pos][0] = 	0.5 * (cs[0][x + 1][y + 1] + cs[0][x][y + 1]) + 0.25 * dtodx * (fd[0][x][y + 1] - fd[0][x + 1][y + 1]);
			u_lax_xhalf->cellsgrid[pos][2] = 	(0.5 * (cs[1][x + 1][y + 1] + cs[1][x][y + 1]) + 0.25 * dtodx * (fd[1][x][y + 1] - fd[1][x + 1][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];
			u_lax_xhalf->cellsgrid[pos][4] = 	(0.5 * (cs[2][x + 1][y + 1] + cs[2][x][y + 1]) + 0.25 * dtodx * (fd[2][x][y + 1] - fd[2][x + 1][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];
			u_lax_xhalf->cellsgrid[pos][3] = 	0.5 * (cs[3][x + 1][y + 1] + cs[3][x][y + 1]) + 0.25 * dtodx * (fd[3][x][y + 1] - fd[3][x + 1][y + 1]);
			u_lax_xhalf->cellsgrid[pos][5] = 	0.5 * (cs[4][x + 1][y + 1] + cs[4][x][y + 1]) + 0.25 * dtodx * (fd[4][x][y + 1] - fd[4][x + 1][y + 1]);

			u_lax_xhalf->cellsgrid[pos][1] = constants->ct * pow(u_lax_xhalf->cellsgrid[pos][0], constants->gamma);
			y=y+1;

		}
	}
	// Eck-Randwerte werden mit Nachbarpunkten gleichgesetzt
	/*pos = (grid->grid_size_total[0] - grid->orderofgrid) + (grid->grid_size_total[1] - grid->orderofgrid) * grid->grid_size_total[0];

	for(int i=0;i<6;i++){
		u_lax_yhalf->cellsgrid[pos][i] = u_lax_yhalf->cellsgrid[pos-1][i];
		u_lax_xhalf->cellsgrid[pos][i] = u_lax_xhalf->cellsgrid[pos-1][i];

	}*/

	//cout <<"x half: "<< u_lax_xhalf->cellsgrid[0][0] << " vs "<< u_lax_xhalf->cellsgrid[1][0]<<endl;
	//cout <<"y half: "<<  u_lax_yhalf->cellsgrid[0][0] << " vs "<< u_lax_yhalf->cellsgrid[1][0]<<endl;

	computation->compute_f_2d(fd, u_lax_yhalf);
	computation->compute_g_2d(gd, u_lax_xhalf);

	for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
			for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
				pos = x + y * grid->grid_size_total[0];

			u_lax_allhalf->cellsgrid[pos][0] =
					0.25 * (cs[0][x + 1][y + 1] + cs[0][x + 1][y] + cs[0][x][y + 1] + cs[0][x][y])
					+ 0.5 * dtodx * (fd[0][x][y] - fd[0][x + 1][y]) + 0.5 * dtody * (gd[0][x][y] - gd[0][x][y + 1]);
			u_lax_allhalf->cellsgrid[pos][1] = constants->ct * pow(u_lax_allhalf->cellsgrid[pos][0], constants->gamma);
			u_lax_allhalf->cellsgrid[pos][2] =
					(0.25 * (cs[1][x + 1][y + 1] + cs[1][x + 1][y] + cs[1][x][y + 1] + cs[1][x][y])
					+ 0.5 * dtodx * (fd[1][x][y] - fd[1][x + 1][y]) + 0.5 * dtody * (gd[1][x][y] - gd[1][x][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];
			u_lax_allhalf->cellsgrid[pos][4] =
					(0.25 * (cs[2][x + 1][y + 1] + cs[2][x + 1][y] + cs[2][x][y + 1] + cs[2][x][y])
					+ 0.5 * dtodx * (fd[2][x][y] - fd[2][x + 1][y]) + 0.5 * dtody * (gd[2][x][y] - gd[2][x][y + 1])) / u_lax_xhalf->cellsgrid[pos][0];
			u_lax_allhalf->cellsgrid[pos][3] =
					0.25 * (cs[3][x + 1][y + 1] + cs[3][x + 1][y] + cs[3][x][y + 1] + cs[3][x][y])
					+ 0.5 * dtodx * (fd[3][x][y] - fd[3][x + 1][y]) + 0.5 * dtody * (gd[3][x][y] - gd[3][x][y + 1]);
			u_lax_allhalf->cellsgrid[pos][5] =
					0.25 * (cs[4][x + 1][y + 1] + cs[4][x + 1][y] + cs[4][x][y + 1] + cs[4][x][y])
					+ 0.5 * dtodx * (fd[4][x][y] - fd[4][x + 1][y]) + 0.5 * dtody * (gd[4][x][y] - gd[4][x][y + 1]);

				/*u_lax_allhalf->cellsgrid[pos][0] =
						0.25 * (cs[0][x + 1][y + 1] + cs[0][x + 1][y] + cs[0][x][y + 1] + cs[0][x][y])
						+ 0.5 * dtodx * (fd[0][x][y] - fd[0][x - 1][y]) + 0.5 * dtody * (gd[0][x][y] - gd[0][x][y - 1]);
				u_lax_allhalf->cellsgrid[pos][1] = constants->ct * pow(u_lax_allhalf->cellsgrid[pos][0], constants->gamma);
				u_lax_allhalf->cellsgrid[pos][2] =
						(0.25 * (cs[1][x + 1][y + 1] + cs[1][x + 1][y] + cs[1][x][y + 1] + cs[1][x][y])
						+ 0.5 * dtodx * (fd[1][x][y] - fd[1][x - 1][y]) + 0.5 * dtody * (gd[1][x][y] - gd[1][x][y - 1])) / u_lax_xhalf->cellsgrid[pos][0];
				u_lax_allhalf->cellsgrid[pos][4] =
						(0.25 * (cs[2][x + 1][y + 1] + cs[2][x + 1][y] + cs[2][x][y + 1] + cs[2][x][y])
						+ 0.5 * dtodx * (fd[2][x][y] - fd[2][x - 1][y]) + 0.5 * dtody * (gd[2][x][y] - gd[2][x][y - 1])) / u_lax_xhalf->cellsgrid[pos][0];
				u_lax_allhalf->cellsgrid[pos][3] =
						0.25 * (cs[3][x + 1][y + 1] + cs[3][x + 1][y] + cs[3][x][y + 1] + cs[3][x][y])
						+ 0.5 * dtodx * (fd[3][x][y] - fd[3][x - 1][y]) + 0.5 * dtody * (gd[3][x][y] - gd[3][x][y - 1]);
				u_lax_allhalf->cellsgrid[pos][5] =
						0.25 * (cs[4][x + 1][y + 1] + cs[4][x + 1][y] + cs[4][x][y + 1] + cs[4][x][y])
						+ 0.5 * dtodx * (fd[4][x][y] - fd[4][x - 1][y]) + 0.5 * dtody * (gd[4][x][y] - gd[4][x][y - 1]);*/

			/*
			// Setzt Randwerte in X_max
			if (x == grid->grid_size_total[0] - grid->orderofgrid - 1) {
				for (int i = 0; i < 6; i++) {
					u_lax_allhalf->cellsgrid[pos + 1][i] = u_lax_allhalf->cellsgrid[pos][i];
				}
			}
			// Setzt Randwerte in X_min
			else if (x == grid->orderofgrid) {
				for (int i = 0; i < 6; i++) {
					u_lax_allhalf->cellsgrid[pos - 1][i] = u_lax_allhalf->cellsgrid[pos][i];
				}
			}
			// Setzt Randwerte in Y_max
			if (y == grid->grid_size_total[1] - grid->orderofgrid - 1) {
				for (int i = 0; i < 6; i++) {
					u_lax_allhalf->cellsgrid[pos + grid->grid_size_total[0]][i] = u_lax_allhalf->cellsgrid[pos][i];
				}
			}
			// Setzt Randwerte in Y_min
			else if (y == grid->orderofgrid) {
				for (int i = 0; i < 6; i++) {
					u_lax_allhalf->cellsgrid[pos - grid->grid_size_total[0]][i] = u_lax_allhalf->cellsgrid[pos][i];
				}
			}
			 */
			/*for (int i = 0; i < 6; i++) {

				if (isnan(u_lax_allhalf->cellsgrid[pos][i])) {
					cout << i << " all: " << u_lax_allhalf->cellsgrid[pos][i] << endl;
					cout << "all cs " << cs[0][x + 1][y + 1] + cs[0][x + 1][y] + cs[0][x][y + 1] + cs[0][x][y] << " | all fd "
							<< fd[0][x][y + 1] - fd[0][x + 1][y + 1] << " | all gd " << gd[0][x + 1][y] - gd[0][x + 1][y + 1] << endl;
				}
			}*/

		}
	}

	// Setzen der Werte an 4-Eckpunkten
	/*pos = 0 + (grid->grid_size_total[1] - grid->orderofgrid) * grid->grid_size_total[0];
	for (int i = 0; i < 6; i++) {
		u_lax_allhalf->cellsgrid[pos][i] = u_lax_allhalf->cellsgrid[pos + 1][i];
	}
	pos = (grid->grid_size_total[0] - grid->orderofgrid) + 0 * grid->grid_size_total[0];
	for (int i = 0; i < 6; i++) {
		u_lax_allhalf->cellsgrid[pos][i] = u_lax_allhalf->cellsgrid[pos - 1][i];
	}
	pos = 0;
	for (int i = 0; i < 6; i++) {
		u_lax_allhalf->cellsgrid[pos][i] = u_lax_allhalf->cellsgrid[pos + grid->grid_size_total[0]][i];
	}
	pos = (grid->grid_size_total[0] - grid->orderofgrid) + (grid->grid_size_total[1] - grid->orderofgrid) * grid->grid_size_total[0];
	for (int i = 0; i < 6; i++) {
		u_lax_allhalf->cellsgrid[pos][i] = u_lax_allhalf->cellsgrid[pos - 1][i];
	}*/

	//cout <<"all half: "<<  u_lax_allhalf->cellsgrid[0][0] << " vs "<< u_lax_allhalf->cellsgrid[1][0]<<endl;

	computation->compute_f_2d(fd, u_lax_allhalf);
	computation->compute_g_2d(gd, u_lax_allhalf);

	//cout <<"f: "<<  fd[0][0][0] << " vs "<< fd[0][1][1]<<endl;
	//cout <<"g: "<<  gd[0][0][0] << " vs "<< gd[0][1][1]<<endl;

	//Lax-Wendroff Update Schritt
	cout << "Lax-Wendroff update with dtodx=" << dtodx << " and dtody=" << dtody << endl;
	for (int y = order; y < grid->grid_size_total[1] - grid->orderofgrid+1; y++) {
		for (int x = order; x < grid->grid_size_total[0] - grid->orderofgrid+1; x++) {
			pos = x + y * grid->grid_size_total[0];

			d = grid->cellsgrid[pos][0];
			ux = grid->cellsgrid[pos][2];
			uy = grid->cellsgrid[pos][4];
			uxr = grid->cellsgrid[pos][3];
			uyr = grid->cellsgrid[pos][5];

			uxd = ux * d;
			uyd = uy * d;

			/*if (pos==grid->grid_size_total[0] - grid->orderofgrid+1*grid->grid_size_total[0] || pos==grid->grid_size_total[0] - grid->orderofgrid*2+2*grid->grid_size_total[0]){
				cout <<"d pre: "<< d<<endl;
				cout <<"d f pre: "<< fd[0][x - 1][y] - fd[0][x][y] + fd[0][x - 1][y - 1] - fd[0][x][y - 1]<< "|"<<fd[0][x - 1][y]<< "|"<<-fd[0][x][y]<< "|"<<fd[0][x - 1][y - 1]<< "|"<<-fd[0][x][y - 1]<<endl;
				cout <<"d g pre: "<< gd[0][x][y - 1] - gd[0][x][y] + gd[0][x - 1][y - 1] - gd[0][x - 1][y]<<endl;
			}*/

			d = d     + 0.5 * dtodx * (fd[0][x - 1][y] - fd[0][x][y] + fd[0][x - 1][y - 1] - fd[0][x][y - 1])
					  + 0.5 * dtody * (gd[0][x][y - 1] - gd[0][x][y] + gd[0][x - 1][y - 1] - gd[0][x - 1][y]);
			/*if (pos==1+grid->grid_size_total[0])
							cout <<"d post: "<< d<<endl;*/
			uxd = uxd + 0.5 * dtodx * (fd[1][x - 1][y] - fd[1][x][y] + fd[1][x - 1][y - 1] - fd[1][x][y - 1])
					  + 0.5 * dtody * (gd[1][x][y - 1] - gd[1][x][y] + gd[1][x - 1][y - 1] - gd[1][x - 1][y]);

			uyd = uyd + 0.5 * dtodx * (fd[2][x - 1][y] - fd[2][x][y] + fd[2][x - 1][y - 1] - fd[2][x][y - 1])
					  + 0.5 * dtody * (gd[2][x][y - 1] - gd[2][x][y] + gd[2][x - 1][y - 1] - gd[2][x - 1][y]);

			uxr = uxr + 0.5 * dtodx * (fd[3][x - 1][y] - fd[3][x][y] + fd[3][x - 1][y - 1] - fd[3][x][y - 1])
					  + 0.5 * dtody * (gd[3][x][y - 1] - gd[3][x][y] + gd[3][x - 1][y - 1] - gd[3][x - 1][y]);

			uyr = uyr + 0.5 * dtodx * (fd[4][x - 1][y] - fd[4][x][y] + fd[4][x - 1][y - 1] - fd[4][x][y - 1])
					  + 0.5 * dtody * (gd[4][x][y - 1] - gd[4][x][y] + gd[4][x - 1][y - 1] - gd[4][x - 1][y]);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][1] = constants->ct * pow(grid->cellsgrid[pos][0], constants->gamma);
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;
		}
	}

	delete u_lax_xhalf;
	delete u_lax_yhalf;
	delete u_lax_allhalf;


}
