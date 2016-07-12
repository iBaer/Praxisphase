#include "force.h"

using namespace std;

/**
 * Konstruktor der FORCE Methode.
 * Ruft den Konstruktor der geerbten Klasse auf.
 * @param const_in Dateiname wo Konstanten gespeichert sind.
 * @param formel_in Dateinamen-Kern für die Formeln.
 * @param save_in Dateiname wo für das Laden eine Speicherstands die Plots gespeichert sind.
 */
Force::Force(Constants *constants, Computation *computation, Grid *grid) :
		Solver("FORCE", constants, computation, grid) {

	// temporär
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
	f_lax = new double**[neqs];
	f_rie = new double**[neqs];
	g_lax = new double**[neqs];
	g_rie = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[size_total[0]];
		fd[i] = new double*[size_total[0]];
		gd[i] = new double*[size_total[0]];
		f_lax[i] = new double*[size_total[0]];
		f_rie[i] = new double*[size_total[0]];
		g_lax[i] = new double*[size_total[0]];
		g_rie[i] = new double*[size_total[0]];

		for (int j = 0; j < size_total[0]; j++) {
			cs[i][j] = uall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			fd[i][j] = fall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			gd[i][j] = gall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_lax[i][j] = f_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			f_rie[i][j] = f_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_lax[i][j] = g_laxall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
			g_rie[i][j] = g_rieall + (i * size_total[0] * size_total[1]) + (j * size_total[1]);
		}
	}

	f_force = new double*[dimension];
	for (int i=0;i<dimension;i++){
		f_force[i] = new double[neqs * (size_m1[0]) * (size_m1[1])];
	}

}

Force::~Force() {
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
		delete[] f_lax[i];
		delete[] f_rie[i];
		delete[] g_lax[i];
		delete[] g_rie[i];
	}
	delete[] cs;
	delete[] fd;
	delete[] gd;
	delete[] f_lax;
	delete[] f_rie;
	delete[] g_lax;
	delete[] g_rie;

	for (int i = 0; i < dimension; i++) {
		delete[] f_force[i];
	}
	delete[] f_force;

}

/**
 *****************************************************************************************
 *  Aktualisiert alle zelle mithilfe des berechneten Flusses.
 *****************************************************************************************/

/**
 * Berechnung des FORCE Flusses.
 * @return 4 Dimensionaler Vektor. Zusammenstellung: Gleichung, x-Position, y-Position , dimension
 */
void Force::calc_method_flux(double dt, Grid * grid) {
	cout << "Berechne FORCE Fluss..." << endl;

	//f_force = new double[neqs * (size_m1[0]) * (size_m1[1])];
	this->grid = grid;

	switch (dimension) {
	case (1): {
		solve_1d(dt);
		break;
	}
	case (2): {
		//2D ACHTUNG - EVENTUELL NOCH FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!
		if (split_method == 0) {
			solve_2d_unsplit(dt);
		}
		else{
			solve_2d_split(dt);
		}
		break;

	}
	}
	time_calculation->set_new_time(dt);


}

void Force::solve_1d(double dt){
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double d, ux, uxd, uxr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;

	computation->compute_u_1d(cs, grid);
	computation->compute_f_1d(fd, grid);
	//LaxFriedrichFluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {
			f_lax[k][i][0] = 0.5 * (fd[k][i][0] + fd[k][i + 1][0]) + 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
		}
	}

	//LaxFriedrichFluss berechnen
	/*for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
	 for (unsigned int k = 0; k < gs->u.size(); k++) {
	 f_lax[k][i][0] = 0.5 * (fd[k][i][0] + fd[k][i + 1][0])
	 + 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
	 }
	 }*/

	//Richtmyer Fluss berechnen
	Grid *u_rie = new Grid(size_total[0]);
	for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {
		u_rie->cellsgrid[i][0] = 0.5 * (cs[0][i][0] + cs[0][i + 1][0]) + 0.5 * (dt / dx) * (fd[0][i][0] - fd[0][i + 1][0]);
		u_rie->cellsgrid[i][2] = (0.5 * (cs[1][i][0] + cs[1][i + 1][0]) + 0.5 * (dt / dx) * (fd[1][i][0] - fd[1][i + 1][0])) / u_rie->cellsgrid[i][0];
		u_rie->cellsgrid[i][3] = 0.5 * (cs[2][i][0] + cs[2][i + 1][0]) + 0.5 * (dt / dx) * (fd[2][i][0] - fd[2][i + 1][0]);
		u_rie->cellsgrid[i][1] = constants->ct * pow(u_rie->cellsgrid[i][0], constants->gamma);
	}

	//Richtmyer Fluss berechnen
	computation->compute_f_1d(f_rie, u_rie);

	delete u_rie;

	//FORCE Fluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int i = 0; i < size_total[0] - grid->orderofgrid; i++) {

			f_force[0][0 + (size_m1[1]) * i + (size_m1[1]) * (size_m1[0]) * k] = 0.5 * (f_lax[k][i][0] + f_rie[k][i][0]);
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

		d = d + dtodx * (f_force[0][index_x1 + index_x3 * (0)] - f_force[0][index_x2 + index_x3 * (0)]);
		uxd = uxd + dtodx * (f_force[0][index_x1 + index_x3 * (1)] - f_force[0][index_x2 + index_x3 * (1)]);
		uxr = uxr + dtodx * (f_force[0][index_x1 + index_x3 * (2)] - f_force[0][index_x2 + index_x3 * (2)]);

		grid->cellsgrid[i][0] = d;
		grid->cellsgrid[i][1] = constants->ct * pow(grid->cellsgrid[i][0], constants->gamma);
		grid->cellsgrid[i][2] = uxd / d;
		grid->cellsgrid[i][3] = uxr;
	}
}

void Force::solve_2d_unsplit(double dt){
	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

	computation->compute_u_2d(cs, grid);
	computation->compute_f_2d(fd, grid);
	computation->compute_g_2d(gd, grid);
	//Lax-Friedrich Fluss berechnen
	for (int k = 0; k < computation->neqs; k++) {
		for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
				f_lax[k][x][y] = 0.5 * (fd[k][x][y] + fd[k][x + 1][y]) + 0.25 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);

				g_lax[k][x][y] = 0.5 * (gd[k][x][y] + gd[k][x][y + 1]) + 0.25 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);
			}
		}
	}


	//Richtmyer Fluss berechnen
	Grid* u_rie_f = new Grid(grid->grid_size_total[0], grid->grid_size_total[1], constants);
	Grid* u_rie_g = new Grid(grid->grid_size_total[0], grid->grid_size_total[1], constants);

	for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
		for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
			pos = x + y * size_total[0];

			u_rie_f->cellsgrid[pos][0] = 0.5 * (cs[0][x][y] + cs[0][x + 1][y]) + (dt / 2 * dx) * (fd[0][x][y] - fd[0][x + 1][y]);
			u_rie_f->cellsgrid[pos][2] = (0.5 * (cs[1][x][y] + cs[1][x + 1][y]) + (dt / 2 * dx) * (fd[1][x][y] - fd[1][x + 1][y]))
					/ u_rie_f->cellsgrid[pos][0];
			u_rie_f->cellsgrid[pos][4] = (0.5 * (cs[2][x][y] + cs[2][x + 1][y]) + (dt / 2 * dx) * (fd[2][x][y] - fd[2][x + 1][y]))
					/ u_rie_f->cellsgrid[pos][0];
			u_rie_f->cellsgrid[pos][3] = 0.5 * (cs[3][x][y] + cs[3][x + 1][y]) + (dt / 2 * dx) * (fd[3][x][y] - fd[3][x + 1][y]);
			u_rie_f->cellsgrid[pos][5] = 0.5 * (cs[4][x][y] + cs[4][x + 1][y]) + (dt / 2 * dx) * (fd[4][x][y] - fd[4][x + 1][y]);
			u_rie_f->cellsgrid[pos][1] = constants->ct * pow(u_rie_f->cellsgrid[pos][0], constants->gamma);

			u_rie_g->cellsgrid[pos][0] = 0.5 * (cs[0][x][y] + cs[0][x][y + 1]) + (dt / 2 * dy) * (gd[0][x][y] - gd[0][x][y + 1]);
			u_rie_g->cellsgrid[pos][2] = (0.5 * (cs[1][x][y] + cs[1][x][y + 1]) + (dt / 2 * dy) * (gd[1][x][y] - gd[1][x][y + 1]))
					/ u_rie_g->cellsgrid[pos][0];
			u_rie_g->cellsgrid[pos][4] = (0.5 * (cs[2][x][y] + cs[2][x][y + 1]) + (dt / 2 * dy) * (gd[2][x][y] - gd[2][x][y + 1]))
					/ u_rie_g->cellsgrid[pos][0];
			u_rie_g->cellsgrid[pos][3] = 0.5 * (cs[3][x][y] + cs[3][x][y + 1]) + (dt / 2 * dy) * (gd[3][x][y] - gd[3][x][y + 1]);
			u_rie_g->cellsgrid[pos][5] = 0.5 * (cs[4][x][y] + cs[4][x][y + 1]) + (dt / 2 * dy) * (gd[4][x][y] - gd[4][x][y + 1]);
			u_rie_g->cellsgrid[pos][1] = constants->ct * pow(u_rie_g->cellsgrid[pos][0], constants->gamma);

			/*cout << "d x "<< u_rie_g->cellsgrid[pos][0] << " | y "<< u_rie_g->cellsgrid[pos][0] <<endl;
			cout << "ux x "<< u_rie_g->cellsgrid[pos][2] << " | y "<< u_rie_g->cellsgrid[pos][2] <<endl;
			cout << "ux cs "<< cs[1][x + 1][y + 1] + cs[1][x][y + 1] << " | fd "<< fd[1][x][y + 1] - fd[1][x + 1][y + 1] <<endl;*/
		}
	}


	//Richtmyer Fluss berechnen
	computation->compute_f_2d(f_rie, u_rie_f);
	computation->compute_g_2d(g_rie, u_rie_g);

	//FORCE Fluss berechnen
	for (int k = 0; k < computation->neqs; k++) {
		for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
				f_force[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] = 0.5 * (f_lax[k][x][y] + f_rie[k][x][y]);
				f_force[1][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] = 0.5 * (g_lax[k][x][y] + g_rie[k][x][y]);
			}
		}
	}

	delete u_rie_f;
	delete u_rie_g;

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

			d = d + dtodx * (f_force[0][index_x1 + index_x3 * (0)] - f_force[0][index_x2 + index_x3 * (0)])
					+ dtody * (f_force[1][index_y1 + index_y3 * (0)] - f_force[1][index_y2 + index_y3 * (0)]);
			uxd = uxd + dtodx * (f_force[0][index_x1 + index_x3 * (1)] - f_force[0][index_x2 + index_x3 * (1)])
					+ dtody * (f_force[1][index_y1 + index_y3 * (1)] - f_force[1][index_y2 + index_y3 * (1)]);
			uyd = uyd + dtodx * (f_force[0][index_x1 + index_x3 * (2)] - f_force[0][index_x2 + index_x3 * (2)])
					+ dtody * (f_force[1][index_y1 + index_y3 * (2)] - f_force[1][index_y2 + index_y3 * (2)]);
			uxr = uxr + dtodx * (f_force[0][index_x1 + index_x3 * (3)] - f_force[0][index_x2 + index_x3 * (3)])
					+ dtody * (f_force[1][index_y1 + index_y3 * (3)] - f_force[1][index_y2 + index_y3 * (3)]);
			uyr = uyr + dtodx * (f_force[0][index_x1 + index_x3 * (4)] - f_force[0][index_x2 + index_x3 * (4)])
					+ dtody * (f_force[1][index_y1 + index_y3 * (4)] - f_force[1][index_y2 + index_y3 * (4)]);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][1] = constants->ct * pow(grid->cellsgrid[pos][0], constants->gamma);
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;
		}
	}
}

void Force::solve_2d_split(double dt){

	int order = grid->orderofgrid;
	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;
	int width_m1 = grid->grid_size_total[0] - 1;
	int height_m1 = grid->grid_size_total[1] - 1;
	int pos;

	computation->compute_u_2d(cs, grid);
	computation->compute_f_2d(fd, grid);
	//Lax-Friedrich Fluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
				f_lax[k][x][y] = 0.5 * (fd[k][x][y] + fd[k][x + 1][y]) + 0.25 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);

			}
		}
	}

	//Richtmyer Fluss berechnen
	Grid* u_rie_f = new Grid(size_total[0], size_total[1], constants);

	for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
		for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {

			pos = x + y * size_total[0];

			u_rie_f->cellsgrid[pos][0] = 0.5 * (cs[0][x][y] + cs[0][x + 1][y]) + (dt / 2 * dx) * (fd[0][x][y] - fd[0][x + 1][y]);
			u_rie_f->cellsgrid[pos][2] = (0.5 * (cs[1][x][y] + cs[1][x + 1][y]) + (dt / 2 * dx) * (fd[1][x][y] - fd[1][x + 1][y]))
					/ u_rie_f->cellsgrid[pos][0];
			u_rie_f->cellsgrid[pos][4] = (0.5 * (cs[2][x][y] + cs[2][x + 1][y]) + (dt / 2 * dx) * (fd[2][x][y] - fd[2][x + 1][y]))
					/ u_rie_f->cellsgrid[pos][0];
			u_rie_f->cellsgrid[pos][3] = 0.5 * (cs[3][x][y] + cs[3][x + 1][y]) + (dt / 2 * dx) * (fd[3][x][y] - fd[3][x + 1][y]);
			u_rie_f->cellsgrid[pos][5] = 0.5 * (cs[4][x][y] + cs[4][x + 1][y]) + (dt / 2 * dx) * (fd[4][x][y] - fd[4][x + 1][y]);
			u_rie_f->cellsgrid[pos][1] = constants->ct * pow(u_rie_f->cellsgrid[pos][0], constants->gamma);

		}
	}

	//Richtmyer Fluss berechnen
	computation->compute_f_2d(f_rie, u_rie_f);

	//FORCE Fluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
				f_force[0][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] = 0.5 * (f_lax[k][x][y] + f_rie[k][x][y]);
			}
		}
	}

	delete u_rie_f;

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

			d = d + dtodx * (f_force[0][index_x1 + index_x3 * (0)] - f_force[0][index_x2 + index_x3 * (0)]);
			uxd = uxd + dtodx * (f_force[0][index_x1 + index_x3 * (1)] - f_force[0][index_x2 + index_x3 * (1)]);
			uyd = uyd + dtodx * (f_force[0][index_x1 + index_x3 * (2)] - f_force[0][index_x2 + index_x3 * (2)]);
			uxr = uxr + dtodx * (f_force[0][index_x1 + index_x3 * (3)] - f_force[0][index_x2 + index_x3 * (3)]);
			uyr = uyr + dtodx * (f_force[0][index_x1 + index_x3 * (4)] - f_force[0][index_x2 + index_x3 * (4)]);

			grid->cellsgrid[pos][0] = d;
			//TODO: Muss hier sein?
			//grid->cellsgrid[pos][1] = constants->ct * pow(d, constants->gamma);
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;
		}
	}

	grid->apply_boundary_conditions();

	computation->compute_u_2d(cs, grid);
	computation->compute_g_2d(gd, grid);
	//Lax-Friedrich Fluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
				g_lax[k][x][y] = 0.5 * (gd[k][x][y] + gd[k][x][y + 1]) + 0.25 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);
			}
		}
	}

	//Richtmyer Fluss berechnen
	Grid* u_rie_g = new Grid(size_total[0], size_total[1], constants);

	for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
		for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {
			pos = x + y * size_total[0];

			u_rie_g->cellsgrid[pos][0] = 0.5 * (cs[0][x][y] + cs[0][x][y + 1]) + (dt / 2 * dy) * (gd[0][x][y] - gd[0][x][y + 1]);
			u_rie_g->cellsgrid[pos][2] = (0.5 * (cs[1][x][y] + cs[1][x][y + 1]) + (dt / 2 * dy) * (gd[1][x][y] - gd[1][x][y + 1]))
					/ u_rie_g->cellsgrid[pos][0];
			u_rie_g->cellsgrid[pos][4] = (0.5 * (cs[2][x][y] + cs[2][x][y + 1]) + (dt / 2 * dy) * (gd[2][x][y] - gd[2][x][y + 1]))
					/ u_rie_g->cellsgrid[pos][0];
			u_rie_g->cellsgrid[pos][3] = 0.5 * (cs[3][x][y] + cs[3][x][y + 1]) + (dt / 2 * dy) * (gd[3][x][y] - gd[3][x][y + 1]);
			u_rie_g->cellsgrid[pos][5] = 0.5 * (cs[4][x][y] + cs[4][x][y + 1]) + (dt / 2 * dy) * (gd[4][x][y] - gd[4][x][y + 1]);
			u_rie_g->cellsgrid[pos][1] = constants->ct * pow(u_rie_g->cellsgrid[pos][0], constants->gamma);
		}
	}

	//Richtmyer Fluss berechnen
	computation->compute_g_2d(g_rie, u_rie_g);

	//FORCE Fluss berechnen
	for (int k = 0; k < neqs; k++) {
		for (int x = 0; x < size_total[0] - grid->orderofgrid; x++) {
			for (int y = 0; y < size_total[1] - grid->orderofgrid; y++) {

				f_force[1][y + (size_m1[1]) * x + (size_m1[1]) * (size_m1[0]) * k] = 0.5 * (g_lax[k][x][y] + g_rie[k][x][y]);

			}
		}
	}

	delete u_rie_g;

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


			d = d + dtody * (f_force[1][index_y1 + index_y3 * (0)] - f_force[1][index_y2 + index_y3 * (0)]);
			uxd = uxd + dtody * (f_force[1][index_y1 + index_y3 * (1)] - f_force[1][index_y2 + index_y3 * (1)]);
			uyd = uyd + dtody * (f_force[1][index_y1 + index_y3 * (2)] - f_force[1][index_y2 + index_y3 * (2)]);
			uxr = uxr + dtody * (f_force[1][index_y1 + index_y3 * (3)] - f_force[1][index_y2 + index_y3 * (3)]);
			uyr = uyr + dtody * (f_force[1][index_y1 + index_y3 * (4)] - f_force[1][index_y2 + index_y3 * (4)]);

			grid->cellsgrid[pos][0] = d;
			grid->cellsgrid[pos][1] = constants->ct * pow(d, constants->gamma);
			grid->cellsgrid[pos][2] = uxd / d;
			grid->cellsgrid[pos][4] = uyd / d;
			grid->cellsgrid[pos][3] = uxr;
			grid->cellsgrid[pos][5] = uyr;

		}
	}
}
