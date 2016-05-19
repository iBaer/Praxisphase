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

	neqs = computation->neqs;

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
 * Berechnung des FORCE Flusses.
 * @return 4 Dimensionaler Vektor. Zusammenstellung: Gleichung, x-Position, y-Position , dimension
 */
double* Force::calc_method_flux(double dt, int dir) {
	cout << "Berechne FORCE Fluss..." << endl;

	//f_force = new double[neqs * (size_m1[0]) * (size_m1[1])];


	switch (dimension) {
	case (1): {
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
		//break;
		return f_force[0];
	}
	case (2): {
		//2D ACHTUNG - EVENTUELL NOCH FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!
		/*if (dir == 0) {
			computation->compute_u_2d(cs, grid);
			computation->compute_f_2d(fd, grid);
			computation->compute_g_2d(gd, grid);
			cout << "neqs " << computation->neqs << " gs->u.size() " << computation->u.size() << endl;
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
			Grid* u_rie_f = new Grid(grid->grid_size_total[0], grid->grid_size_total[1]);
			Grid* u_rie_g = new Grid(grid->grid_size_total[0], grid->grid_size_total[1]);
			int pos;

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
				}
			}

			//Richtmyer Fluss berechnen
			computation->compute_f_2d(f_rie, u_rie_f);
			computation->compute_g_2d(g_rie, u_rie_g);

			//FORCE Fluss berechnen
			for (int k = 0; k < computation->neqs; k++) {
				for (int x = 0; x < grid->grid_size_total[0] - grid->orderofgrid; x++) {
					for (int y = 0; y < grid->grid_size_total[1] - grid->orderofgrid; y++) {
						fiarray[0 + (dimension * y) + (dimension * (size_m1[1]) * x) + (dimension * (size_m1[1]) * (size_m1[0]) * k)] = 0.5
								* (f_lax[k][x][y] + f_rie[k][x][y]);
						fiarray[1 + (dimension * y) + (dimension * (size_m1[1]) * x) + (dimension * (size_m1[1]) * (size_m1[0]) * k)] = 0.5
								* (g_lax[k][x][y] + g_rie[k][x][y]);
					}
				}
			}

			delete u_rie_f;
			delete u_rie_g;
		}*/
		if (dir == 1) {
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
			Grid* u_rie_f = new Grid(size_total[0], size_total[1]);

			int pos;
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
			return f_force[0];
		} else {
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
			Grid* u_rie_g = new Grid(size_total[0], size_total[1]);
			int pos;

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
			return f_force[1];
		}
		break;

	}
	}

	//return fiarray;
	return f_force[0];

}
