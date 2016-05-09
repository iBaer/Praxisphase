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
	width = grid->getwidth();
	height = grid->getheight();
	neqs = gs->neqs;

	uall = new double[neqs * width * height];
	fall = new double[neqs * width * height];
	gall = new double[neqs * width * height];
	f_laxall = new double[neqs * width * height];
	f_rieall = new double[neqs * width * height];
	g_laxall = new double[neqs * width * height];
	g_rieall = new double[neqs * width * height];

	cs = new double**[neqs];
	fd = new double**[neqs];
	gd = new double**[neqs];
	f_lax = new double**[neqs];
	f_rie = new double**[neqs];
	g_lax = new double**[neqs];
	g_rie = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[width];
		fd[i] = new double*[width];
		gd[i] = new double*[width];
		f_lax[i] = new double*[width];
		f_rie[i] = new double*[width];
		g_lax[i] = new double*[width];
		g_rie[i] = new double*[width];

		for (int j = 0; j < width; j++) {
			cs[i][j] = uall + (i * width * height) + (j * height);
			fd[i][j] = fall + (i * width * height) + (j * height);
			gd[i][j] = gall + (i * width * height) + (j * height);
			f_lax[i][j] = f_laxall + (i * width * height) + (j * height);
			f_rie[i][j] = f_rieall + (i * width * height) + (j * height);
			g_lax[i][j] = g_laxall + (i * width * height) + (j * height);
			g_rie[i][j] = g_rieall + (i * width * height) + (j * height);
		}
	}
	fiarray = new double[neqs * (width-1) * (height-1) * dimension];

}

Force::~Force() {
	delete uall;
	delete fall;
	delete gall;
	delete f_laxall;
	delete f_rieall;
	delete g_laxall;
	delete g_rieall;

	for (int i = 0; i < neqs; i++) {
		delete cs[i];
		delete fd[i];
		delete gd[i];
		delete f_lax[i];
		delete f_rie[i];
		delete g_lax[i];
		delete g_rie[i];
	}
	delete cs;
	delete fd;
	delete gd;
	delete f_lax;
	delete f_rie;
	delete g_lax;
	delete g_rie;

	delete fiarray;
}

// ACHTUNG, DIESE FURCHTBAREN VEKTOREN RAUS

/**
 * Berechnung des FORCE Flusses.
 * @return 4 Dimensionaler Vektor. Zusammenstellung: Gleichung, x-Position, y-Position , dimension
 */
double* Force::calc_method_flux(double dt,
		int dir) {
	cout << "Berechne FORCE Fluss..." << endl;

	switch (grid->getdim()) {
	case (1): {
		gs->compute_u_1d(cs, grid, CELLS, ordnung);
		gs->compute_f_1d(fd, grid, CELLS, ordnung);
		cout << "neqs " << gs->neqs << " gs->u.size() " << gs->u.size() << endl;
		//LaxFriedrichFluss berechnen
		for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
			for (int k = 0; k < gs->neqs; k++) {
				f_lax[k][i][0] = 0.5 * (fd[k][i][0] + fd[k][i + 1][0])
						+ 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
			}
		}

		//LaxFriedrichFluss berechnen
		for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
			for (unsigned int k = 0; k < gs->u.size(); k++) {
				f_lax[k][i][0] = 0.5 * (fd[k][i][0] + fd[k][i + 1][0])
						+ 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
			}
		}

		//Richtmyer Fluss berechnen
		/*Grid u_rie(grid->getwidth());
		for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
			u_rie.zelle[i].d = 0.5 * (cs[0][i][0] + cs[0][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[0][i][0] - fd[0][i + 1][0]);
			u_rie.zelle[i].ux = (0.5 * (cs[1][i][0] + cs[1][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[1][i][0] - fd[1][i + 1][0]))
					/ u_rie.zelle[i].d;
			u_rie.zelle[i].uxr = 0.5 * (cs[2][i][0] + cs[2][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[2][i][0] - fd[2][i + 1][0]);
			u_rie.zelle[i].p = konstanten->ct
					* pow(u_rie.zelle[i].d, konstanten->g);
		}*/

		Grid u_rie(grid->getwidth());
		for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
			u_rie.cellsgrid[i][0] = 0.5 * (cs[0][i][0] + cs[0][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[0][i][0] - fd[0][i + 1][0]);
			u_rie.cellsgrid[i][2] = (0.5 * (cs[1][i][0] + cs[1][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[1][i][0] - fd[1][i + 1][0]))
					/ u_rie.cellsgrid[i][0];
			u_rie.cellsgrid[i][3] = 0.5 * (cs[2][i][0] + cs[2][i + 1][0])
					+ 0.5 * (dt / dx) * (fd[2][i][0] - fd[2][i + 1][0]);
			u_rie.cellsgrid[i][1] = konstanten->ct
					* pow(u_rie.cellsgrid[i][0], konstanten->g);
		}

		//Richtmyer Fluss berechnen
		gs->compute_f_1d(f_rie, &u_rie, CELLS, ordnung);

		//FORCE Fluss berechnen
		for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
			for (int k = 0; k < gs->neqs; k++) {

				fiarray[0 + (dimension*0) + (dimension*(height-1)*i) + (dimension*(height-1)*(width-1)*k)]=
				//f_force.at(k).at(i).at(0).at(0) =
						0.5
						* (f_lax[k][i][0] + f_rie[k][i][0]);
			}
		}
		break;
	}
	case (2): {
		//2D ACHTUNG - EVENTUELL NOCH FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!

		gs->compute_u_2d(cs, grid, CELLS, ordnung);
		gs->compute_f_2d(fd, grid, CELLS, ordnung);
		gs->compute_g_2d(gd, grid, CELLS, ordnung);
		cout << "neqs " << gs->neqs << " gs->u.size() " << gs->u.size() << endl;
		//Lax-Friedrich Fluss berechnen
		for (int x = 0; x < CELLS[0] + ordnung + 1; x++) {
			for (int y = 0; y < CELLS[1] + ordnung + 1; y++) {
				for (int k = 0; k < gs->neqs; k++) {
					f_lax[k][x][y] = 0.5 * (fd[k][x][y] + fd[k][x + 1][y])
							+ 0.25 * (dx / dt)
									* (cs[k][x][y] - cs[k][x + 1][y]);

					g_lax[k][x][y] = 0.5 * (gd[k][x][y] + gd[k][x][y + 1])
							+ 0.25 * (dy / dt)
									* (cs[k][x][y] - cs[k][x][y + 1]);
				}
			}
		}

		//Richtmyer Fluss berechnen
		Grid u_rie_f(grid->getwidth(), grid->getheight());
		Grid u_rie_g(grid->getwidth(), grid->getheight());
		int pos;

		for (int x = 0; x < CELLS[0] + ordnung + 1; x++) {
			for (int y = 0; y < CELLS[1] + ordnung + 1; y++) {
				pos = x + y * width;
				u_rie_f.zelle[pos].d = 0.5 * (cs[0][x][y] + cs[0][x + 1][y])
						+ (dt / 2 * dx) * (fd[0][x][y] - fd[0][x + 1][y]);
				u_rie_f.zelle[pos].ux = (0.5 * (cs[1][x][y] + cs[1][x + 1][y])
						+ (dt / 2 * dx) * (fd[1][x][y] - fd[1][x + 1][y]))
						/ u_rie_f.zelle[pos].d;
				u_rie_f.zelle[pos].uy = (0.5 * (cs[2][x][y] + cs[2][x + 1][y])
						+ (dt / 2 * dx) * (fd[2][x][y] - fd[2][x + 1][y]))
						/ u_rie_f.zelle[pos].d;
				u_rie_f.zelle[pos].uxr = 0.5 * (cs[3][x][y] + cs[3][x + 1][y])
						+ (dt / 2 * dx) * (fd[3][x][y] - fd[3][x + 1][y]);
				u_rie_f.zelle[pos].uyr = 0.5 * (cs[4][x][y] + cs[4][x + 1][y])
						+ (dt / 2 * dx) * (fd[4][x][y] - fd[4][x + 1][y]);
				u_rie_f.zelle[pos].p = konstanten->ct
						* pow(u_rie_f.zelle[pos].d, konstanten->g);

				u_rie_g.zelle[pos].d = 0.5 * (cs[0][x][y] + cs[0][x][y + 1])
						+ (dt / 2 * dy) * (gd[0][x][y] - gd[0][x][y + 1]);
				u_rie_g.zelle[pos].ux = (0.5 * (cs[1][x][y] + cs[1][x][y + 1])
						+ (dt / 2 * dy) * (gd[1][x][y] - gd[1][x][y + 1]))
						/ u_rie_g.zelle[pos].d;
				u_rie_g.zelle[pos].uy = (0.5 * (cs[2][x][y] + cs[2][x][y + 1])
						+ (dt / 2 * dy) * (gd[2][x][y] - gd[2][x][y + 1]))
						/ u_rie_g.zelle[pos].d;
				u_rie_g.zelle[pos].uxr = 0.5 * (cs[3][x][y] + cs[3][x][y + 1])
						+ (dt / 2 * dy) * (gd[3][x][y] - gd[3][x][y + 1]);
				u_rie_g.zelle[pos].uyr = 0.5 * (cs[4][x][y] + cs[4][x][y + 1])
						+ (dt / 2 * dy) * (gd[4][x][y] - gd[4][x][y + 1]);
				u_rie_g.zelle[pos].p = konstanten->ct
						* pow(u_rie_g.zelle[pos].d, konstanten->g);
			}
		}

		//Richtmyer Fluss berechnen
		gs->compute_f_2d(f_rie, &u_rie_f, CELLS, ordnung);
		gs->compute_g_2d(g_rie, &u_rie_g, CELLS, ordnung);

		//FORCE Fluss berechnen
		for (int x = 0; x < CELLS[0] + ordnung + 1; x++) {
			for (int y = 0; y < CELLS[1] + ordnung + 1; y++) {
				for (int k = 0; k < gs->neqs; k++) {
					fiarray[0 + (dimension*y) + (dimension*(height-1)*x) + (dimension*(height-1)*(width-1)*k)]=
					//f_force.at(k).at(x).at(y).at(0) =
							0.5
							* (f_lax[k][x][y] + f_rie[k][x][y]);
					fiarray[1 + (dimension*y) + (dimension*(height-1)*x) + (dimension*(height-1)*(width-1)*k)]=

					//f_force.at(k).at(x).at(y).at(1) =
							0.5
							* (g_lax[k][x][y] + g_rie[k][x][y]);
				}
			}
		}

		break;
	}
	}

	return fiarray;
}
