#include "LaxFriedrichMethod.h"

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor von der Klasse LaxFriedrichMethod.
 * Ruft einfach den Konstrukter von der geerbten Klasse auf.
 *****************************************************************************************/
LaxFriedrichMethod::LaxFriedrichMethod(Constants *constants,
		Computation *computation, Grid *grid) :
		Solver("Lax-Friedrich", constants, computation, grid) {
	width = grid->getwidth();
	height = grid->getheight();
	neqs = gs->neqs;

	uall = new double[neqs * width * height];
	fall = new double[neqs * width * height];
	gall = new double[neqs * width * height];

	cs = new double**[neqs];
	f = new double**[neqs];
	g = new double**[neqs];

	for (int i = 0; i < neqs; i++) {
		cs[i] = new double*[width];
		f[i] = new double*[width];
		g[i] = new double*[width];
		for (int j = 0; j < width; j++) {
			cs[i][j] = uall + (i * width * height) + (j * height);
			f[i][j] = fall + (i * width * height) + (j * height);
			g[i][j] = gall + (i * width * height) + (j * height);

		}
	}
	fiarray = new double[neqs * (width-1) * (height-1) * dimension];

}

LaxFriedrichMethod::~LaxFriedrichMethod() {
	delete uall;
	delete fall;
	delete gall;
	for (int i = 0; i < neqs; i++) {
		delete cs[i];
		delete f[i];
		delete g[i];
	}
	delete cs;
	delete f;
	delete g;

	delete fiarray;

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
double* LaxFriedrichMethod::calc_method_flux(
		double dt, int dir) {
	cout << "Lax-Friedrich Fluss berechnen..." << endl;

	switch (dimension) {
	// Eine Dimension
	case (1): {
		gs->compute_u_1d(cs, grid, CELLS, ordnung);
		gs->compute_f_1d(f, grid, CELLS, ordnung);

		//Berechne Lax-Friedrich-Fluss
		for (int k = 0; k < gs->neqs; k++) {
			for (int i = 0; i < CELLS[0] + ordnung + 1; i++) {
				fiarray[0 + (dimension*0) + (dimension*(height-1)*i) + (dimension*(height-1)*(width-1)*k)]=
				//fi.at(k).at(i).at(0).at(0) =
						0.5 * (f[k][i][0] + f[k][i + 1][0])
						+ 0.5 * (dx / dt) * (cs[k][i][0] - cs[k][i + 1][0]);
			}
		}

		break;
	}

		// Zwei Dimensionen
	case (2): {

		//Berechne U, F und G
		gs->compute_u_2d(cs, grid, CELLS, ordnung);
		gs->compute_f_2d(f, grid, CELLS, ordnung);
		gs->compute_g_2d(g, grid, CELLS, ordnung);

		// Berechne Lax-Friedrich Flüsse
		// Faktor 0.25 nach Formel für unsplitting,
		// Faktor 0.5 für die Reproduktion der 1-d Ergebnisse und splitting
		//int index = 0;

		for (int k = 0; k < gs->neqs; k++) {
			for (int x = 0; x < CELLS[0] + ordnung + 1; x++) {
				for (int y = 0; y < CELLS[1] + ordnung + 1; y++) {
					fiarray[0 + (dimension*y) + (dimension*(height-1)*x) + (dimension*(height-1)*(width-1)*k)]=
					//fi.at(k).at(x).at(y).at(0) =
							0.5
							* (f[k][x][y] + f[k][x + 1][y])

					//+ 0.25*(dx/dt)*(cs[k][x][y] - cs[k][x+1][y]);
							+ 0.5 * (dx / dt) * (cs[k][x][y] - cs[k][x + 1][y]);

					fiarray[1 + (dimension*y) + (dimension*(height-1)*x) + (dimension*(height-1)*(width-1)*k)]=

					//fi.at(k).at(x).at(y).at(1) =
							0.5
							* (g[k][x][y] + g[k][x][y + 1])
					//+ 0.25*(dy/dt)*(cs[k][x][y] - cs[k][x][y+1]);
							+ 0.5 * (dy / dt) * (cs[k][x][y] - cs[k][x][y + 1]);


					/*
					for (int dim = 0; dim <=1;dim++){
						cout << "Calculated Index: "<< index <<endl;
						//cout << "Index: "<< dim + (dimension*y) + (dimension*(CELLS[1] + ordnung)*x) + (dimension*(CELLS[1] + ordnung)*(CELLS[0] + ordnung)*k) <<endl;
						cout << "Index: "<< dim + (dimension*y) + (dimension*(height-1)*x) + (dimension*(height-1)*(width-1)*k) <<endl;
						//cout << "["<<k<<"]"<< "["<<x<<"]"<< "["<<y<<"]"<< "["<<dim<<"]" <<endl;

						index++;
					}
					*/


				}
			}
		}

		/*cout << "Array Size: "<< neqs * (width-1) * (height-1) * dimension<<endl;
		cout << "Array Size: "<< neqs * width * height * dimension<<endl;
		cout << "Cells0+o "<<CELLS[0] + ordnung<< " Cells1+o "<<CELLS[1] + ordnung<<endl;
		cout << "Width "<<width<< " height "<< height<<endl;*/


		break;
	}
	}

	return fiarray;
}
