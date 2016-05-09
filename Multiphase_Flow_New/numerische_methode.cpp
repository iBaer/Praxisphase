#include "numerische_methode.h"
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

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor.
 * @param dim Setzt Dimension für die Berechnung.
 * @param ordn Ordnung der numerischen Methode.
 * @param cells Anzahl der Zellen.
 * @param method Name der Methode für unterscheidung bei Output.
 *******************************************1**********************************************/
numerische_methode::numerische_methode(Solver* solver, Constants* constants,
		Computation *computation, Grid *grid) {

	konstanten = constants;
	gs = computation;
	raster = grid;
	name = solver->name;
	this->solver = solver;
	ordnung = 1;

	dimension = konstanten->dimension;
	CELLS = new int(dimension);

	mor = konstanten->mor;
	mol = konstanten->mol;
	mur = konstanten->mur; //2D
	mul = konstanten->mul; //2D
	timeou = konstanten->timeou;
	steps = 0;
	maxnt = konstanten->maxnt;
	teilerend = konstanten->teilerend;
	teiler = konstanten->teiler;
	variante = (int) konstanten->variante;

	CELLS[0] = konstanten->CELLSX;
	CELLS[1] = konstanten->CELLSY;
	g = konstanten->g;
	dx = (mor - mol) / (double) CELLS[0];
	dy = (mur - mul) / (double) CELLS[1];

	dt = 0.0;

	ct = konstanten->ct;

	splitting = 1;
	if (dimension == 2) {
		std::cout << "Wahl!" << endl << " 1: unsplitting, 2: splitting:";
		std::cin >> splitting;
	}
}

/**
 *****************************************************************************************
 * Startet die Berechnung.
 *****************************************************************************************/
void numerische_methode::start_method() {
	double time = 0.0, timetol = 0.000001, timedif = 1.0;
	int step_output = 0;

	cout << "Stepwise Output[1: true :: 0: false]: ";
	cin >> step_output;

	if (dimension == 2) {
		if (splitting == 1)
			cout << "do unsplitting updates" << endl;
		else
			cout << "do splitting updates" << endl;
	}

	if (step_output == 1)
		write();

	for (int n = 1; n <= maxnt && timedif > timetol; n++) {
		cout << n << " : " << maxnt << endl;

		// set boundary conditions
		raster->bcondi(CELLS, ordnung);
		if (step_output == 1)
			write();

		// compute time step
		time = cflcon(n, time);

		// update
		if (splitting == 1)
			// TODO: DER INDEX IN CALC_METHOD_FLUX IST NOCH OHNE BEDEUTUNG,
			// SONDERN ES WIRD IMMER F UND G BERECHNET, MUSS NOCH GEAENDERT WERDEN
			update(solver->calc_method_flux(dt, 1), 1);

		if (splitting == 2) {
			update(solver->calc_method_flux(dt, 1), 1);
			raster->bcondi(CELLS, ordnung);

			// ACHTUNG: HIER SOLLTE UEBERPRUEFT WERDEN, OB DER
			// ZEITSCHRITT NICHT ZU GROSS IST FUER DIE 2 RICHTUNG MIT
			// DEN NEUEN WERTEN, SONST KANN EINEM DAS SYSTEM DIVERGIEREN!
			update(solver->calc_method_flux(dt, 2), 2);
		}
		timedif = fabs(time - timeou);
		steps = n;

	}
	write();
}

/**
 *****************************************************************************************
 * CFL Bedingung anwenden und neue Zeit berechnen.
 * @param n aktueller Zeitschritt.
 * @param time aktuelle Zeit.
 * @return neue Zeit.
 *****************************************************************************************/
double numerische_methode::cflcon(int n, double time) {
	double cref = konstanten->cref;
	double cfl = konstanten->cfl;
	double ccl = konstanten->ccl;
	double done = konstanten->done;
	double gi = 1.0 / g;

	double maxd = 0.0, maxu = 0.0, maxur = 0.0, maxuy = 0.0, maxuyr = 0.0;

	double smax = 0.0, maxs = 0.0;

	int n_eqns;  // number of equations, up in a first line in "formeln....in"

	double uone;
	double utwo;
	double uthree;
	double ufour;
	double ufive;

	if (dimension == 1)
		n_eqns = 3;
	else
		n_eqns = 5;

	string line;

	double p = 0.0, ux = 0.0, d = 0.0, uxr = 0.0, dtwo = 0.0, uy = 0.0, uyr =
			0.0;

	switch (dimension) {
	case (1): {
		//1D
		if ((int) konstanten->calceigv == 1)
		//1 = true, smax wird über Eigenwerte bestimmt
				{
			double values[9];
			LaVectorDouble real(3);
			LaVectorDouble img(3);
			LaGenMatDouble vr(3, 3);

			//Schritt 1: Maxima finden
			for (int x = 0; x < CELLS[0] + 2 * ordnung + 1; x++) {

				// ACHTUNG, NUR VARIANTE 1 IST IN DEN GLEICHUNGEN IMPLEMENTIERT!
				switch (variante) {
				case (1):

					p = ct * pow(raster->cellsgrid[x][0], g);

					raster->set_Zelle_p(p, x);
					dtwo = pow((p / cref), gi);
					break;
					/*case(2):
					 dtwo = ccl/((1/raster->cellsgrid[x][0])-((1-ccl)/done));
					 p = cref * pow(dtwo,g);
					 raster->set_Zelle_p(p,x);
					 break;
					 case(3):
					 dtwo = ccl/((1/raster->cellsgrid[x][0])-((1-ccl)/done));
					 p = ct* pow(raster->cellsgrid[x][0],g);
					 raster->set_Zelle_p(p,x);
					 break;*/
				}

				uone = raster->cellsgrid[x][0];
				utwo = raster->cellsgrid[x][2] * uone;
				uthree = raster->cellsgrid[x][3];

				//Schritt 2: Einsetzen in die Jacobi-Matrix
				matrix_1d(values, n, uone, utwo, uthree, p, done, dtwo, ccl, g,
						ct, cref, variante);

				//Schritt 3: Berechnen der Eigenwerte
				LaGenMatDouble A(values, 3, 3, true);
				LaEigSolve(A, real, img, vr);

				//Schritt 4: Höchsten Eigenwert suchen
				for (int n = 0; n < 3; n++) {
					if (smax < fabs(real(n)))
						smax = fabs(real(n));
				}
			}

			printf("using eig, computed in all cells: %10.4e\n", smax);

			dt = cfl * dx / smax;

			if (n <= teilerend)
				dt = dt * teiler;

			if ((time + dt) > timeou)
				dt = timeou - time;
			time = time + dt;
			cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t"
					<< time << endl;

		}

		else {
			//dt über Näherung berechnen
			for (int i = 0; i < CELLS[0] + 2 * ordnung + 1; i++) {

				d = raster->cellsgrid[i][0];
				uxr = raster->cellsgrid[i][3];
				ux = raster->cellsgrid[i][2];

				switch (variante) {
				case (1): {
					p = ct * pow(d, g);
					raster->set_Zelle_p(p, i);

					dtwo = pow((p / cref), gi);
					break;
				}
				case (2): {
					dtwo = ccl / ((1 / d) - ((1 - ccl) / done));
					p = cref * pow(dtwo, g);
					raster->set_Zelle_p(p, i);

					break;
				}
				case (3): {
					dtwo = ccl / ((1 / d) - ((1 - ccl) / done));
					p = ct * pow(d, g);
					raster->set_Zelle_p(p, i);

					break;
				}
				}

				maxs = fabs(ux)+ AM_1d(g,p,d)+XS_1d(d,dtwo,done,ccl,uxr);

				if (maxs > smax)
					smax = maxs;
			}
			printf("using lambda exact %10.4e\n", smax);

			dt = cfl * dx / smax;

			if (n <= teilerend)
				dt = dt * teiler;

			if ((time + dt) > timeou)
				dt = timeou - time;
			time = time + dt;
			cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t"
					<< time << endl;
		}

		break;

	}
	case (2): {
		//2 Dimensionen

		double smax1 = 0, smax2 = 0;

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, smax wird über Eigenwerte bestimmt

		if ((int) konstanten->calceigv == 0) {
			//über Näherung
			for (int x = 0; x < CELLS[0] + 2 * ordnung + 1; x++) {
				for (int y = 0; y < CELLS[1] + 2 * ordnung + 1; y++) {

					int index = x + y * raster->width;
					d = raster->cellsgrid[index][0];
					uxr = raster->cellsgrid[index][3];
					ux = raster->cellsgrid[index][2];
					uy = raster->cellsgrid[index][4];
					uyr = raster->cellsgrid[index][5];
					p = ct * pow(d, g);
					raster->set_Zelle_p(p, x, y);
					dtwo = pow((p / cref), gi);

					maxs =
							fabs(
									ux)+fabs(uy)+AM_2d(g,p,d)+XS_2d(d,dtwo,done,ccl,uxr,uxy);

					if (maxs > smax) {
						smax = maxs;
					}
				}
			}

			cout << "Größte geratene Eigenwerte: " << smax << endl;

			dt = cfl * min(dx, dy) / smax;

			if (n <= teilerend)
				dt = dt * teiler;

			if ((time + dt) > timeou)
				dt = timeou - time;
			time = time + dt;
			cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t"
					<< time << endl;
		}

		// 0 heisst, 1-d analytische Loseung verwenden,
		// 1 heisst, smax wird über Eigenwerte bestimmt
		else {
			//ueber Eigenwerte

			int n_eqns2 = n_eqns * n_eqns;
			double values1[n_eqns2];
			double values2[n_eqns2];

			/**** Eventuell als Unteraufruf, um dt überprüfen zu können bei Splitting ****/

			//Schritt 1: Daten holen
			for (int x = 0; x < CELLS[0] + 2 * ordnung + 1; x++) {
				for (int y = 0; y < CELLS[1] + 2 * ordnung + 1; y++) {

					maxd = raster->cellsgrid[x][0];
					maxu = raster->cellsgrid[x][2];
					maxuy = raster->cellsgrid[x][4];
					maxur = raster->cellsgrid[x][3];
					maxuyr = raster->cellsgrid[x][5];

					p = ct * pow(raster->cellsgrid[x+y*raster->width][0], g);
					raster->set_Zelle_p(p, x, y);
					dtwo = pow((p / cref), gi);

					maxu = maxu * maxd;
					maxuy = maxuy * maxd;

					//Schritt 2: Einsetzen
					uone = maxd;
					utwo = maxu;
					uthree = maxuy;
					ufour = maxur;
					ufive = maxuyr;
					matrix_2d(values1, values2, n_eqns2, uone, utwo, uthree,
							ufour, ufive, p, done, dtwo, ccl, g, ct, cref);

					//Schritt 3: Berechnen der Eigenwerte
					LaGenMatDouble A(values1, n_eqns, n_eqns, true);
					LaGenMatDouble B(values2, n_eqns, n_eqns, true);
					LaVectorDouble real(n_eqns);
					LaVectorDouble real_b(n_eqns);
					LaVectorDouble img(n_eqns);
					LaVectorDouble img_b(n_eqns);
					LaGenMatDouble vr(n_eqns, n_eqns);
					LaGenMatDouble vr_b(n_eqns, n_eqns);

					LaEigSolve(A, real, img, vr);
					LaEigSolve(B, real_b, img_b, vr_b);

					//Schritt 4: Höchsten Eigenwert suchen
					for (int n = 0; n < n_eqns; n++) {
						if (smax1 < fabs(real(n)))
							smax1 = fabs(real(n));
						if (smax2 < fabs(real_b(n)))
							smax2 = fabs(real_b(n));
					}
				}
			}

			/*****/
			if (splitting == 1) {
				dt = cfl / (smax1 / dx + smax2 / dy);
			} else {
				dt = cfl / max(smax1 / dx, smax2 / dy);
			}

			if (n <= teilerend)
				dt = dt * teiler;

			if ((time + dt) > timeou)
				dt = timeou - time;
			time = time + dt;
			cout << "Größte Eigenwerte: " << smax1 << " und " << smax2 << endl;
			cout << "Neues delta t ist: \t" << dt << " Zeit insgesamt: \t"
					<< time << endl;

		}

		break;
	}

	}
	return time;
}

/**
 *****************************************************************************************
 *  Aktualisiert alle zelle mithilfe des berechneten Flusses.
 *****************************************************************************************/
void numerische_methode::update(double* fi, int dir) {
	cout << "Zellen updaten..." << endl;

	double dtodx = dt / dx;
	double dtody = dt / dy;
	double d, ux, uy, uxd, uyd, uxr, uyr;

	int width = raster->getwidth();
	int height = raster->getheight();

	// update in 1-d
	if (dimension == 1) {
		for (int i = ordnung; i < CELLS[0] + ordnung + 1; i++) {
			d = raster->zelle[i].d;
			ux = raster->zelle[i].ux;
			uxd = d * ux;
			uxr = raster->zelle[i].uxr;

			d = d
					+ dtodx
							* (fi[0 + (dimension * 0)
									+ (dimension * (height - 1) * (i - 1))
									+ (dimension * (height - 1) * (width - 1)
											* 0)]
									- fi[0 + (dimension * 0)
											+ (dimension * (height - 1) * i)
											+ (dimension * (height - 1)
													* (width - 1) * 0)]);
			uxd = uxd
					+ dtodx
							* (fi[0 + (dimension * 0)
									+ (dimension * (height - 1) * (i - 1))
									+ (dimension * (height - 1) * (width - 1)
											* 1)]
									- fi[0 + (dimension * 0)
											+ (dimension * (height - 1) * i)
											+ (dimension * (height - 1)
													* (width - 1) * 1)]);
			uxr = uxr
					+ dtodx
							* (fi[0 + (dimension * 0)
									+ (dimension * (height - 1) * (i - 1))
									+ (dimension * (height - 1) * (width - 1)
											* 2)]
									- fi[0 + (dimension * 0)
											+ (dimension * (height - 1) * i)
											+ (dimension * (height - 1)
													* (width - 1) * 2)]);
			/*d = d + dtodx*(fi.at(0).at(i-1).at(0).at(0)-fi.at(0).at(i).at(0).at(0));
			 uxd = uxd + dtodx*(fi.at(1).at(i-1).at(0).at(0)-fi.at(1).at(i).at(0).at(0));
			 uxr = uxr + dtodx*(fi.at(2).at(i-1).at(0).at(0)-fi.at(2).at(i).at(0).at(0));*/

			raster->zelle[i].d = d;
			raster->zelle[i].ux = uxd / d;
			raster->zelle[i].uxr = uxr;

			raster->cellsgrid[i][0] = d;
			raster->cellsgrid[i][2] = uxd / d;
			raster->cellsgrid[i][3] = uxr;
		}
	}

	// update in 2-d
	if (dimension == 2) {
		if (splitting == 2) {
			if (dir == 1)
				cout << "update x mit dtodx=" << dtodx << endl;
			else
				cout << "update y mit dtody=" << dtody << endl;
		} else {
			cout << "update nonsplitting mit dtodx=" << dtodx << " und dtody="
					<< dtody << endl;
		}
		int pos;

		if (splitting == 1) {
			for (int x = ordnung; x < CELLS[0] + ordnung + 1; x++) {
				for (int y = ordnung; y < CELLS[1] + ordnung + 1; y++) {
					pos = x + y * width;
					d = raster->zelle[pos].d;
					ux = raster->zelle[pos].ux;
					uy = raster->zelle[pos].uy;
					uxr = raster->zelle[pos].uxr;
					uyr = raster->zelle[pos].uyr;
					uxd = ux * d;
					uyd = uy * d;

					d =
							d
									+ dtodx
											* (fi[0 + (dimension * y)
													+ (dimension * (height - 1)
															* (x - 1))
													+ (dimension * (height - 1)
															* (width - 1) * 0)]
													- fi[0 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 0)])
									+ dtody
											* (fi[1 + (dimension * (y - 1))
													+ (dimension * (height - 1)
															* (x))
													+ (dimension * (height - 1)
															* (width - 1) * 0)]
													- fi[1 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 0)]);
					uxd =
							uxd
									+ dtodx
											* (fi[0 + (dimension * y)
													+ (dimension * (height - 1)
															* (x - 1))
													+ (dimension * (height - 1)
															* (width - 1) * 1)]
													- fi[0 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 1)])
									+ dtody
											* (fi[1 + (dimension * (y - 1))
													+ (dimension * (height - 1)
															* (x))
													+ (dimension * (height - 1)
															* (width - 1) * 1)]
													- fi[1 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 1)]);
					uyd =
							uyd
									+ dtodx
											* (fi[0 + (dimension * y)
													+ (dimension * (height - 1)
															* (x - 1))
													+ (dimension * (height - 1)
															* (width - 1) * 2)]
													- fi[0 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 2)])
									+ dtody
											* (fi[1 + (dimension * (y - 1))
													+ (dimension * (height - 1)
															* (x))
													+ (dimension * (height - 1)
															* (width - 1) * 2)]
													- fi[1 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 2)]);
					uxr =
							uxr
									+ dtodx
											* (fi[0 + (dimension * y)
													+ (dimension * (height - 1)
															* (x - 1))
													+ (dimension * (height - 1)
															* (width - 1) * 3)]
													- fi[0 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 3)])
									+ dtody
											* (fi[1 + (dimension * (y - 1))
													+ (dimension * (height - 1)
															* (x))
													+ (dimension * (height - 1)
															* (width - 1) * 3)]
													- fi[1 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 3)]);
					uyr =
							uyr
									+ dtodx
											* (fi[0 + (dimension * y)
													+ (dimension * (height - 1)
															* (x - 1))
													+ (dimension * (height - 1)
															* (width - 1) * 4)]
													- fi[0 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 4)])
									+ dtody
											* (fi[1 + (dimension * (y - 1))
													+ (dimension * (height - 1)
															* (x))
													+ (dimension * (height - 1)
															* (width - 1) * 4)]
													- fi[1 + (dimension * y)
															+ (dimension
																	* (height
																			- 1)
																	* (x))
															+ (dimension
																	* (height
																			- 1)
																	* (width - 1)
																	* 4)]);

					/*d = d + dtodx*(fi.at(0).at(x-1).at(y).at(0)-
					 fi.at(0).at(x).at(y).at(0)) +
					 dtody*(fi.at(0).at(x).at(y-1).at(1)
					 -fi.at(0).at(x).at(y).at(1));
					 uxd =  uxd + dtodx*(fi.at(1).at(x-1).at(y).at(0)-
					 fi.at(1).at(x).at(y).at(0)) +
					 dtody*(fi.at(1).at(x).at(y-1).at(1)-
					 fi.at(1).at(x).at(y).at(1));
					 uyd = uyd + dtodx*(fi.at(2).at(x-1).at(y).at(0)-
					 fi.at(2).at(x).at(y).at(0)) +
					 dtody*(fi.at(2).at(x).at(y-1).at(1)-fi.at(2).at(x).at(y).at(1));
					 uxr = uxr + dtodx*(fi.at(3).at(x-1).at(y).at(0)-
					 fi.at(3).at(x).at(y).at(0)) +
					 dtody*(fi.at(3).at(x).at(y-1).at(1)-fi.at(3).at(x).at(y).at(1));
					 uyr = uyr + dtodx*(fi.at(4).at(x-1).at(y).at(0)-
					 fi.at(4).at(x).at(y).at(0)) +
					 dtody*(fi.at(4).at(x).at(y-1).at(1)-fi.at(4).at(x).at(y).at(1));*/

					/*cout << "d " << d <<endl;
					 cout << "uxd " << uxd<<endl;
					 cout << "uyd " << uyd <<endl;
					 cout << "uxr " << uxr <<endl;
					 cout << "uyr " << uyr <<endl;*/

					raster->zelle[pos].d = d;
					raster->zelle[pos].ux = uxd / d;
					raster->zelle[pos].uy = uyd / d;
					raster->zelle[pos].uxr = uxr;
					raster->zelle[pos].uyr = uyr;

					raster->cellsgrid[pos][0] = d;
					raster->cellsgrid[pos][2] = uxd / d;
					raster->cellsgrid[pos][4] = uyd / d;
					raster->cellsgrid[pos][3] = uxr;
					raster->cellsgrid[pos][5] = uyr;
				}
			}
		}

		if (splitting == 2) {
			if (dir == 1) {
				for (int x = ordnung; x < CELLS[0] + ordnung + 1; x++) {
					for (int y = ordnung; y < CELLS[1] + ordnung + 1; y++) {
						pos = x + y * raster->width;
						d = raster->zelle[pos].d;
						ux = raster->zelle[pos].ux;
						uy = raster->zelle[pos].uy;
						uxr = raster->zelle[pos].uxr;
						uyr = raster->zelle[pos].uyr;

						uxd = ux * d;
						uyd = uy * d;
						/*
						 d = d + dtodx*(
						 fi[0 + (dimension*y) + (dimension*(height-1)*(x-1)) + (dimension*(height-1)*(width-1)*0)]-
						 fi[0 + (dimension*y) + (dimension*(height-1)*(x)) + (dimension*(height-1)*(width-1)*0)])
						 +dtody*(
						 fi[1 + (dimension*(y-1)) + (dimension*(height-1)*(x)) + (dimension*(height-1)*(width-1)*0)]-
						 fi[1 + (dimension*y) + (dimension*(height-1)*(x)) + (dimension*(height-1)*(width-1)*0)]);
						 */

						d = d
								+ dtodx
										* (fi[0 + (dimension * y)
												+ (dimension * (height - 1)
														* (x - 1))
												+ (dimension * (height - 1)
														* (width - 1) * 0)]
												- fi[0 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 0)]);
						uxd = uxd
								+ dtodx
										* (fi[0 + (dimension * y)
												+ (dimension * (height - 1)
														* (x - 1))
												+ (dimension * (height - 1)
														* (width - 1) * 1)]
												- fi[0 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 1)]);
						uyd = uyd
								+ dtodx
										* (fi[0 + (dimension * y)
												+ (dimension * (height - 1)
														* (x - 1))
												+ (dimension * (height - 1)
														* (width - 1) * 2)]
												- fi[0 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 2)]);
						uxr = uxr
								+ dtodx
										* (fi[0 + (dimension * y)
												+ (dimension * (height - 1)
														* (x - 1))
												+ (dimension * (height - 1)
														* (width - 1) * 3)]
												- fi[0 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 3)]);
						uyr = uyr
								+ dtodx
										* (fi[0 + (dimension * y)
												+ (dimension * (height - 1)
														* (x - 1))
												+ (dimension * (height - 1)
														* (width - 1) * 4)]
												- fi[0 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 4)]);

						/*d = d + dtodx*(fi.at(0).at(x-1).at(y).at(0)-
						 fi.at(0).at(x).at(y).at(0));
						 uxd = uxd + dtodx*(fi.at(1).at(x-1).at(y).at(0)-
						 fi.at(1).at(x).at(y).at(0));
						 uyd = uyd + dtodx*(fi.at(2).at(x-1).at(y).at(0)-
						 fi.at(2).at(x).at(y).at(0));
						 uxr = uxr + dtodx*(fi.at(3).at(x-1).at(y).at(0)-
						 fi.at(3).at(x).at(y).at(0));
						 uyr = uyr + dtodx*(fi.at(4).at(x-1).at(y).at(0)-
						 fi.at(4).at(x).at(y).at(0));*/
						raster->zelle[pos].d = d;
						raster->zelle[pos].ux = uxd / d;
						raster->zelle[pos].uy = uyd / d;
						raster->zelle[pos].uxr = uxr;
						raster->zelle[pos].uyr = uyr;

						raster->cellsgrid[pos][0] = d;
						raster->cellsgrid[pos][2] = uxd / d;
						raster->cellsgrid[pos][4] = uyd / d;
						raster->cellsgrid[pos][3] = uxr;
						raster->cellsgrid[pos][5] = uyr;
					}
				}
			} else {
				for (int x = ordnung; x < CELLS[0] + ordnung + 1; x++) {
					for (int y = ordnung; y < CELLS[1] + ordnung + 1; y++) {
						pos = x + y * raster->width;
						d = raster->zelle[pos].d;
						ux = raster->zelle[pos].ux;
						uy = raster->zelle[pos].uy;
						uxr = raster->zelle[pos].uxr;
						uyr = raster->zelle[pos].uyr;

						uxd = ux * d;
						uyd = uy * d;

						/*				dtody*(
						 fi[1 + (dimension*(y-1)) + (dimension*(height-1)*(x)) + (dimension*(height-1)*(width-1)*0)]-
						 fi[1 + (dimension*y) + (dimension*(height-1)*(x)) + (dimension*(height-1)*(width-1)*0)])
						 */

						d = d
								+ dtody
										* (fi[1 + (dimension * (y - 1))
												+ (dimension * (height - 1)
														* (x))
												+ (dimension * (height - 1)
														* (width - 1) * 0)]
												- fi[1 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 0)]);
						uxd = uxd
								+ dtody
										* (fi[1 + (dimension * (y - 1))
												+ (dimension * (height - 1)
														* (x))
												+ (dimension * (height - 1)
														* (width - 1) * 1)]
												- fi[1 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 1)]);
						uyd = uyd
								+ dtody
										* (fi[1 + (dimension * (y - 1))
												+ (dimension * (height - 1)
														* (x))
												+ (dimension * (height - 1)
														* (width - 1) * 2)]
												- fi[1 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 2)]);
						uxr = uxr
								+ dtody
										* (fi[1 + (dimension * (y - 1))
												+ (dimension * (height - 1)
														* (x))
												+ (dimension * (height - 1)
														* (width - 1) * 3)]
												- fi[1 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 3)]);
						uyr = uyr
								+ dtody
										* (fi[1 + (dimension * (y - 1))
												+ (dimension * (height - 1)
														* (x))
												+ (dimension * (height - 1)
														* (width - 1) * 4)]
												- fi[1 + (dimension * y)
														+ (dimension
																* (height - 1)
																* (x))
														+ (dimension
																* (height - 1)
																* (width - 1)
																* 4)]);

						/*d = d + dtody*(fi.at(0).at(x).at(y-1).at(1)-
						 fi.at(0).at(x).at(y).at(1));
						 uxd = uxd + dtody*(fi.at(1).at(x).at(y-1).at(1)-
						 fi.at(1).at(x).at(y).at(1));
						 uyd = uyd + dtody*(fi.at(2).at(x).at(y-1).at(1)-
						 fi.at(2).at(x).at(y).at(1));
						 uxr = uxr + dtody*(fi.at(3).at(x).at(y-1).at(1)-
						 fi.at(3).at(x).at(y).at(1));
						 uyr = uyr + dtody*(fi.at(4).at(x).at(y-1).at(1)-
						 fi.at(4).at(x).at(y).at(1));*/

						raster->zelle[pos].d = d;
						raster->zelle[pos].ux = uxd / d;
						raster->zelle[pos].uy = uyd / d;
						raster->zelle[pos].uxr = uxr;
						raster->zelle[pos].uyr = uyr;

						raster->cellsgrid[pos][0] = d;
						raster->cellsgrid[pos][2] = uxd / d;
						raster->cellsgrid[pos][4] = uyd / d;
						raster->cellsgrid[pos][3] = uxr;
						raster->cellsgrid[pos][5] = uyr;
					}
				}
			}
		}
	}
}

/**
 *****************************************************************************************
 *  Schreibt Ergebnisse in Dateien u,d,ur,p
 *****************************************************************************************/
void numerische_methode::write() {
	double xpos = 0.0;
	double ypos = 0.0;
	double p = 0.0;
	string added;

	if (dimension == 1) {
		added = to_string(CELLS[0]) + "_" + name + "_" + to_string(dimension)
				+ "d_" + to_string(variante) + ".variant_" + "div"
				+ to_string(teiler) + "till" + to_string((int) teilerend) + "_"
				+ to_string(steps) + "Steps";

	} else {
		added = to_string(CELLS[0]) + "x" + to_string(CELLS[1]) + "_" + name
				+ "_" + to_string(dimension) + "d_split" + to_string(splitting)
				+ "_IC" + to_string(raster->choice) + "_div"
				+ to_string(int(1. / teiler)) + "till"
				+ to_string((int) teilerend) + "_" + to_string(steps) + "Steps";
	}

	switch (raster->getdim()) {
	case (1): {
		string d_path = "d_" + added;
		string uxr_path = "uxr_" + added;
		string ux_path = "ux_" + added;
		string p_path = "p_" + added;

		ofstream d_out(d_path.c_str());
		ofstream uxr_out(uxr_path.c_str());
		ofstream ux_out(ux_path.c_str());
		ofstream p_out(p_path.c_str());

		for (int i = ordnung; i < CELLS[0] + ordnung + 1; i++) {
			xpos = (mol + (mor - mol) * ((double) (i - ordnung) / CELLS[0]));

			p = ct * pow(raster->cellsgrid[i][0], g);


			d_out << fixed << setprecision(8) << xpos << " \t"
					<< setprecision(10) << raster->cellsgrid[i][0] << "\n";
			uxr_out << fixed << setprecision(8) << xpos << " \t"
					<< setprecision(10) << raster->cellsgrid[i][3] << "\n";
			ux_out << fixed << setprecision(8) << xpos << " \t"
					<< setprecision(10) << raster->cellsgrid[i][2] << "\n";
			p_out << fixed << setprecision(8) << xpos << " \t" << scientific
					<< setprecision(10) << p << "\n";

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
		for (int x = ordnung; x < CELLS[0] + ordnung + 1; x++) {
			for (int y = ordnung; y < CELLS[1] + ordnung + 1; y++) {
				xpos = mol
						+ (mor - mol)
								* (((double) x - ordnung) / (double) CELLS[0]);
				ypos = mul
						+ (mur - mul)
								* (((double) y - ordnung) / (double) CELLS[1]);

				int index = x + y * raster->width;
				p = ct * pow(raster->cellsgrid[index][0], g);

				d_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << setprecision(10)
						<< raster->cellsgrid[index][0] << "\n";
				uxr_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << setprecision(10)
						<< raster->cellsgrid[index][3] << "\n";
				ux_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << setprecision(10)
						<< raster->cellsgrid[index][2] << "\n";
				uyr_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << setprecision(10)
						<< raster->cellsgrid[index][5] << "\n";
				uy_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << setprecision(10)
						<< raster->cellsgrid[index][4] << "\n";
				p_out << fixed << setprecision(8) << xpos << " " << fixed
						<< setprecision(8) << ypos << " \t" << scientific
						<< setprecision(10) << p << "\n";
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

		double d = 0.0, xout, uout, urout, ux, uy, sign;

		for (int x = ordnung; x < CELLS[0] + ordnung + 1; x++) {
			for (int y = ordnung; y < CELLS[1] + ordnung + 1; y++) {
				if (raster->choice < 3) {
					xpos = mol
							+ (mor - mol)
									* (((double) x - ordnung)
											/ (double) CELLS[0]);
					ypos = mul
							+ (mur - mul)
									* (((double) y - ordnung)
											/ (double) CELLS[1]);

					int index = x + y * raster->width;

					p = ct * pow(raster->cellsgrid[index][0], g);

					if (raster->choice == 0) {
						d = fabs(ypos);
					}
					if (raster->choice == 1) {
						d = fabs(xpos + ypos);
					}
					if (raster->choice == 2) {
						d = fabs(xpos + 2 * ypos);
					}

					if (d < 0.0001) {
						sign = xpos < 0 ? -1.0 : 1.0;
						xout = sign * sqrt(xpos * xpos + ypos * ypos);
						ux = raster->cellsgrid[index][2];
						uy = raster->cellsgrid[index][4];
						sign = ux < 0 ? -1.0 : 1.0;
						uout = sign * sqrt(ux * ux + uy * uy);
						ux = raster->cellsgrid[index][3];
						uy = raster->cellsgrid[index][5];
						sign = ux < 0 ? -1.0 : 1.0;
						urout = sign * sqrt(ux * ux + uy * uy);

						d_out_d1 << fixed << setprecision(8) << xout << " \t"
								<< setprecision(10) << raster->cellsgrid[index][0]
								<< "\n";
						u_out_d1 << fixed << setprecision(8) << xout << " \t"
								<< setprecision(10) << uout << "\n";
						ur_out_d1 << fixed << setprecision(8) << xout << " \t"
								<< setprecision(10) << urout << "\n";
						p_out_d1 << fixed << setprecision(8) << xout << " \t"
								<< scientific << setprecision(10) << p << "\n";
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

/**
 *****************************************************************************************
 *  Jacobi-Matrix für den 1-d Fall
 *****************************************************************************************/

void numerische_methode::matrix_1d(double * values, int n, double uone,
		double utwo, double uthree, double p, double done, double dtwo,
		double ccl, double g, double ct, double cref, int variante) {

	switch (variante) {
	case (1):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(utwo * utwo) / (uone * uone)
				+ ccl * (1 - ccl) * uthree * uthree
				+ g * ct * pow(uone, g - 1.0);
		values[4] = (2 * utwo) / uone;
		values[5] = uone * ccl * (1.0 - ccl) * 2.0 * uthree;
		values[6] = -utwo / (uone * uone) * uthree
				+ pow(cref / ct, 1.0 / g) * g * ct * pow(uone, g - 2.0)
				- (g * ct * pow(uone, g - 1.0)) / done;
		values[7] = uthree / uone;
		values[8] = (utwo / uone) + (1.0 - 2.0 * ccl) * uthree;
		break;
	case (2):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(utwo * utwo) / (uone * uone)
				+ ccl * (1 - ccl) * uthree * uthree
				+ cref * g * pow(done, g - 1.0) * (ccl * done * done)
						/ pow(done + (ccl - 1.0) * uone, 2.0);
		values[4] = (2 * utwo) / uone;
		values[5] = uone * ccl * (1.0 - ccl) * 2.0 * uthree;
		values[6] = -utwo / (uone * uone) * uthree
				+ cref
						* (g * pow(dtwo, g - 2.0)
								- (g * pow(dtwo, g - 1.0)) / done)
						* (ccl * done * done)
						/ pow(done + (ccl - 1.0) * uone, 2.0);
		values[7] = uthree / uone;
		values[8] = (utwo / uone) + (1.0 - 2.0 * ccl) * uthree;
		break;
	case (3):
		values[0] = 0;
		values[1] = 1;
		values[2] = 0;
		values[3] = -(utwo * utwo) / (uone * uone)
				+ ccl * (1 - ccl) * uthree * uthree
				+ g * ct * pow(uone, g - 1.0);
		values[4] = (2.0 * utwo) / uone;
		values[5] = uone * ccl * (1.0 - ccl) * 2.0 * uthree;
		values[6] = -utwo / (uone * uone) * uthree
				+ (1.0 / ccl) * ((1.0 / uone) - (1.0 / done)) * ct * g
						* pow(uone, g - 1.0);
		values[7] = uthree / uone;
		values[8] = (utwo / uone) + (1.0 - 2.0 * ccl) * uthree;
	}
}

/**
 *****************************************************************************************
 *  Jacobi-Matrix für den 2-d Fall in x- und y-Richtung
 *****************************************************************************************/
void numerische_methode::matrix_2d(double * values_x, double * values_y, int n,
		double uone, double utwo, double uthree, double ufour, double ufive,
		double p, double done, double dtwo, double ccl, double g, double ct,
		double cref) {
	values_x[0] = 0;
	values_x[1] = 1;
	values_x[2] = 0;
	values_x[3] = 0;
	values_x[4] = 0;
	values_x[5] = -(utwo * utwo) / (uone * uone)
			+ ccl * (1 - ccl) * ufour * ufour + g * ct * pow(uone, g - 1.0);
	values_x[6] = (2 * utwo) / uone;
	values_x[7] = 0;
	values_x[8] = ccl * (1 - ccl) * 2 * uone * ufour;
	values_x[9] = 0;
	values_x[10] = -(utwo * uthree) / (uone * uone)
			+ ccl * (1 - ccl) * ufour * ufive;
	values_x[11] = uthree / uone;
	values_x[12] = utwo / uone;
	values_x[13] = ccl * (1 - ccl) * uone * ufive;
	values_x[14] = ccl * (1 - ccl) * uone * ufour;
	values_x[15] = -(utwo * ufour) / (uone * uone)
			+ pow(cref / ct, 1.0 / g) * g * ct * pow(uone, g - 2.0)
			- (g * ct * pow(uone, g - 1.0)) / (done);
	values_x[16] = ufour / uone;
	values_x[17] = 0;
	values_x[18] = utwo / uone + (1 - 2 * ccl) * ufour;
	values_x[19] = 0;
	values_x[20] = -(utwo * ufive) / (uone * uone)
			- uthree * ufour / (uone * uone);
	values_x[21] = ufive / uone;
	values_x[22] = ufour / uone;
	values_x[23] = uthree / uone + (1 - 2 * ccl) * ufive;
	values_x[24] = utwo / uone + (1 - 2 * ccl) * ufour;

	values_y[0] = 0;
	values_y[1] = 0;
	values_y[2] = 1;
	values_y[3] = 0;
	values_y[4] = 0;
	values_y[5] = -(utwo * uthree) / (uone * uone)
			+ ccl * (1 - ccl) * ufour * ufive;
	values_y[6] = uthree / uone;
	values_y[7] = utwo / uone;
	values_y[8] = ccl * (1 - ccl) * uone * ufive;
	values_y[9] = ccl * (1 - ccl) * uone * ufour;
	values_y[10] = -(uthree * uthree) / (uone * uone)
			+ ccl * (1 - ccl) * ufive * ufive + g * ct * pow(uone, g - 1.0);
	values_y[11] = 0;
	values_y[12] = (2 * uthree) / uone;
	values_y[13] = 0;
	values_y[14] = ccl * (1 - ccl) * 2 * uone * ufive;
	values_y[15] = -(utwo * ufive) / (uone * uone)
			- uthree * ufour / (uone * uone);
	values_y[16] = ufive / uone;
	values_y[17] = ufour / uone;
	values_y[18] = uthree / uone + (1 - 2 * ccl) * ufive;
	values_y[19] = utwo / uone + (1 - 2 * ccl) * ufour;
	values_y[20] = -(uthree * ufive) / (uone * uone)
			+ pow(cref / ct, 1.0 / g) * g * ct * pow(uone, g - 2.0)
			- (g * ct * pow(uone, g - 1.0)) / (done);
	values_y[21] = 0;
	values_y[22] = ufive / uone;
	values_y[23] = 0;
	values_y[24] = uthree / uone + (1 - 2 * ccl) * ufive;
}
