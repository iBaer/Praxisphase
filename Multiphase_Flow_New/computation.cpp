#include <fstream>

#include "computation.h"
#include "equations.h"

using namespace std;

Computation::Computation(Constants *constants) {
	this->constants = constants;

	cref = constants->cref;
	done = constants->rho_one;
	ccl = constants->ccl;
	gamma = constants->gamma;

	ccl12 = 1.0 - 2.0 * ccl;
	ccl12h = 0.5 * (1.0 - 2.0 * ccl);
	gcinv = 1.0 / gamma;
	gcdgcm1 = gamma / (gamma - 1.0);
	gcm1dgc = (gamma - 1.0) / gamma;
	cclm1 = ccl * (1.0 - ccl);
	powcref = pow(cref, gcinv);

	int dim = constants->dimension;

	if (dim == 1)
		neqs = NEQS_1D;
	if (dim == 2)
		neqs = NEQS_2D;
	//TODO: exprtk oder Header Files müssen neqs setzen

}

Computation& Computation::instance(Constants *constants) {
	static Computation computation(constants);
	return computation;
}

/**
 * Berechnet die Werte U in 1-D
 * @param u, 3-d Feld für die U-Werte
 * @param raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_u_1d(double *** u, Grid * grid) {
	double d, ux, uxr;

	for (int i = 0; i < grid->grid_size_total[0]; i++) {

		d = grid->cellsgrid[i][0];
		ux = grid->cellsgrid[i][2];
		uxr = grid->cellsgrid[i][3];

		u[0][i][0] = COMP_U_1;
		u[1][i][0] = COMP_U_2;
		u[2][i][0] = COMP_U_3;
	}
}

/**
 * Berechnet die Lösungen der Formel für den Fluss F in 1-D
 * @param f, 3-d Feld für die F-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_f_1d(double *** f, Grid * grid) {
	double d, ux, uxr, p;

	for (int i = 0; i < grid->grid_size_total[0]; i++) {

		d = grid->cellsgrid[i][0];
		ux = grid->cellsgrid[i][2];
		uxr = grid->cellsgrid[i][3];
		p = grid->cellsgrid[i][1];

		f[0][i][0] = COMP_F_1;
		f[1][i][0] = COMP_F_2;
		f[2][i][0] = COMP_F_3;
	}
}

/**
 * Berechnet die Werte U in 2-D
 * @param u, 3-d Feld für die U-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_u_2d(double *** u, Grid * grid) {
	int pos;
	double d, ux, uy, uxr, uyr;

	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {
			pos = x + y * grid->grid_size_total[0];
			d = grid->cellsgrid[pos][0];
			ux = grid->cellsgrid[pos][2];
			uy = grid->cellsgrid[pos][4];
			uxr = grid->cellsgrid[pos][3];
			uyr = grid->cellsgrid[pos][5];

			u[0][x][y] = COMP_U_1;
			u[1][x][y] = COMP_U_2;
			u[2][x][y] = COMP_U_4;
			u[3][x][y] = COMP_U_3;
			u[4][x][y] = COMP_U_5;

			for (int i = 0;i<3;i++){
				if(grid->cellsgrid[pos][i]==0){
					//cout <<"U-"<<i<<" out x "<<x<<" | y "<<y<<endl;
				}
				if(isnan(u[i][x][y])){
					cout<<"comp u["<<i<<","<<x<<","<<y<<"] = "<< d<<","<<ux<<","<<uxr<<","<<uy<<","<<uyr<<endl;
				}
			}
		}
	}
}

/**
 * Berechnet die Lösungen der Formel für den Fluss F in 2-D
 * @param f, 3-d Feld für die F-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_f_2d(double *** f, Grid * grid) {
	int pos;
	double d, ux, uy, uxr, uyr, p;

	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {
			pos = x + y * grid->grid_size_total[0];

			d = grid->cellsgrid[pos][0];
			ux = grid->cellsgrid[pos][2];
			uy = grid->cellsgrid[pos][4];
			uxr = grid->cellsgrid[pos][3];
			uyr = grid->cellsgrid[pos][5];
			p = grid->cellsgrid[pos][1];

			/*if(pos == 0 + 0 * grid->grid_size_total[0] || pos == 1 + 1 * grid->grid_size_total[0] || pos == 2 + 2 * grid->grid_size_total[0]){
				cout <<"F d"<<d<<" | ux "<<ux<<endl;
			}*/

			f[0][x][y] = COMP_F_1;
			f[1][x][y] = COMP_F_2;
			f[2][x][y] = COMP_F_4;
			f[3][x][y] = COMP_F_3;
			f[4][x][y] = COMP_F_5;

			for (int i = 0;i<3;i++){
				if(grid->cellsgrid[pos][i]==0){
					//cout <<"F-"<<i<<" out x "<<x<<" | y "<<y<<endl;
				}
				if(isnan(f[i][x][y])){
					cout<<"comp f["<<i<<","<<x<<","<<y<<"] = "<< d<<","<<p<<","<<ux<<","<<uxr<<","<<uy<<","<<uyr<<endl;
				}
			}
		}
	}
}

/**
 * Berechnet die Lösungen der Formel für den Fluss G in 2-D
 * @param g, 3-d Feld für die G-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_g_2d(double *** g, Grid * grid) {
	int pos;
	double d, ux, uy, uxr, uyr, p;

	for (int x = 0; x < grid->grid_size_total[0]; x++) {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {
			pos = x + y * grid->grid_size_total[0];

			d = grid->cellsgrid[pos][0];
			ux = grid->cellsgrid[pos][2];
			uy = grid->cellsgrid[pos][4];
			uxr = grid->cellsgrid[pos][3];
			uyr = grid->cellsgrid[pos][5];
			p = grid->cellsgrid[pos][1];

			g[0][x][y] = COMP_G_1;
			g[1][x][y] = COMP_G_2;
			g[2][x][y] = COMP_G_4;
			g[3][x][y] = COMP_G_3;
			g[4][x][y] = COMP_G_5;

			for (int i = 0;i<3;i++){
				if(grid->cellsgrid[pos][i]==0){
					//cout <<"G-"<<i<<" out x "<<x<<" | y "<<y<<endl;
				}
				if(isnan(g[i][x][y])){
					cout<<"comp g["<<i<<","<<x<<","<<y<<"] = "<< d<<","<<p<<","<<ux<<","<<uxr<<","<<uy<<","<<uyr<<endl;
				}
			}

		}
	}
}
