#include "grid.h"

#include <math.h>
#include <string>
#include <fstream>

#include "constants.h"

using namespace std;

#define KREIS(x,y)  sqrt((x-radius)*(x-radius)+y*y)
#define KREIS_F    "sqrt((x-radius)*(x-radius)+y*y)"

#define WINKEL_60(x,dx) 2*(x-0.01-dx)
#define WINKEL_60F     "2*(x-0.01-dx)"
#define ALPHA_60 60

#define WINKEL_45(x,dx) xpos-0.1*dx
#define WINKEL_45F     "xpos-0.1*dx"
#define ALPHA_45 45

/**
 *****************************************************************************************
 * Konstruktor für ein leeres eindimensionales Raster.
 *****************************************************************************************/
Grid::Grid(int x) {
	/*** currently unused ***/
	this->constants = 0;
	this->choice = 0;
	orderofgrid = 0;
	boundary_conditions = 0;
	/************************/
	dimension = 1;
	grid_size = new int[dimension];
	grid_size_total = new int[dimension];

	grid_size_total[0] = x;
	//grid_size_total[1] = 1;


	cellsgrid_size = grid_size_total[0];

	cellsgrid = new double*[cellsgrid_size];

	for (int i = 0; i < cellsgrid_size; i++) {

		cellsgrid[i] = new double[4];
		for (int j = 0; j < 4; j++) {
			cellsgrid[i][j] = 0;
		}
	}
}

/**
 *****************************************************************************************
 * Konstruktor für ein leeres zweidimensionales Raster.
 *****************************************************************************************/
Grid::Grid(int x, int y) {
	/*** currently unused ***/
	this->constants = 0;
	this->choice = 0;
	orderofgrid = 0;
	boundary_conditions = 0;
	/************************/

	dimension = 2;
	grid_size = new int[dimension];
	grid_size_total = new int[dimension];

	grid_size_total[0] = x;
	grid_size_total[1] = y;

	cellsgrid_size = grid_size_total[0] * grid_size_total[1];

	cellsgrid = new double*[cellsgrid_size];

	for (int i = 0; i < cellsgrid_size; i++) {

		cellsgrid[i] = new double[6];
		for (int j = 0; j < 6; j++) {
			cellsgrid[i][j] = 0;
		}
	}
}

/**
 *****************************************************************************************
 * Konstruktor des Rasters.
 * @param const_in Pfad zur Datei, welche die initialisierungs Parameter enthält.
 * @param function_in Pfad zur Datei, in der die Initierungsfunktion steht.
 *****************************************************************************************/
Grid::Grid(Constants *constants, string save_in, int choice) {
	this->constants = constants;
	this->choice = choice;
	dimension = (int) constants->dimension;
	grid_size = new int[dimension];
	grid_size_total = new int[dimension];
	boundary_conditions = new int[dimension*2];
	
	boundary_conditions[0] = constants->bc_x_min;
	boundary_conditions[1] = constants->bc_x_max;
	if(dimension>=2){
		boundary_conditions[2] = constants->bc_y_min;
		boundary_conditions[3] = constants->bc_y_max;
		if(dimension>=3){
			boundary_conditions[4] = 0; //TODO: constants->bc_z_min;
			boundary_conditions[5] = 0; //constants->bc_z_max;
		}
	}
	
	orderofgrid = (int) constants->order;

	grid_size[0] = (int) constants->grid_size_x;
	grid_size_total[0] = grid_size[0] + 2 * orderofgrid + 1;

	//grid_size[1] = 1;
	if (dimension == 2) {
		grid_size[1] = (int) constants->grid_size_y;
		grid_size_total[1] = grid_size[1] + 2 * orderofgrid + 1;
	} else if (dimension == 3) {
		grid_size[2] = 1; // TODO: (int) constants->grid_size_z
		grid_size_total[2] = grid_size[2] + 2 * orderofgrid + 1;
	}

	cellsgrid_size = 1;
	for (int i = 0; i < dimension; i++) {
		cellsgrid_size *= grid_size_total[i];
	}

	//choice = 0;
	cellsgrid = new double*[cellsgrid_size];
	// TODO: DEPTH

	// d,p,ux,uxr,uy,uyr,uz,uzr (... ?)
	int number_of_variables = 1;	//TODO: Von neqs(+1) abhängig?
	if (dimension == 1) {
		number_of_variables = 4;
	} else if (dimension == 2) {
		number_of_variables = 6;
	} else if (dimension == 3) {
		number_of_variables = 8;
	}
	for (int i = 0; i < cellsgrid_size; i++) {
		cellsgrid[i] = new double[number_of_variables];
		for (int j = 0; j < number_of_variables; j++) {
			cellsgrid[i][j] = 0;
		}
	}

	switch (dimension) {
	case (1):
		//0 = Riemann-Problem; 1 = 1D Schockwelle;
		switch (choice) {
		case (0): {
			init_1d_rarefraction();
			break;
		}
		case (1): {
			init_1d_shockwave();
			break;
		}
		}
		break;
	case (2):
		//0,1,2 -> 2D Riemann-Probleme
		//0 = Rarefraction wave in x or y direction; 1 = Rarefraction wave in 60 degrees direction; 2 = Rarefraction wave in 45 degrees direction
		//3 = Schockwelle-auf-Blase-Simulation; 4 = Konfiguration laden
		switch (choice) {
		case (0): {
			init_2d_rarefraction_0();
			break;
		}
		case (1): {
			init_2d_rarefraction_60();
			break;
		}
		case (2): {
			init_2d_rarefraction_45();
			break;
		}
		case (3): {
			init_2d_shockwave_bubble();
			break;
		}
		case (4): {
			init_2d_load_file(save_in);
			break;
		}
		}
		break;
	case (3):
		break;
	}
}

/**
 *****************************************************************************************
 * Destruktor
 *****************************************************************************************/
Grid::~Grid() {
	//delete[] zelle;

	for (int i = 0; i < cellsgrid_size; i++) {
		delete[] cellsgrid[i];
	}
	delete[] cellsgrid;
	delete[] grid_size_total;
	delete[] grid_size;
}

/**
 *****************************************************************************************
 * boundary condition
 *****************************************************************************************/
void Grid::bcondi() {

	// TODO: ACHTUNG, FUER 2. ORDNUNG (und höher) NOCHMAL ALLE GRENZEN KONTROLLIEREN

	switch (dimension) {
	case (1): {
		if (boundary_conditions[0] == 0) {
			cellsgrid[0][0] = cellsgrid[orderofgrid][0];
			cellsgrid[0][2] = cellsgrid[orderofgrid][2];
			cellsgrid[0][3] = cellsgrid[orderofgrid][3];
		} else {
			cellsgrid[0][0] = cellsgrid[orderofgrid][0];
			cellsgrid[0][2] = -cellsgrid[orderofgrid][2];
			cellsgrid[0][3] = cellsgrid[orderofgrid][3];
		}

		int index_dst = grid_size_total[0] - orderofgrid;
		int index_src = index_dst - 1;

		if (boundary_conditions[1] == 0) {
			cellsgrid[index_dst][0] = cellsgrid[index_src][0];
			cellsgrid[index_dst][2] = cellsgrid[index_src][2];
			cellsgrid[index_dst][3] = cellsgrid[index_src][3];
		} else {
			cellsgrid[index_dst][0] = cellsgrid[index_src][0];
			cellsgrid[index_dst][2] = -cellsgrid[index_src][2];
			cellsgrid[index_dst][3] = cellsgrid[index_src][3];
		}
		break;
	}
	case (2):
		//2D
		//TODO: For-Loop Performance!?
		for (int n = 0; n < 2; n++) {
			for (int x = 0; x < grid_size_total[0]; x++) {
				for (int y = 0; y < grid_size_total[1]; y++) {
					if (x == 0) {

						int index_dst = 0 + y * grid_size_total[0];
						int index_src = index_dst + orderofgrid;

						if (boundary_conditions[0] == 0) {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						} else {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = -cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = -cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						}

					} else if (x == (grid_size_total[0] - orderofgrid)) {

						int index_dst = grid_size_total[0] - orderofgrid + y * grid_size_total[0];
						int index_src = index_dst - 1;

						if (boundary_conditions[1] == 0) {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						} else {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = -cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = -cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						}
					}
					if (y == 0) {

						int index_dst = x;
						int index_src = x + orderofgrid * grid_size_total[0];

						if (boundary_conditions[2] == 0) {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						} else {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = -cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = -cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						}
					} else if (y == (grid_size_total[1] - orderofgrid)) {

						int index_dst = x + (grid_size_total[1] - orderofgrid) * grid_size_total[0];
						int index_src = x + (grid_size_total[1] - orderofgrid - 1) * grid_size_total[0];

						if (boundary_conditions[3] == 0) {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						} else {
							cellsgrid[index_dst][0] = cellsgrid[index_src][0];
							cellsgrid[index_dst][2] = -cellsgrid[index_src][2];
							cellsgrid[index_dst][3] = cellsgrid[index_src][3];
							cellsgrid[index_dst][4] = -cellsgrid[index_src][4];
							cellsgrid[index_dst][5] = cellsgrid[index_src][5];
						}
					}

				}
			}
		}
		break;
	}

}

/**
 *****************************************************************************************
 * Gitter Initialisierungsmethoden
 *****************************************************************************************/

void Grid::init_1d_rarefraction() {
	double xpos = 0.0;
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;
	double rhol = constants->rhol;
	double vl = constants->vl;
	double vrl = constants->vrl;
	double rhor = constants->rhor;
	double vr = constants->vr;
	double vrr = constants->vrr;

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		xpos = mol + (mor - mol) * ((double) (n - orderofgrid) / grid_size[0]);
		if (xpos <= 0.0) {
			cellsgrid[n][0] = rhol;
			cellsgrid[n][2] = vl;
			cellsgrid[n][3] = vrl;
		} else {
			cellsgrid[n][0] = rhor;
			cellsgrid[n][2] = vr;
			cellsgrid[n][3] = vrr;
		}
	}
	return;
}

void Grid::init_1d_shockwave() {
	double xpos = 0.0;
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		xpos = mol + (mor - mol) * ((double) (n - orderofgrid) / grid_size[0]);
		if (xpos == 0.0) {
			cellsgrid[n][0] = 1323;
			cellsgrid[n][2] = 100000;
			cellsgrid[n][3] = 0;
		} else {
			cellsgrid[n][0] = 1000;
			cellsgrid[n][2] = 0;
			cellsgrid[n][3] = 0;
		}
	}
	return;
}

void Grid::init_2d_rarefraction_0() {
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;
	/*double mur = konstanten->mur; //2D
	 double mul = konstanten->mul; //2D*/

	double xpos = 0.0;
	//double ypos = 0.0;  //2D

	double rhol = constants->rhol;
	double vl = constants->vl;
	double vrl = constants->vrl;
	double vyl = constants->vyl;
	double vyrl = constants->vyrl;
	double rhor = constants->rhor;
	double vr = constants->vr;
	double vrr = constants->vrr;
	double vyr = constants->vyr;
	double vyrr = constants->vyrr;

	cout << "Gerade parallel zur y-Achse" << endl;
	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			xpos = mol + (mor - mol) * (((double) n - orderofgrid) / (double) grid_size[0]);
			/*ypos = mul
			 + (mur - mul)
			 * (((double) m - orderofgrid)
			 / (double) grid_size[1]);*/

			int index = n + m * grid_size_total[0];

			if (xpos <= 0.0) {
				cellsgrid[index][0] = rhol;
				cellsgrid[index][2] = vl;
				cellsgrid[index][3] = vrl;
				cellsgrid[index][4] = vyl;
				cellsgrid[index][5] = vyrl;

			} else {
				cellsgrid[index][0] = rhor;
				cellsgrid[index][2] = vr;
				cellsgrid[index][3] = vrr;
				cellsgrid[index][4] = vyr;
				cellsgrid[index][5] = vyrr;

			}
		}
	}
	return;
}

void Grid::init_2d_rarefraction_90() {
	//double mor = constants->mor;
	//double mol = constants->mol;
	double mur = constants->pos_y_max; //2D
	double mul = constants->pos_y_min; //2D

	//double xpos = 0.0;
	double ypos = 0.0;  //2D

	double rhol = constants->rhol;
	double vl = constants->vl;
	double vrl = constants->vrl;
	double vyl = constants->vyl;
	double vyrl = constants->vyrl;
	double rhor = constants->rhor;
	double vr = constants->vr;
	double vrr = constants->vrr;
	double vyr = constants->vyr;
	double vyrr = constants->vyrr;

	cout << "Gerade parallel zur y-Achse" << endl;
	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			/*xpos = mol
			 + (mor - mol)
			 * (((double) n - orderofgrid)
			 / (double) grid_size[0]);*/
			ypos = mul + (mur - mul) * (((double) m - orderofgrid) / (double) grid_size[1]);

			int index = n + m * grid_size_total[1];
			if (ypos <= 0.0) {
				cellsgrid[index][0] = rhol;
				cellsgrid[index][2] = vl;
				cellsgrid[index][3] = vrl;
				cellsgrid[index][4] = vyl;
				cellsgrid[index][5] = vyrl;

			} else {
				cellsgrid[index][0] = rhor;
				cellsgrid[index][2] = vr;
				cellsgrid[index][3] = vrr;
				cellsgrid[index][4] = vyr;
				cellsgrid[index][5] = vyrr;

			}
		}
	}
	return;
}

void Grid::init_2d_rarefraction_60() {
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;
	double mur = constants->pos_y_max; //2D
	double mul = constants->pos_y_min; //2D

	double dx = (mor - mol) / (double) grid_size[0];

	double xpos = 0.0;
	double ypos = 0.0;  //2D

	double rhol = constants->rhol;
	double vl = constants->vl;
	double vrl = constants->vrl;
	double vyrl = constants->vyrl;
	double rhor = constants->rhor;
	double vr = constants->vr;
	double vrr = constants->vrr;
	double vyrr = constants->vyrr;

	cout << "Gerade mit der Gleichung ypos=" << WINKEL_60F << " als Grenze " << endl;

	double alpha = ALPHA_60 * M_PI / 180.0;
	double vxl_winkel = -fabs(vl) * sin(alpha);
	double vyl_winkel = fabs(vl) * cos(alpha);
	double vxr_winkel = fabs(vr) * sin(alpha);
	double vyr_winkel = -fabs(vr) * cos(alpha);

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			xpos = mol + (mor - mol) * (((double) n - orderofgrid) / (double) grid_size[0]);
			ypos = mul + (mur - mul) * (((double) m - orderofgrid) / (double) grid_size[1]);

			int index = n + m * grid_size_total[0];

			if (ypos >= WINKEL_60(xpos, dx)) {
				cellsgrid[index][0] = rhol;
				cellsgrid[index][2] = vxl_winkel;
				cellsgrid[index][3] = vrl;
				cellsgrid[index][4] = vyl_winkel;
				cellsgrid[index][5] = vyrl;
			} else {
				cellsgrid[index][0] = rhor;
				cellsgrid[index][2] = vxr_winkel;
				cellsgrid[index][3] = vrr;
				cellsgrid[index][4] = vyr_winkel;
				cellsgrid[index][5] = vyrr;
			}
		}
	}
	return;
}

void Grid::init_2d_rarefraction_45() {
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;
	double mur = constants->pos_y_max; //2D
	double mul = constants->pos_y_min; //2D

	double dx = (mor - mol) / (double) grid_size[0];

	double xpos = 0.0;
	double ypos = 0.0;  //2D

	double rhol = constants->rhol;
	double vl = constants->vl;
	double vrl = constants->vrl;
	double vyrl = constants->vyrl;
	double rhor = constants->rhor;
	double vr = constants->vr;
	double vrr = constants->vrr;
	double vyrr = constants->vyrr;

	cout << "Gerade mit der Gleichung ypos=" << WINKEL_45F << " als Grenze " << endl;

	double alpha = ALPHA_45 * M_PI / 180.0;
	double vxl_winkel = -fabs(vl) * sin(alpha);
	double vyl_winkel = fabs(vl) * cos(alpha);
	double vxr_winkel = fabs(vr) * sin(alpha);
	double vyr_winkel = -fabs(vr) * cos(alpha);

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			xpos = mol + (mor - mol) * (((double) n - orderofgrid) / (double) grid_size[0]);
			ypos = mul + (mur - mul) * (((double) m - orderofgrid) / (double) grid_size[1]);

			int index = n + m * grid_size_total[0];
			if (ypos >= WINKEL_45(xpos, dx)) {
				cellsgrid[index][0] = rhol;
				cellsgrid[index][2] = vxl_winkel;
				cellsgrid[index][3] = vrl;
				cellsgrid[index][4] = vyl_winkel;
				cellsgrid[index][5] = vyrl;
			} else {
				cellsgrid[index][0] = rhor;
				cellsgrid[index][2] = vxr_winkel;
				cellsgrid[index][3] = vrr;
				cellsgrid[index][4] = vyr_winkel;
				cellsgrid[index][5] = vyrr;
			}
		}
	}
	return;
}

void Grid::init_2d_shockwave_bubble() {
	double mor = constants->pos_x_max;
	double mol = constants->pos_x_min;
	double mur = constants->pos_y_max; //2D
	double mul = constants->pos_y_min; //2D

	double dx = (mor - mol) / (double) grid_size[0];

	double xpos = 0.0;
	double ypos = 0.0;  //2D

	double radius = constants->radius;

	cout << "Kreis mit Radius " << radius << " und Gleichung " << KREIS_F << endl;

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			xpos = mol + (mor - mol) * (((double) n - orderofgrid) / (double) grid_size[0]);
			ypos = mul + (mur - mul) * (((double) m - orderofgrid) / (double) grid_size[1]);

			int index = n + m * grid_size_total[0];

			if (radius >= KREIS(xpos, ypos)) {
				cellsgrid[index][0] = 8;

			} else {
				cellsgrid[index][0] = 10;

			}

			cellsgrid[index][2] = 0;
			cellsgrid[index][3] = 0;
			cellsgrid[index][4] = 0;
			cellsgrid[index][5] = 0;
		}
	}

	//int shockwave_pos = mor + 1.8/(mor - mol);
	//int shockwave_zell = 1.8/(mor - mol)*cells[0];
	int shockwave_zell = 0.5 * (mor - mol) / (mor - mol) * grid_size[0];
	double shockwave_pos = mol + shockwave_zell * dx;
	cout << "Schockwelle bei " << shockwave_pos << ", x-Zelle " << shockwave_zell << endl;
	for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
		int index = shockwave_zell + m * grid_size_total[0];
		cellsgrid[index][0] = 20;
		cellsgrid[index][2] = 10;

	}
	return;
}

void Grid::init_2d_load_file(std::string save_in) {

	//TODO: Savefile Auswahl
	//Speicherstand laden
	ifstream save_input(save_in.c_str());
	ifstream save_input_line;
	string line_open;
	string line;

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][0] = atof(line.c_str());
		}
	}
	save_input_line.close();

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][2] = atof(line.c_str());

		}
	}
	save_input_line.close();

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][4] = atof(line.c_str());

		}
	}
	save_input_line.close();

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][3] = atof(line.c_str());

		}
	}
	save_input_line.close();

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][5] = atof(line.c_str());

		}
	}
	save_input_line.close();

	getline(save_input, line_open);
	save_input_line.open(line_open.c_str());

	for (int n = orderofgrid; n < grid_size_total[0] - orderofgrid; n++) {
		for (int m = orderofgrid; m < grid_size_total[1] - orderofgrid; m++) {
			getline(save_input_line, line);
			line = line.substr(line.find("\t") + 1);

			cellsgrid[n + m * grid_size_total[0]][1] = atof(line.c_str());

		}
	}
	save_input_line.close();
	return;
}
