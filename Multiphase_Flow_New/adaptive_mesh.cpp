/*
 * adaptive_mesh.cpp
 *
 *  Created on: 20.06.2016
 *      Author: pascal
 */

#include "adaptive_mesh.h"

#include "grid.h"
#include "solver.h"
#include "constants.h"
#include "cluster_builder.h"

#include <sstream>
#include <set>

//TODO: Datenstruktur ändern?!
struct {
  double** half_flux;
  int dimension;
  int size;
} Coarse_Half, Fine_Half1, Fine_Half2;

Adaptive_Mesh::Adaptive_Mesh(Solver * solver, Grid *grid, Constants *constants, Time_Step_Calculation* time_calculation) {
	this->solver = solver;
	this->grid_main = grid;
	this->constants = constants;

	this->grid_twostep = new Grid(grid_main->grid_size_total[0],grid_main->grid_size_total[1], constants);
	this->grid_doublestep = new Grid(grid_main->grid_size_total[0],grid_main->grid_size_total[1], constants);

	this->time_calculation = time_calculation;
	this->cluster_builder = new Cluster_Builder();

	marked_cells = nullptr;

	dx = (constants->pos_x_max - constants->pos_x_min) / (double) grid->grid_size[0];
	dy = (constants->pos_y_max - constants->pos_y_min) / (double) grid->grid_size[1];

}

Adaptive_Mesh::~Adaptive_Mesh() {
	delete grid_twostep;
	delete grid_doublestep;

}

//TODO: Fancy Output with std::string(level, '-') -> level=global!?
void Adaptive_Mesh::amr() {
	//TODO: Über Koordinaten statt "Zellnummer"
	double time = 0;
	for (int n = 1; n <= 50; n++) {
		cout << "============ " << n << " : " << constants->maxnt << " ============" << endl;
		cout << "============ Time: " << time << " ============" << endl;
		double dt = 0;
		int grid_level = 0;
		Grid * coarse_grid = grid_main;

		// 1. Gröbstes Netz berechnen + Zeitschritt + Doppelschritt
		calculate_grids(coarse_grid, grid_level, dt);
		time += dt;

		// 2. Zellen markieren und prüfen ob Feines Gitter generiert werden muss
		marked_cells = new int[coarse_grid->cellsgrid_size]();	//TODO: Nicht global!? //TODO: Rekursiv anpassen

		int do_refinement = 0;

		double tolerance = 0.1;	// TODO: Toleranz eventuell als Parameter Input Konstante?
		do_refinement = grid_marker(marked_cells,grid_twostep, grid_doublestep, tolerance);

		if (do_refinement) {
			cout << "Adaptive Grid needed! Cluster Generation started!" << endl;
			// 3. "Feineres Gitter erzeugen"
			// 3.1. Cluster erzeugen
			vector<Cluster_Square> clusters;
			cluster_builder->build_clusters(grid_main, marked_cells, clusters);

			// 3.2 Feines Gitter für jedes Cluster der Ebene erzeugen
			for (unsigned int i = 0; i < clusters.size(); i++) {
				cout << "--> Creating fine grid #" << i << endl;
				create_fine_grid(clusters[i], coarse_grid);

				// 4. Rekursiv feinere Gitter berechnen //TODO: Rekursiv
				int level = 1; //TODO: Remove, sollte Teil der Datenstruktur sein
				int dt_steps = pow(2,(level));
				for (int l = 0; l < dt_steps; l++){
					// Flussberechnen auf feinen Gitter

					// Boundary Condtions auf durchlässig -> //TODO: Verbessern
					clusters[i].grid_fine->boundary_conditions[0]= 0;
					clusters[i].grid_fine->boundary_conditions[1]= 0;
					clusters[i].grid_fine->boundary_conditions[2]= 0;
					clusters[i].grid_fine->boundary_conditions[3]= 0;


					solver->calc_method_flux(dt / dt_steps, clusters[i].grid_fine);
					clusters[i].grid_fine->apply_boundary_conditions();
					if(l==dt_steps/2-1){
						solver->get_2d_half_flux(&Fine_Half1.half_flux, Fine_Half1.dimension, Fine_Half1.size);
					}
					else if(l==dt_steps-1){
						solver->get_2d_half_flux(&Fine_Half2.half_flux, Fine_Half2.dimension, Fine_Half2.size);
					}
				}
				// 5. Gröberes Gitter mit feinerem Gitter verbessern
				fine_to_coarse(clusters[i],dt);

				//TODO: Optional, remove und mit write Funtkion ersetzen
				/***/ std::ostringstream oss;
				/***/ oss << "amr_refined" << n << ".txt";
				/***/ string filename = oss.str();
				/***/ grid_to_file(filename, grid_main);
				// -------------

				cout << "<-- Adjusted coarse grid with fine grid #" << i << endl;
			}

			// **** REMOVE FROM HERE *************
			for(int i=0;i<Coarse_Half.dimension;i++){
				delete[] Coarse_Half.half_flux[i];
			}
			delete[] Coarse_Half.half_flux;

			for(int i=0;i<Fine_Half1.dimension;i++){
				delete[] Fine_Half1.half_flux[i];
			}
			delete[] Fine_Half1.half_flux;

			for(int i=0;i<Fine_Half2.dimension;i++){
				delete[] Fine_Half2.half_flux[i];
			}
			delete[] Fine_Half2.half_flux;
			// ******************************

			cout << "Coarse grid adjusted! Next time step." << endl;
			//cin.get();
			// 6. Next
		} else {
			cout << "No further refinement needed! Skipping to next time step!" << endl;
		}
		delete[] marked_cells;
	}
	//----------------------------------------------
}

void Adaptive_Mesh::calculate_grids(Grid* refinable_grid, int grid_level, double & dt){

	// Mit jedem Level werden die Zeitschritte weiter halbiert, also doppelt so viele sind nötig
	int calc_steps = pow(2,grid_level);
	cout << "Calc steps: "<<calc_steps<<endl;

	dt = time_calculation->cfl_condition(refinable_grid);
	double dt_stepped = dt/calc_steps;


	// 1.1 Gröbstes Netz berechnen

	if(refinable_grid->copy_to(grid_twostep) == -1){
		cout << "Failed to copy Main grid";
		exit(EXIT_FAILURE);
	}

	solver->calc_method_flux(dt_stepped/2,grid_twostep);
	grid_twostep->apply_boundary_conditions();
	solver->calc_method_flux(dt_stepped/2,grid_twostep);
	grid_twostep->apply_boundary_conditions();

	// Auf Ebene 0 ist der "Doppelschritt" identisch zum groben gitter zu t=t+dt
	if(grid_level == 0){
		solver->calc_method_flux(dt,refinable_grid);
		refinable_grid->apply_boundary_conditions();
		grid_doublestep = refinable_grid;

	}
	else{
		for(int i = 0; i<calc_steps;i++){
			solver->calc_method_flux(dt_stepped,refinable_grid);
			refinable_grid->apply_boundary_conditions();

			// Doublestep ist Hauptgitter nach einem "level-abhängigen" Zeitschritt
			if(i == 0){
				if(refinable_grid->copy_to(grid_doublestep) == -1){
					cout << "Failed to copy Grid";
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	solver->get_2d_half_flux(&Coarse_Half.half_flux, Coarse_Half.dimension, Coarse_Half.size);
	//TODO: Remove
	if(&Coarse_Half.half_flux != nullptr)
		cout << "Copied half flux"<<endl;
}

void Adaptive_Mesh::fine_to_coarse(Cluster_Square& cluster, double dt){
	Grid* coarse = cluster.parent;
	Grid* fine = cluster.grid_fine;

	int x_half = 0;
	int y_half = 0;
	int pos = 0;
	int fine_pos_00 = 0;
	int fine_pos_01 = 0;
	int fine_pos_10 = 0;
	int fine_pos_11 = 0;

	// TODO: Sollte nie NULL sein, sobald half_flux für alle Solvers implementiert ist
	if(Coarse_Half.half_flux != NULL){
		x_half = 2;
		y_half = 2;
		for (int x = cluster.pos_x_min+1; x <= cluster.pos_x_max-1; x++) {
			y_half=2;
			for (int y = cluster.pos_y_min+1; y <= cluster.pos_y_max-1; y++) {
				pos = (x) + (y) * coarse->grid_size_total[0];
				fine_pos_00 = (x_half) + (y_half) * fine->grid_size_total[0];
				fine_pos_01 = (x_half) + (y_half+1) * fine->grid_size_total[0];
				fine_pos_10 = (x_half+1) + (y_half) * fine->grid_size_total[0];
				fine_pos_11 = (x_half+1) + (y_half+1) * fine->grid_size_total[0];

				for (int n = 0; n<=5; n++){
					coarse->cellsgrid[pos][n] =
							0.25 * (fine->cellsgrid[fine_pos_00][n] + fine->cellsgrid[fine_pos_01][n] + fine->cellsgrid[fine_pos_10][n] + fine->cellsgrid[fine_pos_11][n]);
				}
				y_half=y_half+2;
			}
			x_half = x_half + 2;
		}
		coarse_correction(cluster,dt);
	}
	else{
		x_half = 0;
		y_half = 0;
		for (int x = cluster.pos_x_min; x <= cluster.pos_x_max; x++) {
			y_half=0;
			for (int y = cluster.pos_y_min; y <= cluster.pos_y_max; y++) {
				pos = (x) + (y) * coarse->grid_size_total[0];
				fine_pos_00 = (x_half) + (y_half) * fine->grid_size_total[0];
				fine_pos_01 = (x_half) + (y_half+1) * fine->grid_size_total[0];
				fine_pos_10 = (x_half+1) + (y_half) * fine->grid_size_total[0];
				fine_pos_11 = (x_half+1) + (y_half+1) * fine->grid_size_total[0];

				for (int n = 0; n<=5; n++){
					coarse->cellsgrid[pos][n] =
							0.25 * (fine->cellsgrid[fine_pos_00][n] + fine->cellsgrid[fine_pos_01][n] + fine->cellsgrid[fine_pos_10][n] + fine->cellsgrid[fine_pos_11][n]);
				}
				y_half=y_half+2;
			}
			x_half = x_half + 2;
		}
	}

}

// Korrektur der Geisterzellen auf groben Gitter
void Adaptive_Mesh::coarse_correction(Cluster_Square& cluster, double dt){

	Grid* coarse = cluster.parent;
	Grid* fine = cluster.grid_fine;

	//TODO: Geister über Korrektur
	//TODO: Das ganze ist eine Differenz zum Coarse U N+1 Wert
	double dtodx = dt/dx;
	double dtody = dt/dy;

	double d = 0;
	double p = 0;
	double uxd = 0;
	double uxr = 0;
	double uyd = 0;
	double uyr = 0;

	int index_coarse_1 = 0;
	int index_coarse_2 = 0;
	int index_fine_a1 = 0;
	int index_fine_b1 = 0;
	int index_fine_2 = 0;

	int x_half = 0;
	int y_half = 0;
	int pos = 0;

	int x = cluster.pos_x_min;
	if (x != 0) { // Nur für Coarsest Grid?
		for (int y = cluster.pos_y_min + 1; y <= cluster.pos_y_max - 1; y++) {
			pos = (x) + (y) * coarse->grid_size_total[0];

			index_coarse_1 = y + (coarse->grid_size_total[1] - 1) * (x - 1);
			index_coarse_2 = (coarse->grid_size_total[1] - 1) * (coarse->grid_size_total[0] - 1);
			index_fine_a1 = y_half + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_b1 = y_half + (fine->grid_size_total[1] - 1) * (x_half + 1);
			index_fine_2 = (fine->grid_size_total[1] - 1) * (fine->grid_size_total[0] - 1);

			//TODO: Index über [dim][pos][erhaltungsgröße] !?

			d = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (0)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (0)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (0)]);

			uxd = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (1)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (1)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (1)]);

			uyd = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (2)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (2)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (2)]);

			uxr = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (3)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (3)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (3)]);

			uyr = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (4)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (4)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (4)]);

			coarse->cellsgrid[pos][0] = coarse->cellsgrid[pos][0] - d;
			p = constants->ct * pow(coarse->cellsgrid[pos][0], constants->gamma);
			coarse->cellsgrid[pos][1] = coarse->cellsgrid[pos][1] - p;
			coarse->cellsgrid[pos][2] = coarse->cellsgrid[pos][2] - (uxd / coarse->cellsgrid[pos][0]); // TODO: So richtig? Oder uxd/d ?
			coarse->cellsgrid[pos][4] = coarse->cellsgrid[pos][4] - (uyd / coarse->cellsgrid[pos][0]);
			coarse->cellsgrid[pos][3] = coarse->cellsgrid[pos][3] - uxr;
			coarse->cellsgrid[pos][5] = coarse->cellsgrid[pos][5] - uyr;

			y_half = y_half + 2;
		}
	}

	x = cluster.pos_x_max;
	x_half = (cluster.pos_x_max-cluster.pos_x_min)*2-2; // TODO: -2? Warum
	y_half = 0;
	pos = 0;

	if (x != coarse->grid_size_total[0]-1) { // -1?
		for (int y = cluster.pos_y_min + 1; y <= cluster.pos_y_max - 1; y++) {
			pos = (x) + (y) * coarse->grid_size_total[0];

			index_coarse_1 = y + (coarse->grid_size_total[1] - 1) * (x); // x + 1 ?
			index_coarse_2 = (coarse->grid_size_total[1] - 1) * (coarse->grid_size_total[0] - 1);
			index_fine_a1 = y_half + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_b1 = y_half + (fine->grid_size_total[1] - 1) * (x_half + 1);
			index_fine_2 = (fine->grid_size_total[1] - 1) * (fine->grid_size_total[0] - 1);

			//TODO: Index über [dim][pos][erhaltungsgröße] !?

			d = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (0)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (0)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (0)]);

			uxd = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (1)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (1)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (1)]);

			uyd = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (2)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (2)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (2)]);

			uxr = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (3)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (3)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (3)]);

			uyr = dtodx * Coarse_Half.half_flux[0][index_coarse_1 + index_coarse_2 * (4)]
					- 0.25 * dtodx
							* (Fine_Half1.half_flux[0][index_fine_a1 + index_fine_2 * (4)] + Fine_Half1.half_flux[0][index_fine_b1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[0][index_fine_a1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[0][index_fine_b1 + index_fine_2 * (4)]);

			coarse->cellsgrid[pos][0] = coarse->cellsgrid[pos][0] - d;
			p = constants->ct * pow(coarse->cellsgrid[pos][0], constants->gamma);
			coarse->cellsgrid[pos][1] = coarse->cellsgrid[pos][1] + p;
			coarse->cellsgrid[pos][2] = coarse->cellsgrid[pos][2] + (uxd / coarse->cellsgrid[pos][0]); // TODO: So richtig? Oder uxd/d ?
			coarse->cellsgrid[pos][4] = coarse->cellsgrid[pos][4] + (uyd / coarse->cellsgrid[pos][0]);
			coarse->cellsgrid[pos][3] = coarse->cellsgrid[pos][3] + uxr;
			coarse->cellsgrid[pos][5] = coarse->cellsgrid[pos][5] + uyr;

			y_half = y_half + 2;
		}
	}

	x_half = 0;
	y_half = 0;
	pos = 0;

	int y = cluster.pos_y_min;
	if (y != 0) { // Nur für Coarsest Grid?
		for (int x = cluster.pos_x_min+1; x <= cluster.pos_x_max-1; x++) { // +1 und -1 Warum?
			pos = (x) + (y) * coarse->grid_size_total[0];

			index_coarse_1 = (y-1) + (coarse->grid_size_total[1] - 1) * (x);
			index_coarse_2 = (coarse->grid_size_total[1] - 1) * (coarse->grid_size_total[0] - 1);
			index_fine_a1 = y_half + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_b1 = (y_half+1) + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_2 = (fine->grid_size_total[1] - 1) * (fine->grid_size_total[0] - 1);

			//TODO: Index über [dim][pos][erhaltungsgröße] !?

			d = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (0)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (0)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (0)]);

			uxd = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (1)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (1)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (1)]);

			uyd = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (2)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (2)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (2)]);

			uxr = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (3)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (3)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (3)]);

			uyr = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (4)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (4)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (4)]);

			coarse->cellsgrid[pos][0] = coarse->cellsgrid[pos][0] - d;
			p = constants->ct * pow(coarse->cellsgrid[pos][0], constants->gamma);
			coarse->cellsgrid[pos][1] = coarse->cellsgrid[pos][1] - p;
			coarse->cellsgrid[pos][2] = coarse->cellsgrid[pos][2] - (uxd / coarse->cellsgrid[pos][0]); // TODO: So richtig? Oder uxd/d ?
			coarse->cellsgrid[pos][4] = coarse->cellsgrid[pos][4] - (uyd / coarse->cellsgrid[pos][0]);
			coarse->cellsgrid[pos][3] = coarse->cellsgrid[pos][3] - uxr;
			coarse->cellsgrid[pos][5] = coarse->cellsgrid[pos][5] - uyr;

			x_half = x_half + 2;
		}
	}

	y = cluster.pos_y_max;
	y_half = (cluster.pos_y_max-cluster.pos_y_min)*2-2; // TODO: -2? Warum
	x_half = 0;
	pos = 0;

	if (y != coarse->grid_size_total[1]-1) { // TODO: -1? Warum
		for (int x = cluster.pos_x_min + 1; x <= cluster.pos_x_max - 1; x++) {
			pos = (x) + (y) * coarse->grid_size_total[0];

			index_coarse_1 = (y) + (coarse->grid_size_total[1] - 1) * (x);
			index_coarse_2 = (coarse->grid_size_total[1] - 1) * (coarse->grid_size_total[0] - 1);
			index_fine_a1 = y_half + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_b1 = (y_half+1) + (fine->grid_size_total[1] - 1) * (x_half);
			index_fine_2 = (fine->grid_size_total[1] - 1) * (fine->grid_size_total[0] - 1);

			//TODO: Index über [dim][pos][erhaltungsgröße] !?

			d = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (0)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (0)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (0)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (0)]);

			uxd = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (1)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (1)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (1)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (1)]);

			uyd = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (2)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (2)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (2)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (2)]);

			uxr = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (3)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (3)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (3)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (3)]);

			uyr = dtody * Coarse_Half.half_flux[1][index_coarse_1 + index_coarse_2 * (4)]
					- 0.25 * dtody
							* (Fine_Half1.half_flux[1][index_fine_a1 + index_fine_2 * (4)] + Fine_Half1.half_flux[1][index_fine_b1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[1][index_fine_a1 + index_fine_2 * (4)]
									+ Fine_Half2.half_flux[1][index_fine_b1 + index_fine_2 * (4)]);

			coarse->cellsgrid[pos][0] = coarse->cellsgrid[pos][0] - d;
			p = constants->ct * pow(coarse->cellsgrid[pos][0], constants->gamma);
			coarse->cellsgrid[pos][1] = coarse->cellsgrid[pos][1] + p;
			coarse->cellsgrid[pos][2] = coarse->cellsgrid[pos][2] + (uxd / coarse->cellsgrid[pos][0]); // TODO: So richtig? Oder uxd/d ?
			coarse->cellsgrid[pos][4] = coarse->cellsgrid[pos][4] + (uyd / coarse->cellsgrid[pos][0]);
			coarse->cellsgrid[pos][3] = coarse->cellsgrid[pos][3] + uxr;
			coarse->cellsgrid[pos][5] = coarse->cellsgrid[pos][5] + uyr;

			x_half = x_half + 2;
		}
	}
}

void Adaptive_Mesh::create_fine_grid(Cluster_Square& cluster, Grid* parent_grid) {
	// TODO: Parent_grid vorher zuweisen!?
	// TODO: Feinegitter Werte von Schritt vorher?

	// Feines Gitter mit doppelter Anzahl Zellen
	int grid_size_x = (cluster.pos_x_max - cluster.pos_x_min + 1)*2;
	int grid_size_y = (cluster.pos_y_max - cluster.pos_y_min + 1)*2;

	// Feines Gitter erstellen und Cluster zuweisen
	Grid* adaptive_grid = new Grid(grid_size_x,grid_size_y, constants); //TODO: DELETE
	cluster.grid_fine = adaptive_grid;
	cluster.parent = parent_grid;

	// Mittelwerte Zellwände (ecken und kanten)
	double half_steps[grid_size_x+1][grid_size_y+1][6]; // n,m = grade -> Zellmittelpunkt (grob) || ungrade -> Zelleckpunkt (grob)

	// half_steps[0][0] bei erster grober Zelle (Mittelpunkt), [1][1] bei (oben rechts) Zellecke [n][m] bei Zellmittelpunkt "grobe Zelle oben rechts"
	// Ausnahme bei an oberen und rechten Rand
	int x_half = 0;
	int y_half = 0;

	for (int x = cluster.pos_x_min; x < cluster.pos_x_max; x++) {
		y_half=0;
		for (int y = cluster.pos_y_min; y < cluster.pos_y_max; y++) {	// TODO: CHeck, vielleicht <= ?
			for (int n = 0; n<=5; n++){
				// Zellemitte, rechte Wand, obere Wand, oben rechts Ecke
				half_steps[x_half][y_half][n] = parent_grid->cellsgrid[(x) + (y) * parent_grid->grid_size_total[0]][n];

				half_steps[x_half+1][y_half][n] = 0.5 * (half_steps[x_half][y_half][n]
													+ parent_grid->cellsgrid[(x+1) + (y) * parent_grid->grid_size_total[0]][n]);

				half_steps[x_half][y_half+1][n] = 0.5 * (half_steps[x_half][y_half][n]
													+ parent_grid->cellsgrid[(x) + (y+1) * parent_grid->grid_size_total[0]][n]);

				half_steps[x_half+1][y_half+1][n] = 0.25 * (half_steps[x_half][y_half][n]
														+ parent_grid->cellsgrid[(x) + (y+1) * parent_grid->grid_size_total[0]][n]
														+ parent_grid->cellsgrid[(x+1) + (y) * parent_grid->grid_size_total[0]][n]
														+ parent_grid->cellsgrid[(x+1) + (y+1) * parent_grid->grid_size_total[0]][n]);
			}
			y_half=y_half+2;
		}
		x_half=x_half+2;
	}

	// "oberer Rand" nur Zellhälften nach rechts
	x_half = 0;
	y_half = (cluster.pos_y_max-cluster.pos_y_min)*2;
	for (int x = cluster.pos_x_min; x < cluster.pos_x_max; x++) {
		for (int n = 0; n<=5; n++){
			half_steps[x_half][y_half][n] = parent_grid->cellsgrid[(x) + (cluster.pos_y_max) * parent_grid->grid_size_total[0]][n];
			half_steps[x_half+1][y_half][n] = 0.5 * (half_steps[x_half][y_half][n]
												+ parent_grid->cellsgrid[(x+1) + (cluster.pos_y_max) * parent_grid->grid_size_total[0]][n]);
		}
		x_half=x_half+2;
	}
	for (int n = 0; n<=5; n++){
		half_steps[x_half][y_half][n] = parent_grid->cellsgrid[(cluster.pos_x_max) + (cluster.pos_y_max) * parent_grid->grid_size_total[0]][n];
	}

	// "rechter Rand" nur Zellhälften nach oben
	x_half = (cluster.pos_x_max-cluster.pos_x_min)*2;
	y_half = 0;
	for (int y = cluster.pos_y_min; y < cluster.pos_y_max; y++) {
		for (int n = 0; n<=5; n++){
			half_steps[x_half][y_half][n] = parent_grid->cellsgrid[(cluster.pos_x_max) + (y) * parent_grid->grid_size_total[0]][n];
			half_steps[x_half][y_half+1][n] = 0.5 * (half_steps[x_half][y_half][n]
												+ parent_grid->cellsgrid[(cluster.pos_x_max) + (y+1) * parent_grid->grid_size_total[0]][n]);
		}
		y_half=y_half+2;
	}

	// Feine Zellmitten über gemittlete Werte an "Ecken"
	for (int x = 1; x < adaptive_grid->grid_size_total[0]-1; x++) {
		for (int y = 1; y < adaptive_grid->grid_size_total[1]-1; y++) {
			for (int n = 0; n <= 5; n++) {
				adaptive_grid->cellsgrid[(x) + (y) * adaptive_grid->grid_size_total[0]][n] = 0.25 * (half_steps[x-1][y-1][n] + half_steps[x-1][y][n] + half_steps[x][y-1][n] + half_steps[x][y][n]);
			}
		}
	}

	// Äußerste Geisterzellen "unten" setzen auf benachbarten Zellwert
	for (int y = 1; y < adaptive_grid->grid_size_total[1]-1; y++) {
		for (int n = 0; n <= 5; n++) {
			adaptive_grid->cellsgrid[(0) + (y) * adaptive_grid->grid_size_total[0]][n] = adaptive_grid->cellsgrid[(1) + (y) * adaptive_grid->grid_size_total[0]][n];
		}
	}
	// Äußerste Geisterzellen "oben" setzen auf benachbarten Zellwert
	for (int y = 1; y < adaptive_grid->grid_size_total[1]-1; y++) {
		for (int n = 0; n <= 5; n++) {
			adaptive_grid->cellsgrid[(adaptive_grid->grid_size_total[0]-1) + (y) * adaptive_grid->grid_size_total[0]][n] = adaptive_grid->cellsgrid[(adaptive_grid->grid_size_total[0]-2) + (y) * adaptive_grid->grid_size_total[0]][n];
		}
	}
	// Äußerste Geisterzellen "links" setzen auf benachbarten Zellwert
	for (int x = 0; x < adaptive_grid->grid_size_total[0]; x++) {
		for (int n = 0; n <= 5; n++) {
			adaptive_grid->cellsgrid[(x) + (0) * adaptive_grid->grid_size_total[0]][n] = adaptive_grid->cellsgrid[(x) + (1) * adaptive_grid->grid_size_total[0]][n];
		}
	}
	// Äußerste Geisterzellen "rechts" setzen auf benachbarten Zellwert
	for (int x = 0; x < adaptive_grid->grid_size_total[0]; x++) {
		for (int n = 0; n <= 5; n++) {
			adaptive_grid->cellsgrid[(x) + (adaptive_grid->grid_size_total[1]-1) * adaptive_grid->grid_size_total[0]][n] = adaptive_grid->cellsgrid[(x) + (adaptive_grid->grid_size_total[1]-2) * adaptive_grid->grid_size_total[0]][n];
		}
	}

	//TODO: Remove
	 ofstream myfile ("amr_finegrid.txt");
	  if (myfile.is_open())
	  {
		for (int y = 0; y < adaptive_grid->grid_size_total[1]; y++) {
			for (int x = 0; x < adaptive_grid->grid_size_total[0]; x++) {
				int pos = x + y * adaptive_grid->grid_size_total[0];
				myfile << x << " " << y << " " << adaptive_grid->cellsgrid[pos][0] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";
}

int Adaptive_Mesh::grid_marker(int* &marked_cells, Grid* grid_one, Grid* grid_two, double tolerance){
	// Markiert Zellen mit -1
	int pos = 0;
	int is_marked = 0;
	for (int y = 0; y < grid_one->grid_size_total[1]; y++) {
		for (int x = 0; x < grid_one->grid_size_total[0]; x++) {
				pos = x + y * grid_one->grid_size_total[0];
				for(int k=0;k<6;k++){
					if(fabs(grid_one->cellsgrid[pos][k]-grid_two->cellsgrid[pos][k]) > tolerance * fabs(grid_one->cellsgrid[pos][k])){
						//cout << "Checked: "<< x<<","<<y<<" for "<<k<<endl;
						//cout << grid_one->cellsgrid[pos][k] << "-"<< grid_two->cellsgrid[pos][k] <<"="<< fabs(grid_one->cellsgrid[pos][k]-grid_two->cellsgrid[pos][k]) << " > " << tolerance*fabs(grid_one->cellsgrid[pos][k])<<endl;
						for (int l = y-2; l<=y+2;l++){
							for (int k = x-2; k<=x+2;k++){
								pos = k + l * grid_one->grid_size_total[0];
								if(k<0 || k >= grid_one->grid_size_total[0] || l<0 || l >= grid_one->grid_size_total[1]){
									//cout << "Couldn't check, out of bounds: "<< k<<","<<l<<endl;
									continue;
								}
								marked_cells[pos] = -1;
								is_marked = 1;
							}
						}

						break;
					}
				}
			}
	}

	//TODO: Remove
	 ofstream myfile ("amr_marked.txt");
	  if (myfile.is_open())
	  {
		for (int y = 0; y < grid_one->grid_size_total[1]; y++) {
			for (int x = 0; x < grid_one->grid_size_total[0]; x++) {
				pos = x + y * grid_one->grid_size_total[0];
				myfile << x << " " << y << " " << marked_cells[pos] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";

	  return is_marked;
}

// TODO: For testing purposes only
void Adaptive_Mesh::grid_to_file(string filename, Grid* grid){
	 ofstream myfile (filename);
	  if (myfile.is_open())
	  {
		for (int y = 0; y < grid->grid_size_total[1]; y++) {
			for (int x = 0; x < grid->grid_size_total[0]; x++) {
				int pos = x + y * grid->grid_size_total[0];
				myfile << x << " " << y << " " << grid->cellsgrid[pos][3] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";
}
