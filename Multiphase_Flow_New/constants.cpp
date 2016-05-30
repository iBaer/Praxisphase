#include "constants.h"

#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

/**
 *****************************************************************************************
 * Konstruktor.
 * @param input_const String des Pfades zu den festen Konstanten.
 * @param input_calc String des Pfades zu den einmalig zu berechnenden Konstanten.
 *****************************************************************************************/
Constants::Constants() {
	input_const = "gas-liquid.in";
	ifstream input(input_const.c_str());
	string line, name, calc;
	double value;

	if (!input_const.empty()) {
		if (input.is_open()) {
			//einlesen Leerzeichen als Trenzeichen
			while (getline(input, line)) {
				name = line.substr(0, line.find(" "));
				value = atof(
						line.substr(line.find(" ") + 1,
								line.size() - line.find(" ") - 1).c_str());
				if (name == "calceigv") {
					calceigv = value;
					cout << "Wert von calceigv ist " << value << endl;
				}
				else if (name == "variante")
					variante = value;
				else if (name == "teiler")
					teiler = value;
				else if (name == "teilerend")
					teilerend = value;
				else if (name == "timeou")
					timeou = value;
				else if (name == "cfl")
					cfl = value;
				else if (name == "maxnt")
					maxnt = value;
				else if (name == "dimension")
					dimension = value;
				else if (name == "ordnung")
					order = value;
				else if (name == "radius")
					radius = value;
				else if (name == "CELLSX")
					grid_size_x = value;
				else if (name == "CELLSY")
					grid_size_y = value;
				else if (name == "g")
					gamma = value;
				else if (name == "mol")
					pos_x_min = value;
				else if (name == "mor")
					pos_x_max = value;
				else if (name == "mul")
					pos_y_min = value;
				else if (name == "mur")
					pos_y_max = value;
				else if (name == "upbc")
					bc_y_max = value;
				else if (name == "downbc")
					bc_y_min = value;
				else if (name == "leftbc")
					bc_x_min = value;
				else if (name == "rightbc")
					bc_x_max = value;
				else if (name == "cref")
					cref = value;
				else if (name == "done")
					rho_one = value;
				else if (name == "ccl")
					ccl = value;
				else if (name == "rhol")
					rhol = value;
				else if (name == "vl")
					vl = value;
				else if (name == "vrl")
					vrl = value;
				else if (name == "vyl")
					vyl = value;
				else if (name == "vyrl")
					vyrl = value;
				else if (name == "rhor")
					rhor = value;
				else if (name == "vr")
					vr = value;
				else if (name == "vrr")
					vrr = value;
				else if (name == "vyr")
					vyr = value;
				else if (name == "vyrr")
					vyrr = value;

				//const_name.push_back(name);
				//const_value.push_back(value);
			}
			input.close();
		} else {
			cout << "Fehler beim öffnen der Konstanten!";
		}
	}

    double alfll = 1.0 - rhol/rho_one + ccl *(rhol/rho_one);
    double dll = ccl * (rhol/alfll);
    double pll = cref*pow(dll,gamma);

    ct = pll/pow(rhol,gamma);

}

Constants& Constants::instance() {
	static Constants constants;
	return constants;
}
