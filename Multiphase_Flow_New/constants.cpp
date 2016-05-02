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
				if (name == "variante")
					variante = value;
				if (name == "teiler")
					teiler = value;
				if (name == "teilerend")
					teilerend = value;
				if (name == "timeou")
					timeou = value;
				if (name == "cfl")
					cfl = value;
				if (name == "maxnt")
					maxnt = value;
				if (name == "dimension")
					dimension = value;
				if (name == "ordnung")
					ordnung = value;
				if (name == "radius")
					radius = value;
				if (name == "CELLSX")
					CELLSX = value;
				if (name == "CELLSY")
					CELLSY = value;
				if (name == "g")
					g = value;
				if (name == "mol")
					mol = value;
				if (name == "mor")
					mor = value;
				if (name == "mul")
					mul = value;
				if (name == "mur")
					mur = value;
				if (name == "upbc")
					upbc = value;
				if (name == "downbc")
					downbc = value;
				if (name == "leftbc")
					leftbc = value;
				if (name == "rightbc")
					rightbc = value;
				if (name == "cref")
					cref = value;
				if (name == "done")
					done = value;
				if (name == "ccl")
					ccl = value;
				if (name == "rhol")
					rhol = value;
				if (name == "vl")
					vl = value;
				if (name == "vrl")
					;
				vrl = value;
				if (name == "vyl")
					;
				vyl = value;
				if (name == "vyrl")
					vyrl = value;
				if (name == "rhor")
					rhor = value;
				if (name == "vr")
					vr = value;
				if (name == "vrr")
					vrr = value;
				if (name == "vyr")
					vyr = value;
				if (name == "vyrr")
					vyrr = value;

				const_name.push_back(name);
				const_value.push_back(value);
			}
			input.close();
		} else {
			cout << "Fehler beim öffnen der Konstanten!";
		}
	}

}

Constants& Constants::myinstance() {
	static Constants constants;
	return constants;
}
