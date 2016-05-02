#ifndef SOLVER_H_
#define SOLVER_H_

#include "Gleichungssystem.h"

#include "numerische_methode.h"
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>

#include "constants.h"

using namespace std;

/*!
* @class Solver
* Abstrakte Klasse für sämtliche Methoden zur Lösung von Differential Gleichungen
*/


class Solver
{
	public:
		/**
		* Konstruktor.
		* @param dim Setzt Dimension für die Berechnung.
		* @param ordn Ordnung der numerischen Methode.
		* @param cells Anzahl der Zellen.
		* @param method Name der Methode für unterscheidung bei Output.
		*/
		Solver(std::string method, std::string const_in, std::string formel_in, std::string save_in);

		/**
		* Abstrakte Methode zur berechnung des Flusses der jeweiligen numerischen Methode.
		* @return Matrix der Flüsse (1D)
		*/
		virtual std::vector< std::vector< std::vector< std::vector <double> > > >calc_method_flux(int dir) =0;

        /**
        * Konstanten Objekt welches für die berechnungen benötigt wird.
        * @see Konstanten
        */
        Constants *konstanten;
        /**
        * Raster in den gerechnet wird.
        * @see Raster
        */
        //Raster raster;
        /**
        * Gleichungssystem Objekt.
        * @see Gleichungssystem
        */
        Computation *gs;
};


#endif /* SOLVER_H_ */
