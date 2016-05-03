#ifndef SOLVER_H_
#define SOLVER_H_

#include "Gleichungssystem.h"

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
		Solver(string methodname, Constants* constants, Computation *computation, Grid* grid);

		/**
		* Abstrakte Methode zur berechnung des Flusses der jeweiligen numerischen Methode.
		* @return Matrix der Flüsse (1D)
		*/
		virtual std::vector< std::vector< std::vector< std::vector <double> > > >calc_method_flux(double dt, int dir) =0;

        /**
        * Konstanten Objekt welches für die berechnungen benötigt wird.
        * @see Konstanten
        */
        Constants *konstanten;
        /**
        * Raster in den gerechnet wird.
        * @see Raster
        */
        Grid *grid;
        /**
        * Gleichungssystem Objekt.
        * @see Gleichungssystem
        */
        Computation *gs;
        string name;
        int* CELLS;

		/**
		* Zeit Output
		*/
		double timeou;
		/**
		* Dimension in der gerechnet wird.
		*/
		int dimension;
		/**
		* ordnung des Verfahrens.
		*/
		int ordnung;
		/**
		* Anzahl an Schritten die gemacht wurden.
		*/
		int steps;
		/**
		* Gesetztes Maximum, damit die Methode nicht unendlich läuft
		* (falls es einen fehler gibt oder andere umstände).
		*/
		int maxnt;
		/**
		* Delta x.
		*/
		double dx;
	/**
		* teiler, Faktor für die ersten Delta t Schritte
		*/
		double teiler;
	/**
		* teilerend, Ende der Multiplikation der Zeitschritte mit teiler
		*/
		double teilerend;
	/**
		* Delta y für die 2. Dimension
		*/
		double dy;
		/**
		* Linke Grenze.
		*/
		int mol;
		/**
		* Rechte Grenze.
		*/
		int mor;
		/**
		* Obere Grenze.
		*/
		int mul;
		/**
		* Untere Grenze.
		*/
		int mur;
		/**
		* Variante der EOS.
		*/
		int variante;

};


#endif /* SOLVER_H_ */
