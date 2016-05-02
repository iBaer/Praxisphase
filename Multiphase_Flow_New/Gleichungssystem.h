#ifndef COMPUTATION_H
#define COMPUTATION_H
#include "Zelle.h"
#include "Raster.h"
#include <iostream>
#include <vector>

#include "constants.h"

/*!
* @class Computation
* Die Klasse Computation verwaltet die Formelvektoren u,f_u und s_u.
* Die Formelvektoren wurden in der ursprüngliche Version als Strings in 
* Vektoren abgespeichert, und dann mittels Funktionen der exprtk Libiary 
* interpretiert, aber das dauerte viel zu lange
*/

class Computation
{
    private:

        /**
        * Konstruktor der klasse Gleichungssystem.
        * @param path Pfad zur Datei welche die Formeln enthält.
        * @param c Konstanten Objekte welches für berechnungen gebraucht wird.
        */
        Computation(Constants *constants);

    public:

        static Computation& myinstance(Constants *constants);

        /**
        * Vektor von u
        */
        std::vector<std::string> u;
        /**
        * Rückrechnung der erhaltenen Variablen in die physikalischen.
        */
        std::vector<std::string> uback;
        /**
        * Vektor der Formeln von f_u.
        */
        std::vector<std::string> f_u;
        /**
        * Vektor der Formeln von g_u.
        */
        std::vector<std::string> g_u;

	/**
	 * Berechnet die Werte U in 1-D
	 * @param u, 3-d Feld für die U-Werte
	 * @param Raster, welche die Werte zur Berechnung enthält.
	 * @param cells, Ausdehnung des Gitters
	 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
	 */
	void compute_u_1d(double *** u, Raster * raster, int * cells, int ordnung);

	/**
	 * Berechnet die Lösungen der Formel für den Fluss F in 1-D
	 * @param f, 3-d Feld für die F-Werte
	 * @param Raster, welche die Werte zur Berechnung enthält.
	 * @param cells, Ausdehnung des Gitters
	 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
	 */
	void compute_f_1d(double *** f, Raster * raster, int * cells, int ordnung);

        /**
        * Berechnet die Werte U in 2-D
        * @param u, 3-d Feld für die U-Werte
        * @param Raster, welche die Werte zur Berechnung enthält.
        * @param cells, Ausdehnung des Gitters
        * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
        */
	void compute_u_2d(double *** u, Raster * raster, int * cells, int ordnung);
	
        /**
        * Berechnet die Lösungen der Formel für den Fluss F in 2-D
        * @param f, 3-d Feld für die F-Werte
        * @param Raster, welche die Werte zur Berechnung enthält.
        * @param cells, Ausdehnung des Gitters
        * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
        */
	void compute_f_2d(double *** f, Raster * raster, int * cells, int ordnung);

        /**
        * Berechnet die Lösungen der Formel für den Fluss G in 2-D
        * @param Raster, welche die Werte zur Berechnung enthält.
        * @param cells, Ausdehnung des Gitters
        * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
        */
	void compute_g_2d(double *** g, Raster * raster, int * cells, int ordnung);

        /**
        * Anzahl der Gleichungen
        */
 	int neqs;
	
        /**
        * Abgespeichertes Konstanten Objekt.
        * Wird nur einmal beim Konstruktor gesetzt, damit
        * es bei Rechnungen nicht immer neu übergeben werden muss.
        */

 	Constants* c;

	double cref;
	double done;
	double ccl;
	double gc;

	double ccl12;
	double ccl12h;
	double gcinv;
	double gcdgcm1;
	double gcm1dgc;
	double cclm1;
	double powcref;
};

#endif // COMPUTATION_H
