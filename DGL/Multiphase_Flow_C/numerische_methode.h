#ifndef NUMERISCHE_METHODE_H
#define NUMERISCHE_METHODE_H
#include "Raster.h"
#include "Gleichungssystem.h"
#include <iostream>

/*!
* @class numerische_methode
* Die abstrakte Klasse numerische_methode gibt den Rahmen
* für alle erbenden numerischen Methoden vor.
*/

class numerische_methode
{
    public:
        /**
        * Konstruktor.
        * @param dim Setzt Dimension für die Berechnung.
        * @param ordn Ordnung der numerischen Methode.
        * @param cells Anzahl der Zellen.
        * @param method Name der Methode für unterscheidung bei Output.
        */
        numerische_methode(std::string method, std::string const_in, std::string formel_in, std::string save_in);
        /**
        * Startet die Berechnung.
        */
        void start_method();
    protected:
        /**
        * Name der Methode.
        * Wird hauptsächlich dazu benötigt um die
        * die Outputs der verschiedenen Methoden ausseinander
        * halten zu können.
        */
        std::string name;
        /**
        * Array welches die Anzahl der Zellen in der entsprechenden Dimension zeigt.
        */
        int* CELLS;
        /**
        * Gamma Konstante
        */
        double g;
        /**
        * Zeit Output
        */
        double timeou;
        /**
        * K in den Formeln.
        */
        double ct;
        /**
        * Dimension in der gerechnet wird.
        */
        int dimension;
        /**
        * Bei dimensions==2, Integrationsschema splitting oder unsplitting
        */
        int splitting;
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
        * Delta t.
        */
        double dt;
        /**
        * Delta t_fixed.
        */
        double dt_fixed;
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
        /**
        * Konstanten Objekt welches für die berechnungen benötigt wird.
        * @see Konstanten
        */
        Konstanten konstanten;
        /**
        * Raster in den gerechnet wird.
        * @see Raster
        */
        Raster raster;
        /**
        * Gleichungssystem Objekt.
        * @see Gleichungssystem
        */
        Gleichungssystem gs;
        /**
        * CFL Bedingung anwenden und neue Zeit berechnen.
        * @param n aktueller Zeitschritt.
        * @param time aktuelle Zeit.
        * @return neue Zeit.
        */
        double cflcon(int n, double time);
        /**
        * Abstrakte methode zur berechnung des Flusses der jeweiligen numerischen Methode.
        * @return Matrix der Flüsse (1D)
        */
        virtual std::vector< std::vector< std::vector< std::vector <double> > > >calc_method_flux(int dir) =0;
        /**
        * Aktualisiert alle zelle mithilfe des berechneten Flusses.
        */
        void update(std::vector< std::vector< std::vector< std::vector<double> > > >fi, int dir);
        /**
        * Schreibt Ergebnisse in Dateien u,d,ur,p
        */
        void write();
        /**
        * bestimmt die Elemente der Jacobi-Matrix in 3 Varianten, abhängig von der Wahl der EOS
        */
	void matrix_1d(double * values, int n, double uone, double utwo,
		       double uthree, double p, double done, double dtwo,
		       double ccl, double g, double ct, double cref, int variante);

        /**
        * bestimmt die Elemente der Jacobi-Matrix in 3 Varianten, abhängig von der Wahl der EOS
        */
	void matrix_2d(double * values_x, double * values_y, int n, double uone, double utwo,
		       double uthree, double ufour, double ufive, double p, double done, double dtwo,
		       double ccl, double g, double ct, double cref);
};

#endif // NUMERISCHE_METHODE_H
