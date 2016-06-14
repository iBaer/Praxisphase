#ifndef LAX_FRIEDRICH_H
#define LAX_FRIEDRICH_H

#include "numerical_method.h"
#include "solver.h"

/*!
* @class Lax_Friedrich
* Klasse der Lax-Friedrichs Methode.
* Die Berechnung des Lax-Friedrich Flusses wird hier
* spezialisiert, erbt von der Klasse Solver
*/


class Lax_Friedrich : public Solver
{
    public:
        /**
        * Konstruktor von der Klasse Lax_Friedrich.
        * Ruft den Konstruktor der Superklasse Solver auf.
        * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
        * @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird.
        * @param grid Pointer auf Raster-Objekt.
        */
        Lax_Friedrich(Constants *constants, Computation *computation, Grid *grid);

        /**
        * Destruktor
        */
        virtual ~Lax_Friedrich();

        /**
        * Untermethode, die die Berechnung des Lax-Friedrich-Flusses
        * und das Updaten der Zellen für eine Dimension durchführt.
        * @param dt Delta t.
        */
        void solve_1d(double dt);

        /**
        * Untermethode, die die Berechnung des Lax-Friedrich-Flusses
        * und das Updaten der Zellen für zwei Dimension durchführt.
        * @param dt Delta t.
        */
        void solve_2d_unsplit(double dt);

        /**
        * Untermethode, die die Berechnung des Lax-Friedrich-Flusses
        * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
        * @param dt Delta t.
        */
        void solve_2d_split_xtoy(double dt, int with_average, int rerun);
        void solve_2d_split_ytox(double dt, int with_average);


    protected:
        /**
        * Implementierung der virtuellen Methode zur Berechnung des Flußes.
        * Wählt lediglich das gewünschte Schema aus und delegiert die Berechnung weiter.
        * @param dt Delta t.
        * @param dir Unsplitting = 0, Splitting = 1.
        */
        //TODO: aus "dir" -> "splitting"
        void calc_method_flux(double dt,int split_method);

        //int size_total[0];
        //int size_total[1];
        //int size_m1[0];
        //int size_m1[1];

        /**
        * Anzahl der Gleichungen
        */
        int neqs;

        /**
        * 1D Array für den Fluss U, für schnelleren Zugriff
        */
        double *uall;

        /**
        * 1D Array für den Fluss F, für schnelleren Zugriff
        */
        double *fall;

        /**
        * 1D Array für den Fluss G, für schnelleren Zugriff
        */
        double *gall;

        /**
        * 3D Array für den Fluss U, für einfacherere Verwendung
        * Enthält Zeiger auf 1D Array
        * [0] = Gleichung, [1] = Zellen-Position in X, [2] = Zellen-Position in Y
        */
        //TODO: Dritte Dimension unbrauchbar, da dann auch Vierte gebraucht wird
        double *** cs;

        /**
        * 3D Array für den Fluss F, für einfacherere Verwendung
        * Enthält Zeiger auf 1D Array
        * [0] = Gleichung, [1] = Zellen-Position in X, [2] = Zellen-Position in Y
        */
        double *** f;

        /**
        * 3D Array für den Fluss G, für einfacherere Verwendung
        * Enthält Zeiger auf 1D Array
        * [0] = Gleichung, [1] = Zellen-Position in X, [2] = Zellen-Position in Y
        */
        double *** g;

        /**
        * 2D Array für den berechneten Lax-Friedrich-Fluss
        * [0] = Dimension, [1] = Zellen-Position, immer abstrahiert auf eine Dimension, bis zu 3D Richtung und Anzahl der Gleichungen
        */
    	double** f_lax;

    	Grid ** split_grid;

};

#endif // LAX_FRIEDRICH_H
