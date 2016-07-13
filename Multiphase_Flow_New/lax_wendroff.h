#ifndef LAX_WENDROFF_H
#define LAX_WENDROFF_H

#include "numerical_method.h"
#include "solver.h"

/*!
* @class Lax_Wendroff
* Klasse der Lax-Wendroff Methode.
* Die Berechnung des Lax-Wendroff Flusses wird hier
* spezialisiert, erbt von der Klasse Solver
*/


class Lax_Wendroff : public Solver
{
    public:
        /**
        * Konstruktor von der Klasse Lax_Wendroff.
        * Ruft den Konstruktor der Superklasse Solver auf.
        * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
        * @param computation Pointer auf das Objekt, das für die Berechnungen benötigt wird.
        * @param grid Pointer auf Raster-Objekt.
        */
        Lax_Wendroff(Constants *constants, Computation *computation, Grid *grid);

        /**
        * Destruktor
        */
        virtual ~Lax_Wendroff();

        /**
        * Untermethode, die die Berechnung des Lax-Wendroff-Flusses
        * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
        * @param dt Delta t.
        */
        void solve_1d(double dt);

        /**
        * Untermethode, die die Berechnung des Lax-Wendroff-Flusses
        * und das Updaten der Zellen für zwei Dimensionen mit dem konventionellen Splitting-Schema durchführt.
        * @param dt Delta t.
        */
        void solve_2d_unsplit(double dt);

    protected:
        /**
        * Implementierung der virtuellen Methode zur Berechnung des Flußes.
        * Wählt lediglich das gewünschte Schema aus und delegiert die Berechnung weiter.
        * @param dt Delta t.
        * @param dir Unsplitting = 0, Splitting = 1.
        */
        //TODO: aus "dir" -> "splitting"
        void calc_method_flux(double dt, Grid * grid);
    	void allocate_cache(Grid * grid);
    	void delete_cache();

        //int size_total[0];
        //int size_total[1];
        //int size_m1[0];
        //int size_m1[1];

	  double *uall;
	  double *fall;
	  double *gall;
	  double *f_laxall;
	  double *f_rieall;
	  double *g_laxall;
	  double *g_rieall;

	  double *** cs;
	  double *** fd;
	  double *** gd;
	  double *** f_lax_half;
	  double *** f_rie;
	  double *** g_lax_half;
	  double *** g_rie;

	  double* f_laxwen;

};

#endif // LAX_WENDROFF_H
