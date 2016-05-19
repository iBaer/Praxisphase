#ifndef LAX_FRIEDRICH_H
#define LAX_FRIEDRICH_H

#include "numerische_methode.h"
#include "solver.h"

/*!
* @class LaxFriedrichMethod
* Klasse der Lax-Friedrichs Methode.
* Die Berechnung des Lax-Friedrich Flusses wird hier
* spezialisiert, rest wird von der Klasse numerische_methode geerbt.
*/


class Lax_Friedrich : public Solver
{
    public:
        /**
        * Konstruktor von der Klasse LaxFriedrichMethod.
        * Ruft einfach den Konstrukter von der geerbten Klasse auf.
        */
        Lax_Friedrich(Constants *constants, Computation *computation, Grid *grid);
        virtual ~Lax_Friedrich();

    protected:
        /**
        * Berechnung des Lax-Friedrich Flusses.
        * @return Das zurückgelieferte Objekt ist ein Vektor mit 4 Dimensionen 
	* (Formel, Raster x-Koordinate, Raster y-Koordinate, Flussrichtung)
        */
        double* calc_method_flux(double dt,int dir);

        //int size_total[0];
        //int size_total[1];
        //int size_m1[0];
        //int size_m1[1];
        int neqs;

        double *uall;
        double *fall;
        double *gall;

        double *** cs;
        double *** f;
        double *** g;

    	double** f_lax;

};

#endif // LAXFRIEDRICHMETHOD_H
