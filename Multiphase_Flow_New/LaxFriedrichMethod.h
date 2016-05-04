#ifndef LAXFRIEDRICHMETHOD_H
#define LAXFRIEDRICHMETHOD_H

#include "numerische_methode.h"
#include "solver.h"

/*!
* @class LaxFriedrichMethod
* Klasse der Lax-Friedrichs Methode.
* Die Berechnung des Lax-Friedrich Flusses wird hier
* spezialisiert, rest wird von der Klasse numerische_methode geerbt.
*/


class LaxFriedrichMethod : public Solver
{
    public:
        /**
        * Konstruktor von der Klasse LaxFriedrichMethod.
        * Ruft einfach den Konstrukter von der geerbten Klasse auf.
        */
        LaxFriedrichMethod(Constants *constants, Computation *computation, Grid *grid);
        ~LaxFriedrichMethod();

    protected:
        /**
        * Berechnung des Lax-Friedrich Flusses.
        * @return Das zurückgelieferte Objekt ist ein Vektor mit 4 Dimensionen 
	* (Formel, Raster x-Koordinate, Raster y-Koordinate, Flussrichtung)
        */
        double* calc_method_flux(double dt,int dir);

        int width;
        int height;
        int neqs;

        double *uall;
        double *fall;
        double *gall;

        double *** cs;
        double *** f;
        double *** g;

    	double* fiarray;


};

#endif // LAXFRIEDRICHMETHOD_H
