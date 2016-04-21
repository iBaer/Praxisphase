#ifndef LAXFRIEDRICHMETHOD_H
#define LAXFRIEDRICHMETHOD_H

#include "numerische_methode.h"

/*!
* @class LaxFriedrichMethod
* Klasse der Lax-Friedrichs Methode.
* Die Berechnung des Lax-Friedrich Flusses wird hier
* spezialisiert, rest wird von der Klasse numerische_methode geerbt.
*/


class LaxFriedrichMethod : public numerische_methode
{
    public:
        /**
        * Konstruktor von der Klasse LaxFriedrichMethod.
        * Ruft einfach den Konstrukter von der geerbten Klasse auf.
        */
        LaxFriedrichMethod(std::string const_in, std::string formel_in, std::string save_in);
    protected:
        /**
        * Berechnung des Lax-Friedrich Flusses.
        * @return Das zurückgelieferte Objekt ist ein Vektor mit 4 Dimensionen 
	* (Formel, Raster x-Koordinate, Raster y-Koordinate, Flussrichtung)
        */
        std::vector< std::vector< std::vector< std::vector <double> > > > calc_method_flux(int dir);
};

#endif // LAXFRIEDRICHMETHOD_H
