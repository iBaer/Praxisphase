#ifndef FORCE_H
#define FORCE_H

#include "numerische_methode.h"

/*!
* @class FORCE
* Klasse der FORCE Methode.
*/


class FORCE : public numerische_methode
{
    public:
        /**
        * Konstruktor der FORCE Methode.
        * Ruft den Konstruktor der geerbten Klasse auf.
        * @param const_in Dateiname wo Konstanten gespeichert sind.
        * @param formel_in Dateinamen-Kern für die Formeln.
        * @param winkel_in Dateiname wo Funktion der Geraden gespeichert ist.
        * @param kreis_in Dateiname wo Funktion des Kreises gespeichert ist.
        * @param save_in Dateiname wo für das Laden eine Speicherstands die Plots gespeichert sind.
        */
        FORCE(std::string const_in, std::string formel_in, std::string winkel_in, std::string kreis_in, std::string save_in);
    protected:
        /**
        * Berechnung des FORCE Flusses.
        * @return 4 Dimensionaler Vektor. Zusammenstellung: Gleichung, x-Position, y-Position , dimension
        */
        std::vector< std::vector< std::vector< std::vector <double> > > > calc_method_flux();
};

#endif // FORCE_H
