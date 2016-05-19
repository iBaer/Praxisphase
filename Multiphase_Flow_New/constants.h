#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <string>

/*!
* @class Konstanten
* Die Klasse Konstanten dient hauptsächlich dazu sämtliche für
* Berechnungen benötigte Konstanten zusammenzufassen und
* leicht für Rechnungen verfügbar zu machen.
*/

class Constants
{
    private:
        /**
        * Konstruktor. Setzt die Konstanten.
        * @param input_const String des Pfades zu den festen Konstanten.
        * @param input_calc String des Pfades zu den einmalig zu berechnenden Konstanten.
        */
        Constants();

    public:

        static Constants& instance();
        /**
        * Beizeichnungen der Konstanten in einem Vektor.
        */
        //std::vector<std::string> const_name;
        /**
        * Vektor von Werten der Konstanten.
        */
        //std::vector<double> const_value;

    std::string input_const;
	int calceigv;
	int   variante;
	double teiler;
	int teilerend;
	double timeou;
	double cfl;
	int maxnt;
	int dimension;
	int order;
	double radius;
	int grid_size_x;
	int grid_size_y;
	double gamma;
	double pos_x_min;
	double pos_x_max;
	double pos_y_min;
	double pos_y_max;
	int bc_y_max;
	int bc_y_min;
	int bc_x_min;
	int bc_x_max;
	double cref;
	double done;
	double ccl;
	double rhol;
	double vl;
	double vrl;
	double vyl;
	double vyrl;
	double rhor;
	double vr;
	double vrr;
	double vyr;
	double vyrr;
    /**
    * K in den Formeln.
    */
    double ct;
};

#endif // CONSTANTS_H
