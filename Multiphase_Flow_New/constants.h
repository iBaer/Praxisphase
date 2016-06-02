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
	int cfl_with_eig;
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
	double rho_one;
	double ccl;
	double rho_left;
	double v_left;
	double v_r_left;
	double v_y_left;
	double v_r_y_left;
	double rho_right;
	double v_right;
	double v_r_right;
	double v_y_right;
	double v_r_y_right;
    /**
    * K in den Formeln.
    */
    double ct;
};

#endif // CONSTANTS_H
