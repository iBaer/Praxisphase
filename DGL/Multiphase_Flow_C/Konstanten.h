#ifndef KONSTANTEN_H
#define KONSTANTEN_H
#include <iostream>
#include <math.h>
#include <vector>

/*!
* @class Konstanten
* Die Klasse Konstanten dient hauptsächlich dazu sämtliche für
* Berechnungen benötigte Konstanten zusammenzufassen und
* leicht für Rechnungen verfügbar zu machen.
*/V

class Konstanten
{
    public:
        /**
        * Konstruktor. Setzt die Konstanten.
        * @param input_const String des Pfades zu den festen Konstanten.
        * @param input_calc String des Pfades zu den einmalig zu berechnenden Konstanten.
        */
        Konstanten(std::string input_const);
        /**
        * Beizeichnungen der Konstanten in einem Vektor.
        */
        std::vector<std::string> const_name;
        /**
        * Vektor von Werten der Konstanten.
        */
        std::vector<double> const_value;

	int calceigv;
	int   variante;
	double teiler;
	int teilerend;
	double timeou;
	double cfl;
	int maxnt;
	int dimension;
	int ordnung;
	double radius;
	int CELLSX;
	int CELLSY;
	double g;
	double mol;
	double mor;
	double mul;
	double mur;
	int upbc;
	int downbc;
	int leftbc;
	int rightbc;
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
};
#endif // KONSTANTEN_H
