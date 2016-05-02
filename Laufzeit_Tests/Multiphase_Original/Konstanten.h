#ifndef KONSTANTEN_H
#define KONSTANTEN_H
#include <iostream>
#include "exprtk.hpp"
#include <vector>





/*!
* @class Konstanten
* Die Klasse Konstanten dient hauptsächlich dazu sämtliche für
* Berechnungen benötigte Konstanten zusammenzufassen und
* leicht für Rechnungen verfügbar zu machen.
*/

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
        /**
        * Registriert alle Konstanten in der übergebenen symbol_table und liefert sie danach wieder zurück.
        * @param in Eingangs symbol_table.
        * @return Ausgangs symbol_table.
        */
        exprtk::symbol_table<double> register_constants(exprtk::symbol_table<double> in);
        /**
        * Sucht nach einer Konstanten und liefert ihren Wert.
        */
        double search_con(std::string in);

};

#endif // KONSTANTEN_H
