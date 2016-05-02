#ifndef GLEICHUNGSSYSTEM_H
#define GLEICHUNGSSYSTEM_H
#include "exprtk.hpp"
#include "Konstanten.h"
#include "Zelle.h"
#include <iostream>
#include <vector>

/*!
* @class Gleichungssystem
* Die Klasse Gleichungssystem verwaltet die Formelvektoren u,f_u und s_u.
* Die Formelvektoren werden als Strings in Vektoren abgespeichert, und dann
* mittels Funktionen der exprtk Libiary interpretiert um mit ihnen zu rechnen.
*/

class Gleichungssystem
{
    public:
        /**
        * Konstruktor der klasse Gleichungssystem.
        * @param path Pfad zur Datei welche die Formeln enthält.
        * @param c Konstanten Objekte welches für berechnungen gebraucht wird.
        */
        Gleichungssystem(std::string path, Konstanten c);
        /**
        * Vektor von u.
        */
        std::vector<std::string> u;
        /**
        * Rückrechnung der Conserved Vatiables in die Physikalischen.
        */
        std::vector<std::string> uback;
        /**
        * Vektor der Formeln von f_u.
        */
        std::vector<std::string> f_u;
        /**
        * Vektor der Formeln von s_u.
        */
        std::vector<std::string> g_u;
        /**
        * Berechnet die Lösungen der Formel an der entsprechenden Position in u.
        * @param zelle welche die Werte zur Berechnung enthält.
        * @param pos Position der zu benützenden Formel im Formelvektor.
        * @return Lösung.
        */
        double solve_u(Zelle zelle, int pos);
        /**
        * Berechnet die Lösungen der Formel an der entsprechenden Position in f_u.
        * @param zelle welche die Werte zur Berechnung enthält.
        * @param pos Position der zu benützenden Formel im Formelvektor.
        * @return Lösung.
        */
        double solve_f_u(Zelle zelle, int pos);
        /**
        * Berechnet die Lösungen der Formel an der entsprechenden Position in s_u.
        * @param zelle welche die Werte zur Berechnung enthält.
        * @param pos Position der zu benützenden Formel im Formelvektor.
        * @return Lösung.
        */
        double solve_g_u(Zelle zelle, int pos);
	
	exprtk::symbol_table<double> symbol_table_u;
	exprtk::symbol_table<double> symbol_table_f_u;
	exprtk::symbol_table<double> symbol_table_g_u;
	exprtk::parser<double> parser;
        exprtk::expression<double> expression_u[5];
        exprtk::expression<double> expression_f_u[5];
        exprtk::expression<double> expression_g_u[5];
	
	double ux, uxr, uy, uyr, d, p;
	
    private:
        /**
        * Abgespeichertes Konstanten Objekt.
        * Wird nur einmal bei Konstruktor gesetzt, damit
        * es bei Rechnungen nicht immer neu übergeben werden muss.
        */
        Konstanten cons;
};

#endif // GLEICHUNGSSYSTEM_H
