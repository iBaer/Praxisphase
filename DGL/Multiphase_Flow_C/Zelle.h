#ifndef ZELLE_H
#define ZELLE_H
#include <iostream>

/*!
* @class Zelle
* Diese Klasse dient dazu sämtliche Werte einer Zelle in das 
* Objekt Zelle zusammenzufassen.
*/

class Zelle
{
 public:
  double ux; /**< Geschwindigkeit der Zelle in x-Richtung.*/
  double uy; /**< Geschwindigkeit der Zelle in y-Richtung.*/
  double uxr; /**< Relative Geschwindigkeit in der Zelle in x-Richtung.*/
  double uyr; /**< Relative Geschwindigkeit in der Zelle in y-Richtung.*/
  double d; /**< Dichte in der Zelle.*/
  double p; /**< Druck in der Zelle.*/
  
  /**
   * Konstruktor für eine Zelle.
   * Alle Werte der zelle werden zu beginn mit 0.0 initiert.
   */
  Zelle();
};

#endif // ZELLE_H
