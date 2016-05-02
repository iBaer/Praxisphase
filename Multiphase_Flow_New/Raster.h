#ifndef RASTER_H
#define RASTER_H
#include "Zelle.h"
#include <iostream>

#include "constants.h"



/*!
* @class Raster
* Diese Klasse bietet ein Raster aus Zellen. Zellen werden in einem 1-Dimensionalen Array
* abgespeichert, welches dann auch noch falls nötig 2-Dimensional oder 3-Dimensional
* interpretiert wird.
*/

class Raster
{
  friend class LaxFriedrichMethod;
  friend class FORCE;
  friend class numerische_methode;
  friend class Computation;

    public:
        /**
        * Konstruktor des Rasters.
        * @param const_in Pfad zur Datei, welche die initialisierungs Parameter enthält.
        * @param function_in Pfad zur Datei, in der die Initierungsfunktion steht.
        */
        Raster(Constants *constants, std::string save_in);
        /**
        * Konstruktor für ein leeres eindimensionales Raster.
        */
        Raster(int x);
        /**
        * Konstruktor für ein leeres zweidimensionales Raster.
        */
        Raster(int x,int y);
        /**
        * Destruktor des Rasters.
        * Gibt alokierten Speicher wieder frei.
        */
        ~Raster();
        /**
        * Liefert Dimension des Rasters.
        */
        int getdim();
        /**
        * Wendet Randbedingungen an.
        * @param konstanten Konstanten wo festgelegt ist welche 
	* Randbedingungen angewendet werden.
        */
        void bcondi(int* CELLS , int ordnung);
        /**
        * Liefert die Breite des Rasters.
        */
        int getwidth();
        /**
        * Liefert die Höhe des Rasters.
        */
        int getheight();
        /**
        * Liefert Kopie der Zelle an gewünschter Position(1D).
        * @param x Position in x-Richtung.
        * @return Zellen Objekt.
        */
        Zelle get_Zelle(int x);
        /**
        * Liefert Kopie der Zelle an gewünschter Position(2D).
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        * @return Zellen Objekt.
        */
        Zelle get_Zelle(int x, int y);
        /**
        * Liefert Kopie der Zelle an gewünschter Position(3D).
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        * @param z Position in z-Richtung.
        * @return Zellen Objekt.
        */
        Zelle get_Zelle(int x,int y, int z);
        /**
        * Setzt die Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_ux(double in, int x);
        /**
        * Setzt die Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_ux(double in, int x, int y);
        /**
        * Setzt die relative Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_uxr(double in, int x);
        /**
        * Setzt die relative Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_uxr(double in, int x, int y);
                /**
        * Setzt die Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_uy(double in, int x);
        /**
        * Setzt die Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_uy(double in, int x, int y);
        /**
        * Setzt die relative Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_uyr(double in, int x);
        /**
        * Setzt die relative Geschwindigkeit einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_uyr(double in, int x, int y);
        /**
        * Setzt die Dichte einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_d(double in, int x);
        /**
        * Setzt die Dichte einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_d(double in, int x, int y);
        /**
        * Setzt den Druck einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        */
        void set_Zelle_p(double in, int x);
        /**
        * Setzt den Druck einer Zelle.
        * @param in Wert der gesetzt wird.
        * @param x Position in x-Richtung.
        * @param y Position in y-Richtung.
        */
        void set_Zelle_p(double in, int x, int y);
	/**
        * Initialisierung des Rasters.
        */
        int choice;
        /**
        * Konstanten Objekt welches für die berechnungen benötigt wird.
        * @see Konstanten
        */
        Constants *konstanten;
    private:
        /**
        * Zellen werden 1-Dimensional im Raster abgespeichert.
        * Array wird dann je nach Dimension des Rasters interpretiert.
        */
        Zelle* zelle;
        /**
        * Dimension des Rasters (Maximal 3).
        */
        int dimension;
        /**
        * Breite des Rasters.
        */
        int width;
        /**
        * Höhe des Rasters.
        */
        int height;
	/**
	 * Anzahl der Zellen
	 */
	int cells[2];
};

#endif // RASTER_H
