#ifndef GRID_H
#define GRID_H
#include <iostream>
#include "constants.h"

/*!
 * @class Raster
 * Diese Klasse bietet ein Raster aus Zellen. Zellen werden in einem 1-Dimensionalen Array
 * abgespeichert, welches dann auch noch falls nötig 2-Dimensional oder 3-Dimensional
 * interpretiert wird.
 */

class Grid {
	friend class Lax_Friedrich;
	friend class Force;
	friend class numerische_methode;
	friend class Computation;

public:
	/**
	 * Konstruktor des Rasters.
	 * @param const_in Pfad zur Datei, welche die initialisierungs Parameter enthält.
	 * @param function_in Pfad zur Datei, in der die Initierungsfunktion steht.
	 */
	Grid(Constants *constants, std::string save_in, int choice);
	/**
	 * Konstruktor für ein leeres eindimensionales Raster.
	 */
	Grid(int x);
	/**
	 * Konstruktor für ein leeres zweidimensionales Raster.
	 */
	Grid(int x, int y);


	void init_1d_rarefraction();
	void init_1d_shockwave();
	void init_2d_rarefraction_0();
	void init_2d_rarefraction_90();
	void init_2d_rarefraction_60();
	void init_2d_rarefraction_45();
	void init_2d_shockwave_bubble();
	void init_2d_load_file(std::string save_in);

	/**
	 * Destruktor des Rasters.
	 * Gibt alokierten Speicher wieder frei.
	 */
	~Grid();

	/**
	 * Wendet Randbedingungen an.
	 * @param konstanten Konstanten wo festgelegt ist welche
	 * Randbedingungen angewendet werden.
	 */
	void bcondi();

	/**
	 * Liefert Kopie der Zelle an gewünschter Position(1D).
	 * @param x Position in x-Richtung.
	 * @return Zellen Objekt.
	 */

	/**
	 * Initialisierung des Rasters.
	 */
	int choice;
	/**
	 * Konstanten Objekt welches für die berechnungen benötigt wird.
	 * @see Konstanten
	 */
	Constants *constants;

	/**
	 * Zellen werden 1-Dimensional im Raster abgespeichert.
	 * Array wird dann je nach Dimension des Rasters interpretiert.
	 *
	 */
	double** cellsgrid;
	/**
	 * Dimension des Rasters (Maximal 3).
	 */
	int dimension;
	/**
	 * Anzahl der Zellen ohne Grenze [0] height / [1] width / [2] depth
	 */
	int* grid_size;
	/**
	 * Anzahl der Zellen mit Grenze (ordnung * 2)[0] height / [1] width / [2] depth
	 */
	int* grid_size_total;
	int* boundary_conditions;
	int orderofgrid;
	int cellsgrid_size;
};

#endif // GRID_H
