#ifndef GRID_H
#define GRID_H
#include <iostream>
#include "constants.h"

/*!
 * @class Grid
 * Diese Klasse bietet ein Raster aus Zellen. Zellen werden in einem 1-Dimensionalen Array
 * abgespeichert, welches dann auch noch falls nötig 2-Dimensional oder 3-Dimensional
 * interpretiert wird.
 */

class Grid {

public:
	/**
	 * Konstruktor des Rasters.
	 * @param constants Pointer auf das Objekt, welches die Konstanten enthält.
	 * @param save_in Pfad zur Datei mit dem Speicherstand.
	 * @param choice Gibt die zu verwendende Initialisierungsmethode an
	 */
	//TODO: save_in entfernen
	Grid(Constants *constants, std::string save_in, int choice);

	/**
	 * Konstruktor für ein leeres eindimensionales Raster.
	 */
	Grid(int x);

	/**
	 * Konstruktor für ein leeres zweidimensionales Raster.
	 */
	Grid(int x, int y);

	/**
	 * Destruktor des Rasters.
	 * Gibt alokierten Speicher wieder frei.
	 */
	~Grid();

	/**
	 * Initialisierungsmethode für eine Verdünnungswelle auf einem 1D-Raster
	 */
	void init_1d_rarefaction();
	/**
	 * Initialisierungsmethode für eine Schockwelle auf einem 1D-Raster
	 */
	void init_1d_shockwave();
	/**
	 * Initialisierungsmethode für eine Verdünnungswelle in X-Richtung auf einem 2D-Raster
	 */
	void init_2d_rarefaction_0();
	/**
	 * Initialisierungsmethode für eine Verdünnungswelle in Y-Richtung auf einem 2D-Raster
	 */
	void init_2d_rarefaction_90();
	/**
	 * Initialisierungsmethode für eine Verdünnungswelle in einem 60 Grad Winkel auf einem 2D-Raster
	 */
	void init_2d_rarefaction_60();
	/**
	 * Initialisierungsmethode für eine Verdünnungswelle in einem 45 Grad Winkel auf einem 2D-Raster
	 */
	void init_2d_rarefaction_45();
	/**
	 * Initialisierungsmethode für eine Schockwelle-auf-Blase-Simulation auf einem 2D-Raster
	 */
	void init_2d_shockwave_bubble();
	/**
	 * Initialisierungsmethode mit einem Speicherstand auf einem 2D-Raster
	 */
	void init_2d_load_file(std::string save_in);

	/**
	 * Wendet Randbedingungen an.
	 * @param konstanten Konstanten wo festgelegt ist welche
	 * Randbedingungen angewendet werden.
	 */
	void apply_boundary_conditions();

	/**
	 * Wahl der Initialisierungsmethode des Rasters.
	 */
	//TODO: Remove
	int choice;

	/**
	 * Pointer auf ein Objekt mit allen Konstanten, welche für die Berechnungen benötigt werden
	 * @see Constants
	 */
	Constants *constants;

	/**
	 * Zellen werden in einem 2-dimensionalen Array abgespeichert.
	 * Die erste Dimension entspricht der Position und wird dann je nach Dimension des Rasters interpretiert.
	 * Die zweite Dimension entspricht den Zellwerten in der Reihenfolge: d, p, ux, uxr, uy, ...
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
	 * Anzahl der Zellenränder (+ 1) mit Grenze (ordnung * 2)
	 * [0] Breite
	 * [1] Weite
	 * [2] Tiefe
	 */
	int* grid_size_total;
	/**
	 * Dynamisch angelegtes Array je nach Dimension, dass die Art der Boundary Condition speichert
	 * [0] kleinster X-Wert, [1] größter X-Wert, [2] kleinster Y-Wert, ...
	 */
	int* boundary_conditions;
	/**
	 * Ordnung des Rasters
	 * Gibt an wie viele Nachbarwerte (pro Seite) miteinbezogen werden sollen
	 */
	int orderofgrid;
	/**
	 * Anzahl der Zellen
	 */
	int cellsgrid_size;
};

#endif // GRID_H
