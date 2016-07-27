#ifndef FORCE_H
#define FORCE_H

#include "numerical_method.h"
#include "solver.h"
/*!
* @class FORCE
* Klasse der FORCE Methode.
*/


class Force : public Solver
{
    public:
        /**
        * Konstruktor der FORCE Methode.
        * Ruft den Konstruktor der geerbten Klasse auf.
        * @param const_in Dateiname wo Konstanten gespeichert sind.
        * @param formel_in Dateinamen-Kern für die Formeln.
        * @param save_in Dateiname wo für das Laden eine Speicherstands die Plots gespeichert sind.
        */
        Force(Constants *constants, Computation *computation, Grid *grid);
        virtual ~Force();

    protected:
        /**
        * Berechnung des FORCE Flusses.
        * @return 4 Dimensionaler Vektor. Zusammenstellung: Gleichung, x-Position, y-Position , dimension
        */
        void calc_method_flux(double dt, Grid * grid);
        void allocate_cache(Grid * grid);
        void delete_cache();
        void get_2d_half_flux(double*** half_flux, int& dimension, int& size);

  	  //int size_total[0];
  	  //int size_total[1];

  	  double *uall;
  	  double *fall;
  	  double *gall;
  	  double *f_laxall;
  	  double *f_rieall;
  	  double *g_laxall;
  	  double *g_rieall;

  	  double *** cs;
  	  double *** fd;
  	  double *** gd;
  	  double *** f_lax;
  	  double *** f_rie;
  	  double *** g_lax;
  	  double *** g_rie;

  	  double** f_force;

      void solve_1d(double dt);
      void solve_2d_unsplit(double dt);
      void solve_2d_split(double dt);
};

#endif // FORCE_H
