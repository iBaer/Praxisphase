#include <fstream>

#include "computation.h"
using namespace std;

Computation::Computation(Constants *constants)
{
	c = constants;

  cref = c->cref;
  done = c->done;
  ccl = c->ccl;
  gc = c->g;

  ccl12 = 1.0-2.0*ccl;
  ccl12h = 0.5*(1.0-2.0*ccl);
  gcinv = 1.0/gc;
  gcdgcm1 = gc/(gc-1.0);
  gcm1dgc = (gc-1.0)/gc;
  cclm1 = ccl*(1.0-ccl);
  powcref = pow(cref,gcinv);  

  int dim = c->dimension;

  if(dim == 1)
    neqs = 3;
  if(dim == 2)
    neqs = 5;

}

Computation& Computation::instance(Constants *constants) {
	static Computation computation(constants);
	return computation;
}

/**
 * Berechnet die Werte U in 1-D
 * @param u, 3-d Feld für die U-Werte
 * @param raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_u_1d(double *** u, Grid * raster, int * cells, int ordnung)
{
  for(int i = 0 ; i < cells[0]+2*ordnung+1 ; i++)  
    {
      u[0][i][0] = raster->cellsgrid[i][0];
      u[1][i][0] = raster->cellsgrid[i][0]*raster->cellsgrid[i][2];
      u[2][i][0] = raster->cellsgrid[i][3];
    }
}

/**
 * Berechnet die Lösungen der Formel für den Fluss F in 1-D
 * @param f, 3-d Feld für die F-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_f_1d(double *** f, Grid * raster, int * cells, int ordnung)
{
  double d,ux,uxr,p;
  
  for(int i = 0 ; i < cells[0]+2*ordnung+1 ; i++)  
    {
      /*d = raster->zelle[i].d;
      ux = raster->zelle[i].ux;
      uxr = raster->zelle[i].uxr;
      p = raster->zelle[i].p;*/

      d = raster->cellsgrid[i][0];
      ux = raster->cellsgrid[i][2];
      uxr = raster->cellsgrid[i][3];
      p = raster->cellsgrid[i][1];
	  
      f[0][i][0] = d*ux;
      f[1][i][0] = d*ux*ux+d*cclm1*uxr*uxr+p;
      f[2][i][0] = ux*uxr+ccl12h*uxr*uxr+powcref*gcdgcm1*pow(p,gcm1dgc)-(p/done);	    
    }
}

/**
 * Berechnet die Werte U in 2-D
 * @param u, 3-d Feld für die U-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_u_2d(double *** u, Grid * raster, int * cells, int ordnung)
{
  int pos;
  double d,ux,uy,uxr,uyr;

  for(int x = 0 ; x < cells[0]+2*ordnung+1 ; x++)
    {
      for(int y = 0 ; y < cells[1]+2*ordnung+1 ; y++)
	{
	  pos = x+y*raster->getwidth();
	  d = raster->cellsgrid[pos][0];
	  ux = raster->cellsgrid[pos][2];
	  uy = raster->cellsgrid[pos][4];
	  uxr = raster->cellsgrid[pos][3];
	  uyr = raster->cellsgrid[pos][5];

	  u[0][x][y] = d;
	  u[1][x][y] = d*ux;
	  u[2][x][y] = d*uy;
	  u[3][x][y] = uxr;
	  u[4][x][y] = uyr;
	}
    }
}

/**
 * Berechnet die Lösungen der Formel für den Fluss F in 2-D
 * @param f, 3-d Feld für die F-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_f_2d(double *** f, Grid * raster, int * cells, int ordnung)
{
  int pos;
  double d,ux,uy,uxr,uyr,p;

  for(int x = 0 ; x < cells[0]+2*ordnung+1 ; x++)
    {
      for(int y = 0 ; y < cells[1]+2*ordnung+1 ; y++)
	{
	  pos = x+y*raster->getwidth();
	  d = raster->zelle[pos].d;
	  ux = raster->zelle[pos].ux;
	  uy = raster->zelle[pos].uy;
	  uxr = raster->zelle[pos].uxr;
	  uyr = raster->zelle[pos].uyr;
	  p = raster->zelle[pos].p;

		
	  f[0][x][y]=d*ux;
	  f[1][x][y]=d*ux*ux+d*cclm1*uxr*uxr+p;
	  f[2][x][y]=d*ux*uy+d*cclm1*uxr*uyr;
	  f[3][x][y]=ux*uxr+ccl12h*uxr*uxr+
	    powcref*gcdgcm1*pow(p,gcm1dgc)-(p/done);
	  f[4][x][y]=ux*uyr+uy*uxr+ccl12*uxr*uyr;
	}
    }
}

/**
 * Berechnet die Lösungen der Formel für den Fluss G in 2-D
 * @param g, 3-d Feld für die G-Werte
 * @param Raster, welche die Werte zur Berechnung enthält.
 * @param cells, Ausdehnung des Gitters
 * @param ordnung des Algorithmus, wichtig für die Ausdehnung des Gitters
 */
void Computation::compute_g_2d(double *** g, Grid * raster, int * cells, int ordnung)
{
  int pos;
  double d,ux,uy,uxr,uyr,p;

  for(int x = 0 ; x < cells[0]+2*ordnung+1 ; x++)
    {
      for(int y = 0 ; y < cells[1]+2*ordnung+1 ; y++)
	{
	  pos = x+y*raster->getwidth();
	  d = raster->zelle[pos].d;
	  ux = raster->zelle[pos].ux;
	  uy = raster->zelle[pos].uy;
	  uxr = raster->zelle[pos].uxr;
	  uyr = raster->zelle[pos].uyr;
	  p = raster->zelle[pos].p;


	  g[0][x][y] = d*uy;
	  g[1][x][y] = d*ux*uy+d*cclm1*uxr*uyr;
	  g[2][x][y] = d*uy*uy+d*cclm1*uyr*uyr+p;
	  g[3][x][y] = ux*uyr+uy*uxr+ccl12*uxr*uyr;
	  g[4][x][y] = uy*uyr+ccl12h*uyr*uyr+
	    powcref*gcdgcm1*pow(p,gcm1dgc)-(p/done);
	}
    }
}
