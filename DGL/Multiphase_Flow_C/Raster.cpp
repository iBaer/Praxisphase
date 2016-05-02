#include "Raster.h"
#include "Konstanten.h"
#include <math.h>
#include <string>
#include <fstream>

using namespace std;

#define KREIS(x,y)  sqrt((x-radius)*(x-radius)+y*y)
#define KREIS_F    "sqrt((x-radius)*(x-radius)+y*y)"

#define WINKEL_60(x,dx) 2*(x-0.01-dx)
#define WINKEL_60F     "2*(x-0.01-dx)"
#define ALPHA_60 60

#define WINKEL_45(x,dx) xpos-0.1*dx
#define WINKEL_45F     "xpos-0.1*dx"
#define ALPHA_45 45


/**
 *****************************************************************************************
 * Konstruktor für ein leeres eindimensionales Raster.
 *****************************************************************************************/
Raster::Raster(int x)
{
    width = x;
    height = 1;
    dimension = 1;
    zelle = new Zelle[x];
}

/**
 *****************************************************************************************
 * Konstruktor für ein leeres zweidimensionales Raster.
 *****************************************************************************************/
Raster::Raster(int x,int y)
{
    width = x;
    height = y;
    dimension = 2;
    zelle = new Zelle[x*y];
}

/**
 *****************************************************************************************
 * Konstruktor des Rasters.
 * @param const_in Pfad zur Datei, welche die initialisierungs Parameter enthält.
 * @param function_in Pfad zur Datei, in der die Initierungsfunktion steht.
 *****************************************************************************************/
Raster::Raster(Konstanten konstanten, string save_in)
{

  //Alle Konstanten holen, die gebraucht werden
  dimension = (int)konstanten.dimension;
  cells[0] = (int)konstanten.CELLSX;
  cells[1] = 1;
  if (dimension == 2)
    {
      cells[1] = (int)konstanten.CELLSY;
    }
  int ordnung = (int)konstanten.ordnung;

  width = cells[0]+2*ordnung+1;
  height = cells[1]+2*ordnung+1;
  choice = 0;
  zelle = new Zelle[width*height];
    
  double mor = konstanten.mor;
  double mol = konstanten.mol;
  double mur = konstanten.mur; //2D
  double mul = konstanten.mul; //2D
  
  double dx = (mor-mol)/(double)cells[0];
  
  double xpos = 0.0;
  double ypos = 0.0;  //2D
  
  double rhol = konstanten.rhol;
  double vl = konstanten.vl;
  double vrl = konstanten.vrl;
  double vyl = konstanten.vyl;
  double vyrl = konstanten.vyrl;
  double rhor = konstanten.rhor;
  double vr = konstanten.vr;
  double vrr = konstanten.vrr;
  double vyr = konstanten.vyr;
  double vyrr = konstanten.vyrr;
  
  switch(dimension)
    {
      // 1-dimensionale Fall
    case(1):
      cout << dimension << "-dimensionales Gitter" << endl;
      cout << "Raster Initiierung: - 0 Rarefraction wave - 1 one shock wave !" << endl;
      cout << "Wahl:";
      cin >> choice;

      switch(choice)
        {
        case(0):
	  {
            //1D Riemann-Problem
            for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
                xpos = mol + (mor-mol)*((double)(n-ordnung)/cells[0]);
                if(xpos <= 0.0)
		  {
                    set_Zelle_d(rhol,n);
                    set_Zelle_ux(vl,n);
                    set_Zelle_uxr(vrl,n);
		  }
                else
		  {
                    set_Zelle_d(rhor,n);
                    set_Zelle_ux(vr,n);
                    set_Zelle_uxr(vrr,n);
		  }
	      }
            break;
	  }
	  //1D Schockwelle
        case(1):
	  {
            for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
                xpos = mol + (mor-mol)*((double)(n-ordnung)/cells[0]);
                if(xpos == 0.0)
		  {
                    set_Zelle_d(1323,n);
                    set_Zelle_ux(100000,n);
                    set_Zelle_uxr(0,n);
		  }
                else
		  {
                    set_Zelle_d(1000,n);
                    set_Zelle_ux(0,n);
                    set_Zelle_uxr(0,n);
		  }
	      }
            break;
	  }
        }
      break;
      
      // 2-dimensionale Fall
    case(2):

      cout << dimension << "-dimensionales Gitter" << endl;
      cout << "Raster Initiierung:"  << endl << "0: Rarefraction wave" << endl
	   << "1: 60 Grad Winkel" << endl << "2: 45 Grad Winkel" << endl
	   << "3: Blasen Simulation" << endl << "4: Speicherstand laden!"
	   << endl;
      cout << "Wahl:";
      cin >> choice;
      
      cout << "Simulationsfläche: ["<< mol << "," << mul << "] bis ["
	   << mor  << "," << mur << "]" << endl;
      cout << "Simulationsgitter: ["<< cells[0] << "," << cells[1] << "]" << endl;
      
      
      switch(choice)
        {
	  
	  //2D Riemann-Problem: 0 Rarefraction wave in x or y direction
	case(0):
	  {
	    cout << "Gerade parallel zur y-Achse" << endl;
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    xpos = mol + (mor-mol)*(((double)n-ordnung)/(double)cells[0]);
		    ypos = mul + (mur-mul)*(((double)m-ordnung)/(double)cells[1]);
		    if(xpos <= 0.0)
		      {
			this->set_Zelle_d(rhol,n,m);
			this->set_Zelle_ux(vl,n,m);
			this->set_Zelle_uxr(vrl,n,m);
			this->set_Zelle_uy(vyl,n,m);
			this->set_Zelle_uyr(vyrl,n,m);
		      }
		    else
		      {
			this->set_Zelle_d(rhor,n,m);
			this->set_Zelle_ux(vr,n,m);
			this->set_Zelle_uxr(vrr,n,m);
			this->set_Zelle_uy(vyr,n,m);
			this->set_Zelle_uyr(vyrr,n,m);
		      }
		  }
	      }
	    break;
	  }

	  //2D Riemann-Problem mit Winkel, Rarefraction wave in 60 degrees direction
        case(1):
	  {
	    cout << "Gerade mit der Gleichung ypos=" << WINKEL_60F << " als Grenze " << endl;
	      
	    double alpha=ALPHA_60*M_PI/180.0;
	    double vxl_winkel = - fabs(vl)*sin(alpha);
	    double vyl_winkel = fabs(vl)*cos(alpha);	    
	    double vxr_winkel = fabs(vr)*sin(alpha);
	    double vyr_winkel = -fabs(vr)*cos(alpha);

	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    xpos = mol + (mor - mol)*(((double)n-ordnung)/(double)cells[0]);
		    ypos = mul + (mur - mul)*(((double)m-ordnung)/(double)cells[1]);
		    if(ypos >= WINKEL_60(xpos,dx))
		      {
			this->set_Zelle_d(rhol,n,m);
			this->set_Zelle_ux(vxl_winkel,n,m);
			this->set_Zelle_uy(vyl_winkel,n,m);

			this->set_Zelle_uxr(vrl,n,m);
			this->set_Zelle_uyr(vyrl,n,m);
		      }
		    else
		      {
			this->set_Zelle_d(rhor,n,m);
			this->set_Zelle_ux(vxr_winkel,n,m);
			this->set_Zelle_uy(vyr_winkel,n,m);

			this->set_Zelle_uxr(vrr,n,m);
			this->set_Zelle_uyr(vyrr,n,m);
		      }
		  }
	      }
	    break;
	  }
	  
	  //2D Riemann-Problem mit Winkel, Rarefraction wave in 45 degrees direction
        case(2):
	  {
	    cout << "Gerade mit der Gleichung ypos=" << WINKEL_45F
		 << " als Grenze " << endl;
	      
	    double alpha=ALPHA_45*M_PI/180.0;
	    double vxl_winkel = - fabs(vl)*sin(alpha);
	    double vyl_winkel = fabs(vl)*cos(alpha);
	    double vxr_winkel = fabs(vr)*sin(alpha);
	    double vyr_winkel = -fabs(vr)*cos(alpha);
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    xpos = mol + (mor - mol)*(((double)n-ordnung)/(double)cells[0]);
		    ypos = mul + (mur - mul)*(((double)m-ordnung)/(double)cells[1]); 
		    if(ypos >=  WINKEL_45(xpos,dx))
		      {
			this->set_Zelle_d(rhol,n,m);
			this->set_Zelle_ux(vxl_winkel,n,m);
			this->set_Zelle_uy(vyl_winkel,n,m);
			this->set_Zelle_uxr(vrl,n,m);
			this->set_Zelle_uyr(vyrl,n,m);
		      }
		    else
		      {
			this->set_Zelle_d(rhor,n,m);
			this->set_Zelle_ux(vxr_winkel,n,m);
			this->set_Zelle_uy(vyr_winkel,n,m);
			this->set_Zelle_uxr(vrr,n,m);
			this->set_Zelle_uyr(vyrr,n,m);
		      }
		  }
	      }
	    break;
	  }
	  
	  //Schockwelle-auf-Blase-Simulation
	case(3):
	  {
	    
	    double radius = konstanten.radius;
	    
	    cout << "Kreis mit Radius " << radius << " und Gleichung " << KREIS_F << endl;   
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    xpos = mol + (mor - mol)*(((double)n-ordnung)/(double)cells[0]);
		    ypos = mul + (mur - mul)*(((double)m-ordnung)/(double)cells[1]); 

		    if(radius >= KREIS(xpos,ypos))
		      {
			set_Zelle_d(8,n,m);
		      }
		    else
		      {
			set_Zelle_d(10,n,m);
		      }
		    set_Zelle_ux(0,n,m);
		    set_Zelle_uxr(0,n,m);
		    set_Zelle_uy(0,n,m);
		    set_Zelle_uyr(0,n,m);
		  }
	      }
	    
	    //int shockwave_pos = mor + 1.8/(mor - mol);
	    //int shockwave_zell = 1.8/(mor - mol)*cells[0];
	    int shockwave_zell = 0.5*(mor - mol)/(mor - mol)*cells[0];
	    double shockwave_pos = mol + shockwave_zell*dx;
	    cout << "Schockwelle bei " << shockwave_pos << ", x-Zelle " 
		 << shockwave_zell << endl;
	    for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
	      {
		set_Zelle_d(20,shockwave_zell,m);
		set_Zelle_ux(10,shockwave_zell,m);
	      }
	    break;
	  }
	  
	  // Konfiguration laden
        case(4):
	  {
	    //Speicherstand laden
	    ifstream save_input(save_in.c_str());
	    ifstream save_input_line;
	    string line_open;
	    string line;
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_d(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_ux(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_uy(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_uxr(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_uyr(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    getline(save_input,line_open);
	    save_input_line.open(line_open.c_str());
	    
	    for(int n = ordnung ; n < cells[0]+ordnung+1 ; n++)
	      {
		for(int m = ordnung ; m < cells[1]+ordnung+1 ; m++)
		  {
		    getline(save_input_line,line);
		    line = line.substr(line.find("\t")+1);
		    set_Zelle_p(atof(line.c_str()),n,m);
		  }
	      }
	    save_input_line.close();
	    
	    break;
	  }
        }
      break;
    }
}

/**
 *****************************************************************************************
 * Destruktor
 *****************************************************************************************/
Raster::~Raster()
{
    delete[] zelle;
}

/**
 *****************************************************************************************
 * boundary condition
 *****************************************************************************************/
void Raster::bcondi(Konstanten konstanten, int* CELLS , int ordnung)
{
    double upbc = konstanten.upbc;
    double downbc = konstanten.downbc;
    double rightbc = konstanten.rightbc;
    double leftbc = konstanten.leftbc;
    //double xpos = 0.0 , ypos = 0.0;

    // ACHTUNG, FUER 2. ORDNUNG NOCHMAL ALLE GRENZEN KONTROLLIEREN

    switch(dimension){
    case(1):
        if(leftbc == 0)
        {
            set_Zelle_d( get_Zelle(ordnung).d ,0);
            set_Zelle_ux( get_Zelle(ordnung).ux ,0);
            set_Zelle_uxr( get_Zelle(ordnung).uxr ,0);
        }
        else
        {
            set_Zelle_d( get_Zelle(ordnung).d ,0);
            set_Zelle_ux( -get_Zelle(ordnung).ux ,0); // 0.0 - get...
            set_Zelle_uxr( get_Zelle(ordnung).uxr ,0);
        }
        if(rightbc == 0)
        {
            set_Zelle_d( get_Zelle(CELLS[0]+ordnung).d ,CELLS[0]+ordnung+1);
            set_Zelle_ux( get_Zelle(CELLS[0]+ordnung).ux ,CELLS[0]+ordnung+1);
            set_Zelle_uxr( get_Zelle(CELLS[0]+ordnung).uxr ,CELLS[0]+ordnung+1);
        }
        else
        {
            set_Zelle_d( get_Zelle(CELLS[0]+ordnung).d ,CELLS[0]+ordnung+1);
            set_Zelle_ux( 0.0-get_Zelle(CELLS[0]+ordnung).ux ,CELLS[0]+ordnung+1);
            set_Zelle_uxr( get_Zelle(CELLS[0]+ordnung).uxr , CELLS[0]+ordnung+1);
        }
        break;
    case(2):
        //2D
        for(int n = 0 ; n < 2 ; n++)
        {
            for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
            {
                for(int y = 0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                {
                    if(x==0)
                    {
                        if(leftbc==0)
                        {
                            set_Zelle_d( get_Zelle(ordnung,y).d,0,y);
                            set_Zelle_ux( get_Zelle(ordnung,y).ux ,0,y);
                            set_Zelle_uxr( get_Zelle(ordnung,y).uxr ,0,y);
                            set_Zelle_uy( get_Zelle(ordnung,y).uy ,0,y);
                            set_Zelle_uyr( get_Zelle(ordnung,y).uyr ,0,y);
                        }
                        else
                        {
                            set_Zelle_d( get_Zelle(ordnung,y).d ,0,y);
                            set_Zelle_ux( 0.0-get_Zelle(ordnung,y).ux ,0,y);
                            set_Zelle_uxr( get_Zelle(ordnung,y).uxr ,0,y);
                            set_Zelle_uy( 0.0-get_Zelle(ordnung,y).uy ,0,y);
                            set_Zelle_uyr( get_Zelle(ordnung,y).uyr ,0,y);
                        }

                    }
                    if(x==(CELLS[0]+ordnung+1))
                    {
                        if(rightbc==0)
                        {
                            set_Zelle_d( get_Zelle(CELLS[0]+ordnung,y).d ,CELLS[0]+ordnung+1,y);
                            set_Zelle_ux( get_Zelle(CELLS[0]+ordnung,y).ux ,CELLS[0]+ordnung+1,y);
                            set_Zelle_uxr( get_Zelle(CELLS[0]+ordnung,y).uxr ,CELLS[0]+ordnung+1,y);
                            set_Zelle_uy( get_Zelle(CELLS[0]+ordnung,y).uy ,CELLS[0]+ordnung+1,y);
                            set_Zelle_uyr( get_Zelle(CELLS[0]+ordnung,y).uyr ,CELLS[0]+ordnung+1,y);
                        }
                        else
                        {
                            set_Zelle_d( get_Zelle(CELLS[0]+ordnung,y).d ,CELLS[0]+ordnung+1,y);
                            set_Zelle_ux( 0.0-get_Zelle(CELLS[0]+ordnung,y).ux ,CELLS[0]+ordnung+1,y);
                            set_Zelle_uxr( get_Zelle(CELLS[0]+ordnung,y).uxr , CELLS[0]+ordnung+1,y);
                            set_Zelle_uy( 0.0-get_Zelle(CELLS[0]+ordnung,y).uy ,CELLS[0]+ordnung+1,y);
                            set_Zelle_uyr( get_Zelle(CELLS[0]+ordnung,y).uyr , CELLS[0]+ordnung+1,y);
                        }
                    }
                    if(y==0)
                    {
                        if(downbc==0)
                        {
                            set_Zelle_d( get_Zelle(x,ordnung).d,x,0);
                            set_Zelle_uy( get_Zelle(x,ordnung).uy ,x,0);
                            set_Zelle_uyr( get_Zelle(x,ordnung).uyr ,x,0);
                            set_Zelle_ux( get_Zelle(x,ordnung).ux ,x,0);
                            set_Zelle_uxr( get_Zelle(x,ordnung).uxr ,x,0);
                        }
                        else
                        {
                            set_Zelle_d( get_Zelle(x,ordnung).d,x,0);
                            set_Zelle_uy( 0.0-get_Zelle(x,ordnung).uy ,x,0);
                            set_Zelle_uyr( get_Zelle(x,ordnung).uyr ,x,0);
                            set_Zelle_ux( 0.0-get_Zelle(x,ordnung).ux ,x,0);
                            set_Zelle_uxr( get_Zelle(x,ordnung).uxr ,x,0);
                        }
                    }
                    if(y==(CELLS[1]+ordnung+1))
                    {
                        if(upbc==0)
                        {
                            set_Zelle_d( get_Zelle(x,CELLS[1]+ordnung).d,x,CELLS[1]+ordnung+1);
                            set_Zelle_uy( get_Zelle(x,CELLS[1]+ordnung).uy ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uyr(get_Zelle(x,CELLS[1]+ordnung).uyr ,x,CELLS[1]+ordnung+1);
                            //set_Zelle_ux( get_Zelle(x,CELLS[1]+ordnung).ux ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uxr( get_Zelle(x,CELLS[1]+ordnung).uxr ,x,CELLS[1]+ordnung+1);

			    zelle[x + (CELLS[1]+ordnung+1)*width].ux = zelle[x + (CELLS[1]+ordnung)*width].ux;

			}
                        else
                        {
                            set_Zelle_d( get_Zelle(x,CELLS[1]+ordnung).d,x,CELLS[1]+ordnung+1);
                            set_Zelle_uy( 0.0-get_Zelle(x,CELLS[1]+ordnung).uy ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uyr( get_Zelle(x,CELLS[1]+ordnung).uyr ,x,CELLS[1]+ordnung+1);
			    set_Zelle_ux( 0.0-get_Zelle(x,CELLS[1]+ordnung).ux   ,  x,  CELLS[1]+ordnung+1);
                            set_Zelle_uxr( get_Zelle(x,CELLS[1]+ordnung).uxr ,x,CELLS[1]+ordnung+1);
                        }
                    }

                }
            }
        }


        break;
    }
}

/**
 *****************************************************************************************
 * überflüssige Get- und Setmethoden
 *****************************************************************************************/

int Raster::getwidth()
{
    return width;
}

int Raster::getheight()
{
    return height;
}

int Raster::getdim()
{
    return dimension;
}

Zelle Raster::get_Zelle(int x)
{
    return get_Zelle(x,0,0);
}

Zelle Raster::get_Zelle(int x , int y)
{
    return get_Zelle(x,y,0);
}

Zelle Raster::get_Zelle(int x, int y , int z)
{
    return  zelle[x + y * width + z * width * height];
}

void Raster::set_Zelle_ux(double in, int x)
{
    zelle[x].ux = in;
}

void Raster::set_Zelle_ux(double in, int x, int y)
{
    zelle[x+y*width].ux = in;
}

void Raster::set_Zelle_uxr(double in, int x)
{
    zelle[x].uxr = in;
}

void Raster::set_Zelle_uxr(double in, int x, int y)
{
    zelle[x+y*width].uxr = in;
}

void Raster::set_Zelle_uy(double in, int x)
{
    zelle[x].uy = in;
}

void Raster::set_Zelle_uy(double in, int x, int y)
{
    zelle[x+y*width].uy = in;
}

void Raster::set_Zelle_uyr(double in, int x)
{
    zelle[x].uyr = in;
}

void Raster::set_Zelle_uyr(double in, int x, int y)
{
    zelle[x+y*width].uyr = in;
}

void Raster::set_Zelle_d(double in, int x)
{
    zelle[x].d = in;
}

void Raster::set_Zelle_d(double in, int x, int y)
{
    zelle[x+y*width].d = in;
}

void Raster::set_Zelle_p(double in, int x)
{
    zelle[x].p = in;
}

void Raster::set_Zelle_p(double in, int x, int y)
{
    zelle[x+y*width].p = in;
}
