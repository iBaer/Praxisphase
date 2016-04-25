#include "Raster.h"
#include <math.h>
#include <string>
#include <fstream>
#include "exprtk.hpp"

using namespace std;
using namespace exprtk;

Raster::Raster(int x)
{
    width = x;
    height = 1;
    dimension = 1;
    zelle = new Zelle[x];
}

Raster::Raster(int x,int y)
{
    width = x;
    height = y;
    dimension = 2;
    zelle = new Zelle[x*y];
}

Raster::Raster(Konstanten konstanten, string winkel_in, string kreis_in, string save_in)
{

    //Alle Konstanten holen, die gebraucht werden
    dimension = (int)konstanten.search_con("dimension");
    int ordnung = (int)konstanten.search_con("ordnung");
    int CELLSX = (int)konstanten.search_con("CELLSX");
    int CELLSY = (int)konstanten.search_con("CELLSY");

    int choice = 0;

    width = CELLSX+2*ordnung+1;
    height = CELLSY+2*ordnung+1;
    zelle = new Zelle[width*height];

    double mor = konstanten.search_con("mor");
    double mol = konstanten.search_con("mol");
    double mur = konstanten.search_con("mur"); //2D
    double mul = konstanten.search_con("mul"); //2D

    double dx = (mor-mol)/(double)CELLSX;
    double dy = (mor-mur)/(double)CELLSY;

    double xpos = 0.0;
    double ypos = 0.0;  //2D

    double rhol = konstanten.search_con("rhol");
    double vl = konstanten.search_con("vl");
    double vrl = konstanten.search_con("vrl");
    double vyl = konstanten.search_con("vyl");
    double vyrl = konstanten.search_con("vyrl");
    double rhor = konstanten.search_con("rhor");
    double vr = konstanten.search_con("vr");
    double vrr = konstanten.search_con("vrr");
    double vyr = konstanten.search_con("vyr");
    double vyrr = konstanten.search_con("vyrr");
    double rhoul = konstanten.search_con("rhoul");
    double vul = konstanten.search_con("vul");
    double vrul = konstanten.search_con("vrul");
    double vyul = konstanten.search_con("vyul");
    double vyrul = konstanten.search_con("vyrul");
    double rhour = konstanten.search_con("rhour");
    double vur = konstanten.search_con("vur");
    double vrur = konstanten.search_con("vrur");
    double vyur = konstanten.search_con("vyur");
    double vyrur = konstanten.search_con("vyrur");

    switch(dimension)
    {
    case(1):

        cout << "Raster Initiieren - 0 Blockweise - 1 Schockwelle !" << endl;
        cout << "Wahl:";
        cin >> choice;

        switch(choice)
        {
        case(0):
            {
            //1D Riemann-Problem
            for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
            {
                xpos = mol + (mor-mol)*((double)(n-ordnung)/CELLSX);
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
        case(1):
            {
            //1D Schockwelle
            for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
            {
                xpos = mol + (mor-mol)*((double)(n-ordnung)/CELLSX);
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
    //2D
    case(2):
        ifstream input_winkel(winkel_in.c_str());


        ifstream input_kreis(kreis_in.c_str());
        string winkel_line;
        string kreis_line;

        symbol_table<double> symbol_table;
        expression<double> expression;
        parser<double> parser;

        cout << "Raster Initiieren - 0 Blockweise - 1 Winkel - 2 Blasen Simulation - 3 Speicherstand laden!" << endl;
        cout << "Wahl:";
        cin >> choice;

        switch(choice)
        {
        case(0):
            {
                //2D Riemann-Problem
                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        xpos = mol + (mor - mol)*(((double)n-ordnung)/CELLSX);
                        ypos = mul + (mol - mul)*(((double)m-ordnung)/CELLSY);
                        if(xpos <= 0.0)
                        {
                            if(ypos <= 0.0)
                            {
                                this->set_Zelle_d(rhoul,n,m);
                                this->set_Zelle_ux(vul,n,m);
                                this->set_Zelle_uxr(vrul,n,m);
                                this->set_Zelle_uy(vyul,n,m);
                                this->set_Zelle_uyr(vyrul,n,m);
                            }
                            else
                            {
                                this->set_Zelle_d(rhol,n,m);
                                this->set_Zelle_ux(vl,n,m);
                                this->set_Zelle_uxr(vrl,n,m);
                                this->set_Zelle_uy(vyl,n,m);
                                this->set_Zelle_uyr(vyrl,n,m);
                            }

                        }
                        else
                        {
                            if(ypos <= 0.0)
                            {
                                this->set_Zelle_d(rhour,n,m);
                                this->set_Zelle_ux(vur,n,m);
                                this->set_Zelle_uxr(vrur,n,m);
                                this->set_Zelle_uy(vyur,n,m);
                                this->set_Zelle_uyr(vyrur,n,m);
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
                }
                break;
            }
        case(1):
            {
                //2D Riemann-Problem mit Winkel
                getline(input_winkel,winkel_line);

                symbol_table.add_variable("ypos",ypos);
                symbol_table.add_variable("xpos",xpos);
                symbol_table.add_variable("dx",dx);
                symbol_table.add_variable("dy",dy);

                expression.register_symbol_table(symbol_table);

                parser.compile(winkel_line,expression);

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        xpos = mol + (mor - mol)*(((double)n-ordnung)/CELLSX);
                        ypos = mul + (mol - mul)*(((double)m-ordnung)/CELLSY);
                        if(ypos >= expression.value())
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
        case(2):
            {
                //Schockwelle-auf-Blase-Simulation
                getline(input_kreis,kreis_line);

                double radius = konstanten.search_con("radius");


                symbol_table.add_variable("ypos",ypos);
                symbol_table.add_variable("xpos",xpos);
                symbol_table.add_variable("dx",dx);
                symbol_table.add_variable("dy",dy);

                expression.register_symbol_table(symbol_table);

                parser.compile(kreis_line,expression);

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        xpos = mol + (mor - mol)*(((double)n-ordnung)/CELLSX);
                        ypos = mul + (mol - mul)*(((double)m-ordnung)/CELLSY);
                        if(radius >= expression.value())
                        {
                            set_Zelle_d(1,n,m);
                            set_Zelle_ux(0,n,m);
                            set_Zelle_uxr(0,n,m);
                            set_Zelle_uy(0,n,m);
                            set_Zelle_uyr(0,n,m);
                        }
                        else
                        {
                            set_Zelle_d(1000,n,m);
                            set_Zelle_ux(0,n,m);
                            set_Zelle_uxr(0,n,m);
                            set_Zelle_uy(0,n,m);
                            set_Zelle_uyr(0,n,m);
                        }
                    }
                }
                for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                {
                    set_Zelle_d(1323,CELLSX/2,m);
                    set_Zelle_ux(100000,CELLSX/2,m);
                    set_Zelle_uxr(0,CELLSX/2,m);
                    set_Zelle_uy(0,CELLSX/2,m);
                    set_Zelle_uyr(0,CELLSX/2,m);
                }
                break;
            }
        case(3):
            {
                //Speicherstand laden
                ifstream save_input(save_in.c_str());
                ifstream save_input_line;
                string line_open;
                string line;

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_d(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_ux(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_uy(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_uxr(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_uyr(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                getline(save_input,line_open);
                save_input_line.open(line_open.c_str());

                for(int n = ordnung ; n < CELLSX+ordnung+1 ; n++)
                {
                    for(int m = ordnung ; m < CELLSY+ordnung+1 ; m++)
                    {
                        getline(save_input_line,line);
                        line = line.substr(line.find("\t")+1);
                        set_Zelle_p(atof(line.c_str()),n,m);
                    }
                }
                save_input_line.close();

                break;
                //*/
            }
        }
        break;
    }
}

Raster::~Raster()
{
    delete[] zelle;
}



void Raster::bcondi(Konstanten konstanten, int* CELLS , int ordnung)
{


    double upbc = konstanten.search_con("upbc");
    double downbc = konstanten.search_con("downbc");
    double rightbc = konstanten.search_con("rightbc");
    double leftbc = konstanten.search_con("leftbc");
    //double xpos = 0.0 , ypos = 0.0;

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
            set_Zelle_ux( 0.0-get_Zelle(ordnung).ux ,0);
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
                            set_Zelle_ux( get_Zelle(x,CELLS[1]+ordnung).ux ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uxr( get_Zelle(x,CELLS[1]+ordnung).uxr ,x,CELLS[1]+ordnung+1);
                        }
                        else
                        {
                            set_Zelle_d( get_Zelle(x,CELLS[1]+ordnung).d,x,CELLS[1]+ordnung+1);
                            set_Zelle_uy( 0.0-get_Zelle(x,CELLS[1]+ordnung).uy ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uyr( get_Zelle(x,CELLS[1]+ordnung).uyr ,x,CELLS[1]+ordnung+1);
                            set_Zelle_ux( 0.0-get_Zelle(x,CELLS[1]+ordnung).ux ,x,CELLS[1]+ordnung+1);
                            set_Zelle_uxr( get_Zelle(x,CELLS[1]+ordnung).uxr ,x,CELLS[1]+ordnung+1);

                        }
                    }

                }
            }
        }

        break;
    }
}

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
