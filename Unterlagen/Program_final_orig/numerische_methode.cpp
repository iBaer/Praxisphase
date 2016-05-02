#include "numerische_methode.h"
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <lapackpp/gmd.h>
#include <lapackpp/lavd.h>
#include <lapackpp/laslv.h>



using namespace std;
using namespace exprtk;


numerische_methode::numerische_methode(string method, string const_in, string formel_in, string winkel_in, string kreis_in,  string save_in) : konstanten(const_in) , raster(konstanten, winkel_in, kreis_in, save_in), gs(formel_in, konstanten)
{
    name = method;
    ordnung = 1;
    dimension = konstanten.search_con("dimension");
    CELLS = new int(dimension);

    double cref = konstanten.search_con("cref");
    double done = konstanten.search_con("done");
    double ccl = konstanten.search_con("ccl");

    mor = konstanten.search_con("mor");
    mol = konstanten.search_con("mol");
    mur = konstanten.search_con("mur"); //2D
    mul = konstanten.search_con("mul"); //2D
    timeou = konstanten.search_con("timeou");
    steps = 0;
    maxnt = konstanten.search_con("maxnt");
    CELLS[0] = konstanten.search_con("CELLSX");
    CELLS[1] = konstanten.search_con("CELLSY");
    g = konstanten.search_con("g");
    dx = (mor-mol)/(double)CELLS[0];
    dy = (mor-mur)/(double)CELLS[1];
    dt = 0.0;

    double rhol = konstanten.search_con("rhol");

    double alfll = 1.0 - rhol/done + ccl *(rhol/done);
    double dll = ccl * (rhol/alfll);
    double pll = cref*pow(dll,g);

    ct = pll/pow(rhol,g);


}

void numerische_methode::start_method()
{
    double time = 0.0, timetol = 0.000001, timedif = 1.0;
    int step_output = 0;

    cout << "Stepwise Output[1: true :: 0: false]: ";
    cin >> step_output;

    if(step_output==1)write();
    

    for( int n = 1 ; n <= maxnt && timedif > timetol ; n++)
    {
        cout << n << " : " << maxnt << endl;
        raster.bcondi(konstanten,CELLS,ordnung);
        if(step_output==1)write();
        cout << "Neuen Zeitschritt berechnen..." << endl;
        time = cflcon(n,time);
        update(calc_method_flux());
        timedif = abs(time-timeou);
        steps = n;
    }
    write();
}


double numerische_methode::cflcon(int n, double time)
{
    double cref = konstanten.search_con("cref");
    double cfl = konstanten.search_con("cfl");
    double ccl = konstanten.search_con("ccl");
    double done = konstanten.search_con("done");
    int variante = (int) konstanten.search_con("variante");
    double gi = 1.0/g;

    double maxd = 0.0 , maxu = 0.0 , maxur = 0.0, maxp = 0.0, maxuy = 0.0, maxuyr = 0.0;

    double smax=0.0,maxs=0.0;
    ifstream cfl_in;
    ifstream matrix_in;
    ifstream matrix_in2;

    if(dimension == 1)cfl_in.open("cfl.in");
    else cfl_in.open("cfl2d.in");

    switch(dimension)
    {
    case(1):
        {
            switch(variante)
            {
            case(1):
                matrix_in.open("matrix_variante1.in");
                break;
            case(2):
                matrix_in.open("matrix_variante2.in");
                break;
            case(3):
                matrix_in.open("matrix_variante3.in");
                break;
            }
            break;
        }
    case(2):
        {
            matrix_in.open("matrix2d_1.in");
            matrix_in2.open("matrix2d_2.in");
        }
    }

    string line;

    string am_calc, xs_calc;
    double p=0.0,ux=0.0,d=0.0,uxr=0.0,dtwo=0.0, uy = 0.0, uyr = 0.0;
    getline(cfl_in,am_calc);
    getline(cfl_in,xs_calc);
    cfl_in.close();
    symbol_table<double> st;
    expression<double> expression_am, expression_xs;
    parser<double> parser;

    switch(dimension)
    {
    case(1):
        {
            //1D
            if((int)konstanten.search_con("calceigv") == 1) //1 = true, smax wird über Eigenwerte bestimmt
            {
                string matrix_array[9];
                for(int n = 0 ; n < 9 ; n++)
                {
                    getline(matrix_in,matrix_array[n]);
                }

                double values[9];

                //Schritt 1: Maxima finden
                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                    if(raster.get_Zelle(x).d > maxd)maxd = raster.get_Zelle(x).d;
                    if(raster.get_Zelle(x).ux > maxu)maxu = raster.get_Zelle(x).ux;
                    if(raster.get_Zelle(x).uxr > maxur)maxur = raster.get_Zelle(x).uxr;
                }

                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                    switch(variante)
                    {
                    case(1):
                        p = ct* pow(raster.get_Zelle(x).d,g);
                        raster.set_Zelle_p(p,x);
                        dtwo = pow((p/cref),gi);
                        break;
                    case(2):
                        dtwo = ccl/((1/raster.get_Zelle(x).d)-((1-ccl)/done));
                        p = cref * pow(dtwo,g);
                        raster.set_Zelle_p(p,x);
                        break;
                    case(3):
                        dtwo = ccl/((1/raster.get_Zelle(x).d)-((1-ccl)/done));
                        p = ct* pow(raster.get_Zelle(x).d,g);
                        raster.set_Zelle_p(p,x);
                        break;
                    }
                    if(maxp < p)maxp = p;
                }

                maxu = maxu * maxd;


                //Schritt 2: Einsetzen

                symbol_table<double> st_values;
                st_values.add_variable("ct",ct);
                st_values.add_variable("dtwo",dtwo);
                st_values.add_variable("uone", maxd);
                st_values.add_variable("utwo", maxu);
                st_values.add_variable("uthree", maxur);
                st_values = konstanten.register_constants(st_values);

                expression<double> expressions_values[9];

                for(int n = 0 ; n < 9 ; n++)
                {
                    expressions_values[n].register_symbol_table(st_values);
                    parser.compile(matrix_array[n],expressions_values[n]);
                    values[n] = expressions_values[n].value();
                    //cout << values[n] << endl;
                }

                //Schritt 3: Berechnen der Eigenwerte


                LaGenMatDouble A(values,3,3,true);
                LaVectorDouble real(3);
                LaVectorDouble img(3);
                LaGenMatDouble vr(3,3);

                LaEigSolve(A,real,img,vr);

                //Schritt 4: Höchsten Eigenwert suchen
                for(int n = 0 ; n < 3 ; n++)
                {
                    if(smax < real(n))smax = real(n);
                }


                dt = cfl*dx/smax;

                if(n <= konstanten.search_con("teilerend"))dt = dt * konstanten.search_con("teiler");

                if((time+dt)> timeou)dt=timeou-time;

                cout << "Neues delta t ist: \t" << dt << endl;

                time = time + dt;

            }
            else
            {
                //dt über Näherung berechnen
                st.add_variable("p",p);
                st.add_variable("ux",ux);
                st.add_variable("d",d);
                st.add_variable("uxr",uxr);
                st.add_variable("dtwo",dtwo);
                st = konstanten.register_constants(st);
                expression_am.register_symbol_table(st);
                expression_xs.register_symbol_table(st);
                parser.compile(am_calc,expression_am);
                parser.compile(xs_calc,expression_xs);

                for(int i = 0 ; i < CELLS[0]+2*ordnung+1 ; i++)
                {
                    d = raster.get_Zelle(i).d;
                    uxr = raster.get_Zelle(i).uxr;
                    ux = raster.get_Zelle(i).ux;


                    switch(variante)
                    {
                    case(1):
                        {
                        p = ct* pow(d,g);
                        raster.set_Zelle_p(p,i);
                        dtwo = pow((p/cref),gi);
                        break;
                        }
                    case(2):
                        {
                        dtwo = ccl/((1/d)-((1-ccl)/done));
                        p = cref * pow(dtwo,g);
                        raster.set_Zelle_p(p,i);
                        break;
                        }
                    case(3):
                        {
                        dtwo = ccl/((1/d)-((1-ccl)/done));
                        p = ct* pow(d,g);
                        raster.set_Zelle_p(p,i);
                        break;
                        }
                    }

                    maxs = abs(ux)+expression_am.value()+expression_xs.value();
                    if(maxs > smax)smax = maxs;
                }



                dt = cfl*dx/smax;

                if(n <= konstanten.search_con("teilerend"))dt = dt * konstanten.search_con("teiler");
                //dt = 0,00002;

                if((time+dt)> timeou)dt=timeou-time;

                cout << "Neues delta t ist: \t" << dt << endl;

                time = time + dt;
            }

            break;


        }
    case(2):
        {
         //2D     

            if((int)konstanten.search_con("calceigv") != 1) //1 = true, smax wird über Eigenwerte bestimmt
            {
                //über Näherung
                st.add_variable("p",p);
                st.add_variable("uy",ux);
                st.add_variable("uy",uy);
                st.add_variable("d",d);
                st.add_variable("uxr",uxr);
                st.add_variable("uyr",uyr);
                st.add_variable("dtwo",dtwo);
                st = konstanten.register_constants(st);
                expression_am.register_symbol_table(st);
                expression_xs.register_symbol_table(st);
                parser.compile(am_calc,expression_am);
                parser.compile(xs_calc,expression_xs);

                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                    for(int y = 0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                    {
                        d = raster.get_Zelle(x,y).d;
                        uxr = raster.get_Zelle(x,y).uxr;
                        ux = raster.get_Zelle(x,y).ux;
                        uy = raster.get_Zelle(x,y).uy;
                        uyr = raster.get_Zelle(x,y).uyr;

                        p = ct* pow(d,g);
                        //cout << "p:" << p << endl;
                        raster.set_Zelle_p(p,x,y);

                        dtwo = pow((p/cref),gi);

                        maxs = abs(ux)+expression_am.value()+expression_xs.value();
                        if(maxs > smax)smax = maxs;
                    }
                }

                dt = cfl*max(dx,dy)/smax; //maximum von dy und dx

                //cout << "dt: " << dt << endl;
                //cout << "dx: " << dx << endl;
                //cout << "dy: " << dy << endl;



                if(n <= konstanten.search_con("teilerend"))dt = dt * konstanten.search_con("teiler");

                if((time+dt)> timeou)dt=timeou-time;

                cout << "Neues delta t ist: \t" << dt << endl;

                time = time + dt;
            }
            else
            {
                //über Eigenwerte
                string matrix_array1[25];
                string matrix_array2[25];
                double values1[25];
                double values2[25];

                for(int n = 0 ; n < 9 ; n++)
                {
                    getline(matrix_in,matrix_array1[n]);
                }

                for(int n = 0 ; n < 9 ; n++)
                {
                    getline(matrix_in2,matrix_array2[n]);
                }

                //Schritt 1: Maxima finden
                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                    for(int y = 0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                    {
                        if(raster.get_Zelle(x,y).d > maxd)maxd = raster.get_Zelle(x).d;
                        if(raster.get_Zelle(x,y).ux > maxu)maxu = raster.get_Zelle(x).ux;
                        if(raster.get_Zelle(x,y).uy > maxu)maxuy = raster.get_Zelle(x).uy;
                        if(raster.get_Zelle(x,y).uxr > maxur)maxur = raster.get_Zelle(x).uxr;
                        if(raster.get_Zelle(x,y).uyr > maxuyr)maxuyr = raster.get_Zelle(x).uyr;
                    }
                }

                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                    for(int y = 0 ; x < CELLS[1]+2*ordnung+1 ; x++)
                    {
                        p = ct* pow(raster.get_Zelle(x,y).d,g);
                        raster.set_Zelle_p(p,x,y);
                        dtwo = pow((p/cref),gi);

                        if(maxp < p)maxp = p;
                    }
                }

                maxu = maxu * maxd;
                maxuy = maxuy * maxd;


                //Schritt 2: Einsetzen

                symbol_table<double> st_values;
                st_values.add_variable("p",maxp);
                st_values.add_variable("dtwo",dtwo);
                st_values.add_variable("uone", maxd);
                st_values.add_variable("utwo", maxu);
                st_values.add_variable("uthree", maxuy);
                st_values.add_variable("ufour", maxur);
                st_values.add_variable("ufive", maxuyr);
                st_values = konstanten.register_constants(st_values);

                expression<double> expressions_values[25];
                expression<double> expressions_values_b[25];

                for(int n = 0 ; n < 25 ; n++)
                {
                    expressions_values[n].register_symbol_table(st_values);
                    parser.compile(matrix_array1[n],expressions_values[n]);
                    values1[n] = expressions_values[n].value();

                    expressions_values_b[n].register_symbol_table(st_values);
                    parser.compile(matrix_array2[n],expressions_values_b[n]);
                    values2[n] = expressions_values_b[n].value();
                }

                //for(int n = 0 ; n < 25 ; n++)cout << values1[n] << endl;

                //cout << endl;

                //for(int n = 0 ; n < 25 ; n++)cout << values2[n] << endl;


                //Schritt 3: Berechnen der Eigenwerte

                LaGenMatDouble A(values1,5,5,true);
                LaGenMatDouble B(values2,5,5,true);
                LaVectorDouble real(5);
                LaVectorDouble real_b(5);
                LaVectorDouble img(5);
                LaVectorDouble img_b(5);
                LaGenMatDouble vr(5,5);
                LaGenMatDouble vr_b(5,5);

                LaEigSolve(A,real,img,vr);
                LaEigSolve(B,real_b,img_b,vr_b);


                //Schritt 4: Höchsten Eigenwert suchen
                for(int n = 0 ; n < 5 ; n++)
                {
                    if(smax < real(n))smax = real(n);
                    if(smax < real_b(n))smax = real_b(n);
                }

                dt = cfl*dx/smax;

                if(n <= konstanten.search_con("teilerend"))dt = dt * konstanten.search_con("teiler");
                //dt = 0,00002;

                if((time+dt)> timeou)dt=timeou-time;

                cout << "Neues delta t ist: \t" << dt << endl;

                time = time + dt;

            }

            break;
        }
    }
    return time;
}

void numerische_methode::update(vector< vector <vector< vector<double> > > > fi)
{
    cout << "Zellen updaten..." << endl;

    //******************Dieser Abschnitt für Anpassbare rückrechnung gedacht. hatte nicht funktioniert*********
    //******************Wegen Zeitmangel nicht mehr drauf zurück gekommen************************************
    /*
    Zelle tmp_cell;
    double ux = 0.0;
    double uxr = 0.0;
    double uy = 0.0;
    double uyr = 0.0;
    double d = 0.0;
    double p = 0.0;
    double cs = 0.0;

    

    symbol_table<double> symbol_table;
    symbol_table.add_variable("ux",ux);
    symbol_table.add_variable("uxr",uxr);
    symbol_table.add_variable("uy",uy);
    symbol_table.add_variable("uyr",uyr);
    symbol_table.add_variable("d",d);
    symbol_table.add_variable("p",p);
    symbol_table.add_variable("cs",cs);

    symbol_table = konstanten.register_constants(symbol_table);

    expression<double> expression0;
    expression<double> expression1;
    expression<double> expression2;
    expression<double> expression3;
    expression<double> expression4;

    expression0.register_symbol_table(symbol_table);
    expression1.register_symbol_table(symbol_table);
    expression2.register_symbol_table(symbol_table);
    expression3.register_symbol_table(symbol_table);
    expression4.register_symbol_table(symbol_table);

    parser<double> parser0;
    parser<double> parser1;
    parser<double> parser2;
    parser<double> parser3;
    parser<double> parser4;


    parser0.compile( gs.uback.at(0), expression0);
    parser1.compile( gs.uback.at(1), expression1);
    parser2.compile( gs.uback.at(2), expression2);
    if(dimension == 2)
    {
        parser3.compile( gs.uback.at(3), expression3);
        parser4.compile( gs.uback.at(4), expression4);
    }
    */
    //*********************************************************************************

    double dtodx = dt/dx;
    double dtody = dt/dy;

    double csm[fi.size()][raster.getwidth()][raster.getheight()];

    //cout << fi.size() << " - " << raster.getwidth() << " - " << raster.getheight() << endl;

    if(dimension == 1)
    {
        //U aussrechnen
        for(int i = 0 ; i < CELLS[0]+2*ordnung+1 ; i++)
        {
            for(unsigned int k = 0 ; k < gs.u.size() ; k++)
            {
                csm[k][i][0] = gs.solve_u(raster.get_Zelle(i),k);
            }
        }
        //Update schritt
        for(int i = ordnung ; i < CELLS[0]+ordnung+1 ; i++)
        {
            for(unsigned int k = 0 ; k < fi.size() ; k++)
            {
                csm[k][i][0] = csm[k][i][0] + dtodx*(fi.at(k).at(i-1).at(0).at(0)-fi.at(k).at(i).at(0).at(0));
            }
        }
        //setzen der Zellen
        for(int i = ordnung ; i < CELLS[0]+ordnung+1 ; i++ )
        {
            /*
            tmp_cell = raster.get_Zelle(i);
            d = tmp_cell.d;

            p = tmp_cell.p;
            ux = tmp_cell.ux;
            uxr = tmp_cell.ux_r;
            uy = tmp_cell.uy;
            uyr = tmp_cell.uy_r;


            cs = csm[0][i][0];
            raster.set_Zelle_d(expression0.value(),i);

            cs = csm[1][i][0];
            raster.set_Zelle_ux(expression1.value(),i);

            cs = csm[2][i][0];
            raster.set_Zelle_ux_r(expression2.value(),i);

            raster.set_Zelle_uy(0,i);

            raster.set_Zelle_uy_r(0,i);
            */


            raster.set_Zelle_d(csm[0][i][0],i);
            raster.set_Zelle_ux(csm[1][i][0]/raster.get_Zelle(i).d,i);
            raster.set_Zelle_uxr(csm[2][i][0],i);
        }
    }
    else
    {
        //U ausrechnen
        for(int x = ordnung ; x < CELLS[0]+ordnung+1 ; x++)
        {
            for(int y = ordnung ; y < CELLS[1]+ordnung+1 ; y++)
            {
                for(unsigned int k = 0 ; k < fi.size() ; k++)
                {
                    //cout << "csm: " << csm[k][x][y] << endl;
                    csm[k][x][y] = gs.solve_u(raster.get_Zelle(x,y),k);
                }
            }

        }
        //Update Schritt
        for(int x = ordnung ; x < CELLS[0]+ordnung+1 ; x++)
        {
            for(int y = ordnung ; y < CELLS[1]+ordnung+1 ; y++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    csm[k][x][y] = csm[k][x][y] + dtodx*(fi.at(k).at(x-1).at(y).at(0)-fi.at(k).at(x).at(y).at(0)) + dtody*(fi.at(k).at(x).at(y-1).at(1)-fi.at(k).at(x).at(y).at(1));
                    //if(k==3)cout << csm[k][x][y] << endl;
                    //csm[k][x][y] = 0.25*(csm[k][x+1][y]+csm[k][x-1][y]+csm[k][x][y+1]+csm[k][x][y-1])-0.5*dt/dx*(f[k][x+1][y]-f[k][x-1][y])-0.5*dt/dy*(g[k][x][y+1]-g[k][x][y-1]);
                    //if(k==3)cout << csm[k][x][y] << endl;
                }
            }
        }
        //setzen der Zellen
        for(int x = ordnung ; x < CELLS[0]+ordnung+1 ; x++ )
        {
            for(int y = ordnung ; y < CELLS[1]+ordnung+1 ; y++)
            {
                raster.set_Zelle_d(csm[0][x][y],x,y);
                raster.set_Zelle_ux(csm[1][x][y]/raster.get_Zelle(x,y).d,x,y);
                raster.set_Zelle_uy(csm[2][x][y]/raster.get_Zelle(x,y).d,x,y);
                raster.set_Zelle_uxr(csm[3][x][y],x,y);
                raster.set_Zelle_uyr(csm[4][x][y],x,y);
            }
        }
    }
}


void numerische_methode::write()
{
    double xpos = 0.0;
    double ypos = 0.0;
    double p = 0.0;

    string added = to_string(CELLS[0])+"x"+to_string(CELLS[1])+"_"+name+"_"+to_string(dimension)+"d_"+to_string((int)konstanten.search_con("variante"))+".variant_"+"divider"+to_string(konstanten.search_con("teiler"))+"till"+to_string((int)konstanten.search_con("teilerend"))+"_"+to_string(steps)+"Steps";

    switch(raster.getdim())
    {
    case(1):
        {
            string d_path = "d_"+added;
            string uxr_path = "uxr_"+added;
            string ux_path = "ux_"+added;
            string p_path = "p_"+added;

            ofstream d_out(d_path.c_str());
            ofstream uxr_out(uxr_path.c_str());
            ofstream ux_out(ux_path.c_str());
            ofstream p_out(p_path.c_str());

            for(int i = ordnung ; i < CELLS[0]+ordnung+1 ; i++)
            {
                xpos = (mol + (mor-mol)*((double)(i-ordnung)/CELLS[0]));

                p = ct*pow(raster.get_Zelle(i).d,g);

                d_out << fixed << setprecision(1) << xpos << " \t"  << setprecision(10) << raster.get_Zelle(i).d << "\n";
                uxr_out << fixed << setprecision(1) << xpos << " \t" << setprecision(10) << raster.get_Zelle(i).uxr << "\n";
                ux_out << fixed << setprecision(1) << xpos << " \t" << setprecision(10) << raster.get_Zelle(i).ux << "\n";
                p_out << fixed << setprecision(1) << xpos << " \t" << setprecision(10) << p << "\n";

            }

            d_out.close();
            uxr_out.close();
            ux_out.close();
            p_out.close();
            break;
        }
    case(2):
        {
            //2D

            string d_path = "d_"+added;
            string uxr_path = "uxr_"+added;
            string uyr_path = "uyr_"+added;
            string ux_path = "ux_"+added;
            string uy_path = "uy_"+added;
            string p_path = "p_"+added;

            ofstream d_out(d_path.c_str());
            ofstream uxr_out(uxr_path.c_str());
            ofstream ux_out(ux_path.c_str());
            ofstream uyr_out(uyr_path.c_str());
            ofstream uy_out(uy_path.c_str());
            ofstream p_out(p_path.c_str());

            for(int x = ordnung ; x < CELLS[0]+ordnung+1 ; x++)
            {
                for(int y = ordnung ; y < CELLS[1]+ordnung+1 ; y++)
                {
                    xpos = mol + (mor-mol)*((double)(x-ordnung)/CELLS[0]);
                    ypos = mur + (mor - mur)*((double)(y-ordnung)/CELLS[1]);

                    p = ct*pow(raster.get_Zelle(x,y).d,g);

                    d_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << raster.get_Zelle(x,y).d << "\n";
                    uxr_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << raster.get_Zelle(x,y).uxr << "\n";
                    ux_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << raster.get_Zelle(x,y).ux << "\n";
                    uyr_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << raster.get_Zelle(x,y).uyr << "\n";
                    uy_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << raster.get_Zelle(x,y).uy << "\n";
                    p_out << fixed << setprecision(1) << xpos << " " << fixed << setprecision(1) << ypos << " \t" << setprecision(10) << p << "\n";
                }
            }

            d_out.close();
            uxr_out.close();
            ux_out.close();
            uyr_out.close();
            uy_out.close();
            p_out.close();
            break;
        }
    }
}

