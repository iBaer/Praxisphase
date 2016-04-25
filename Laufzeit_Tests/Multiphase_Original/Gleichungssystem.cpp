#include "Gleichungssystem.h"
#include <fstream>
using namespace std;
using namespace exprtk;

Gleichungssystem::Gleichungssystem(string path, Konstanten c) : cons("")
{
    cons=c;
    int p=0;
    int dim = cons.search_con("dimension");
    int variante = cons.search_con("variante");
    //cout << dim << endl;
    unsigned int n=0;
    string buffer;

    if(dim == 1)
    {
        if(variante == 3)
        {
            path = path + "_variante3.in";
        }
        else
        {
            path = path + ".in";
        }
    }
    else
    {
        if(variante == 3)
        {
            path = path + "2d_variante3.in";
        }
        else
        {
            path = path + "2d.in";
        }
    }

    cout << "Verwendet wird: " << path << endl;

    ifstream input(path.c_str());
    string line;

     if(input.is_open())
    {
        getline(input,line);

        p = line.find(" ",n);
        buffer = line.substr(n, p-n);
        uback.push_back(buffer);
        n = p+1;


        p = line.find(" ",n);
        buffer = line.substr(n, p-n);
        uback.push_back(buffer);
        n = p+1;

        p = line.find(" ",n);
        buffer = line.substr(n, p-n);
        uback.push_back(buffer);

        if(dim == 2)
        {
            n = p+1;
            p = line.find(" ",n);
            buffer = line.substr(n, p-n);
            uback.push_back(buffer);
            n = p+1;

            p = line.find(" ",n);
            buffer = line.substr(n, p-n);
            uback.push_back(buffer);
            n = p+1;
        }

        while(getline(input,line))
        {
            n=0;
            p=0;

            p = line.find(" ",n);
            buffer = line.substr(n, p-n);
            u.push_back(buffer);
            n = p+1;


            p = line.find(" ",n);
            buffer = line.substr(n, p-n);
            f_u.push_back(buffer);


            if(dim == 2)
            {
                n = p+1;
                p = line.find(" ",n);
                buffer = line.substr(n, p-n);
                g_u.push_back(buffer);
            }
        }
    }
    else
    {
        cout << "Fehler beim öffnen!" << endl;
    }
    input.close();


    symbol_table_u.add_variable("ux",ux);
    symbol_table_u.add_variable("uxr",uxr);
    symbol_table_u.add_variable("uy",uy);
    symbol_table_u.add_variable("uyr",uyr);
    symbol_table_u.add_variable("d",d);
    symbol_table_u.add_variable("p",this->p);

    symbol_table_f_u.add_variable("ux",ux);
    symbol_table_f_u.add_variable("uxr",uxr);
    symbol_table_f_u.add_variable("uy",uy);
    symbol_table_f_u.add_variable("uyr",uyr);
    symbol_table_f_u.add_variable("d",d);
    symbol_table_f_u.add_variable("p",this->p);

    symbol_table_g_u.add_variable("ux",ux);
    symbol_table_g_u.add_variable("uxr",uxr);
    symbol_table_g_u.add_variable("uy",uy);
    symbol_table_g_u.add_variable("uyr",uyr);
    symbol_table_g_u.add_variable("d",d);
    symbol_table_g_u.add_variable("p",this->p);

    symbol_table_u = cons.register_constants(symbol_table_u);
    symbol_table_f_u = cons.register_constants(symbol_table_f_u);
    symbol_table_g_u = cons.register_constants(symbol_table_g_u);

    for(int i=0;i<=4;i++){
    	expression_u[i].register_symbol_table(symbol_table_u);
	expression_f_u[i].register_symbol_table(symbol_table_f_u);
	expression_g_u[i].register_symbol_table(symbol_table_g_u);
	parser.compile( u.at(i), expression_u[i]);
	parser.compile( f_u.at(i), expression_f_u[i]);
	parser.compile( g_u.at(i), expression_g_u[i]);

    }

}



double Gleichungssystem::solve_u(Zelle zelle, int pos)
{


    double result=0.0;


    ux = zelle.ux;
    uxr = zelle.uxr;
    uy = zelle.uy;
    uyr = zelle.uyr;
    d = zelle.d;
    p = zelle.p;

    result = expression_u[pos].value();



    //if(pos==1)cout << "u: " << result << "\r";



    return result;
}

double Gleichungssystem::solve_f_u(Zelle zelle, int pos)
{
    double result=0.0;


    ux = zelle.ux;
    uxr = zelle.uxr;
    uy = zelle.uy;
    uyr = zelle.uyr;
    d = zelle.d;
    p = zelle.p;

    result = expression_f_u[pos].value();

    //if(pos==1)cout << "f: " << result << "\r";

    return result;
}

double Gleichungssystem::solve_g_u(Zelle zelle, int pos)
{
    double result=0.0;


    ux = zelle.ux;
    uxr = zelle.uxr;
    uy = zelle.uy;
    uyr = zelle.uyr;
    d = zelle.d;
    p = zelle.p;

    result = expression_g_u[pos].value();

    //if(pos==1)cout << "g: " << result << "\r";

    return result;
}






