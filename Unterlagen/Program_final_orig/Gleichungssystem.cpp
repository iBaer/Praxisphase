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
}

double Gleichungssystem::solve_u(Zelle zelle, int pos)
{
    double result=0.0;
    symbol_table<double> symbol_table;

    symbol_table.add_variable("ux",zelle.ux);
    symbol_table.add_variable("uxr",zelle.uxr);
    symbol_table.add_variable("uy",zelle.uy);
    symbol_table.add_variable("uyr",zelle.uyr);
    symbol_table.add_variable("d",zelle.d);
    symbol_table.add_variable("p",zelle.p);

    symbol_table = cons.register_constants(symbol_table);

    expression<double> expression;
    expression.register_symbol_table(symbol_table);

    parser<double> parser;
    parser.compile( u.at(pos), expression);
    result = expression.value();

    //if(pos==1)cout << "u: " << result << "\r";

    return result;
}

double Gleichungssystem::solve_f_u(Zelle zelle, int pos)
{
    double result=0.0;
    symbol_table<double> symbol_table;

    symbol_table.add_variable("ux",zelle.ux);
    symbol_table.add_variable("uxr",zelle.uxr);
    symbol_table.add_variable("uy",zelle.uy);
    symbol_table.add_variable("uyr",zelle.uyr);
    symbol_table.add_variable("d",zelle.d);
    symbol_table.add_variable("p",zelle.p);

    symbol_table = cons.register_constants(symbol_table);

    expression<double> expression;
    expression.register_symbol_table(symbol_table);

    parser<double> parser;
    parser.compile( f_u.at(pos), expression);
    result = expression.value();

    //if(pos==1)cout << "f: " << result << "\r";

    return result;
}

double Gleichungssystem::solve_g_u(Zelle zelle, int pos)
{
    double result=0.0;
    symbol_table<double> symbol_table;

    symbol_table.add_variable("ux",zelle.ux);
    symbol_table.add_variable("uxr",zelle.uxr);
    symbol_table.add_variable("uy",zelle.uy);
    symbol_table.add_variable("uyr",zelle.uyr);
    symbol_table.add_variable("d",zelle.d);
    symbol_table.add_variable("p",zelle.p);

    symbol_table = cons.register_constants(symbol_table);

    expression<double> expression;
    expression.register_symbol_table(symbol_table);

    parser<double> parser;
    parser.compile( g_u.at(pos), expression);
    result = expression.value();

    //if(pos==1)cout << "g: " << result << "\r";

    return result;
}






