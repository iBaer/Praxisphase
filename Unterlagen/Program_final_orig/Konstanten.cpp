#include "Konstanten.h"
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace exprtk;

Konstanten::Konstanten(string input_const)
{
        ifstream input(input_const.c_str());
        string line, name, calc;
        double value;
        if(!input_const.empty())
        {
            if(input.is_open())
            {
                //einlesen Leerzeichen als Trenzeichen
                while(getline(input,line))
                {
                    name = line.substr(0,line.find(" "));
                    value = atof(line.substr(line.find(" ")+1, line.size()-line.find(" ")-1).c_str());
                    const_name.push_back(name);
                    const_value.push_back(value);
                }
                input.close();
            }else
            {
                cout <<"Fehler beim öffnen der Konstanten!";
            }
        }
}

symbol_table<double> Konstanten::register_constants(symbol_table<double> in)
{
    unsigned int n;
    for( n = 0 ; n < const_name.size() ; n++)in.add_variable(const_name.at(n),const_value.at(n));
    return in;
}

double Konstanten::search_con(string in)
{
    unsigned int n=0;
    bool found = false;

    while(n < const_name.size() && !found)
    {
        if((const_name.at(n).compare(in)) == 0)found = true;
        else n++;
    }
    return const_value.at(n);
}

