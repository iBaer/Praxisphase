#include <sys/time.h>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "exprtk.hpp"

using namespace std;

double start;
struct timeval mytime;
double duration;

template <typename T>
void trig_function()
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;

   std::string expression_string = "clamp(-1.0,sin(2 * pi * x) + cos(x / 2 * pi),+1.0)";

   T x;

   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_constants();

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);

   for (x = T(-5); x <= T(+5); x += T(0.001))
   {
      T y = expression.value();
      printf("%19.15f\t%19.15f\n",x,y);
   }
}

int main()
{
// ***
gettimeofday(&mytime,NULL);
start = (double) (mytime.tv_sec) + ((double) (mytime.tv_usec))/1.0e6;
// ***
   trig_function<double>();
// ***
gettimeofday(&mytime,NULL);
duration = (double) (mytime.tv_sec) + ((double) (mytime.tv_usec))/1.0e6 - start;
// ***
cout << "Time for step: "<<duration<<endl; 
   return 0;
}
