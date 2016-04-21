#include "LaxFriedrichMethod.h"
#include "FORCE.h"

int main()
{

    int d = 1;

    std::cout << "Wahl! 1 Lax-Friedrich - 2 FORCE :";
    std::cin >> d;

    if(d==1)
    {
        LaxFriedrichMethod lax("gas-liquid dia.in","formeln","winkel.in","kreis.in","save.in");
        lax.start_method();
    }
    else
    {
        FORCE force("gas-liquid dia.in","formeln","winkel.in","kreis.in","save.in");
        force.start_method();
    }


    return 0;
}
