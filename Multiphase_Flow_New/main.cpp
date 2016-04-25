#include "LaxFriedrichMethod.h"
#include "FORCE.h"

int main()
{

    int method = 1;
    std::cout << "Multiphase Flow Simulator"<<std::endl;
    std::cout << "Wahl! 1 Lax-Friedrich - 2 FORCE :";
    std::cin >> method;

    if (method == 1)
    {
        LaxFriedrichMethod lax("gas-liquid.in","formeln", "save.in");
        lax.start_method();
    }
    if (method == 2)
    {
        FORCE force("gas-liquid.in","formeln", "save.in");
        force.start_method();
    }

    return 0;

}
