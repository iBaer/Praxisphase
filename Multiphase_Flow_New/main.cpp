#include "LaxFriedrichMethod.h"
#include "FORCE.h"
#include "constants.h"
#include "Gleichungssystem.h"


using namespace std;

int main()
{
	Constants constants = Constants::instance();
	Computation computation = Computation::instance(&constants);
    Grid grid = Grid(&constants, "save.in");

    int method = 1;
    cout << "####################################"<<endl;
    cout << "Multiphase Flow Simulator"<<endl;
    cout << "####################################"<<endl<<endl;

    cout << "Please select the numerical method for flux calculation:"<<endl;
    cout << "1 - Lax-Friedrich"<<endl;
    cout << "2 - FORCE"<<endl;
    cin >> method;

    /* Weitere Wahlmöglichkeiten hier!
     * Splitting/Unsplitting (2D+3D)
     * Schrittweiser Output
     * Raster Initialisierung
     */

    if (method == 1)
    {
        LaxFriedrichMethod lax(&constants,&computation, &grid);
        numerische_methode num_meth(&lax,&constants,&computation, &grid);
        num_meth.start_method();
    }
    if (method == 2)
    {
    	FORCE force(&constants,&computation, &grid);
		numerische_methode num_meth(&force,&constants,&computation, &grid);
		num_meth.start_method();
    }

    return 0;

}
