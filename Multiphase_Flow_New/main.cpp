#include "computation.h"
#include "LaxFriedrichMethod.h"
#include "constants.h"
#include "force.h"


using namespace std;

int main()
{
	Constants constants = Constants::instance();
	Computation computation = Computation::instance(&constants);
    Grid grid = Grid(&constants, "save.in");
    Solver* solver;
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
        solver = new LaxFriedrichMethod(&constants,&computation, &grid);
    }
    else if (method == 2)
    {
    	solver = new Force(&constants,&computation, &grid);
    }
    else{
    	cout << "No valid method selected"<<endl<<"Closing Simulator"<<endl;
    	return 0;
    }

    numerische_methode num_meth(solver,&constants,&computation, &grid);
    num_meth.start_method();
    return 0;

}
