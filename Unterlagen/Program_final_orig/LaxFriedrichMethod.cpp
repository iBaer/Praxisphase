#include "LaxFriedrichMethod.h"

using namespace std;

LaxFriedrichMethod::LaxFriedrichMethod(std::string const_in, std::string formel_in, std::string winkel_in, std::string kreis_in, std::string save_in) : numerische_methode("Lax-Friedrich", const_in, formel_in, winkel_in,kreis_in,save_in)
{

}

vector < vector< vector< vector <double> > > > LaxFriedrichMethod::calc_method_flux()
{
    cout << "Lax-Friedrich Fluss berechnen..." << endl;
    double cs[gs.u.size()][raster.getwidth()][raster.getheight()];
    double f[gs.u.size()][raster.getwidth()][raster.getheight()];
    double g[gs.u.size()][raster.getwidth()][raster.getheight()];

    vector< vector< vector< vector< double > > > > fi (gs.u.size(), vector< vector< vector<double> > > ( raster.getwidth(), vector< vector<double> > ( raster.getheight() , vector<double> (dimension , 0.0) ) ) );
    switch(dimension)
    {
    case(1):
        {
            //1D
            //Berechne U und F
            for(int i = 0 ; i < CELLS[0]+2*ordnung+1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    cs[k][i][0] = gs.solve_u(raster.get_Zelle(i),k);

                }

                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    f[k][i][0] = gs.solve_f_u(raster.get_Zelle(i),k);

                }
            }
            //Berechne Lax-Friedrich-Fluss
            for(int i = 0 ; i < CELLS[0] +ordnung+1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    fi.at(k).at(i).at(0).at(0) = 0.5*(f[k][i][0] + f[k][i+1][0]) + 0.5*(dx/dt)*(cs[k][i][0] - cs[k][i+1][0]);
                    //cout << fi.at(k).at(i).at(0).at(0) << endl;
                }
            }

            break;
        }
        case(2):
            {
                //2D
                //Berechne U, F und G
                for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
                {
                   for(int y = 0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                    {
                        for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                        {
                            cs[k][x][y] = gs.solve_u(raster.get_Zelle(x,y),k);
                            //if(k==1)cout << "u: " << cs[k][x][y] << endl;
                        }

                        for(unsigned int k = 0 ; k < gs.f_u.size() ; k++)
                        {
                            f[k][x][y] = gs.solve_f_u(raster.get_Zelle(x,y),k);
                            //if(k==1)cout << "f: " << f[k][x][y] << endl;
                        }
                        for(unsigned int k = 0 ; k < gs.g_u.size() ; k++)
                        {
                            g[k][x][y] = gs.solve_g_u(raster.get_Zelle(x,y),k);
                            //cout << "g: " << g[k][x][y] << endl;
                            //cin >> blub;
                            //if(k==1)cout << "g: " << g[k][x][y] << endl;
                        }
                    }
                }
                //Berechne Lax-Friedrich Flüsse
                for(int x = 0 ; x < CELLS[0]+ordnung+1 ; x++)
                {
                   for(int y = 0 ; y < CELLS[1]+ordnung+1 ; y++)
                    {
                        for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                        {
                            fi.at(k).at(x).at(y).at(0) = 0.5*(f[k][x][y] + f[k][x+1][y]) + 0.25*(dx/dt)*(cs[k][x][y] - cs[k][x+1][y]);
                            //if(k==0)cout << "fi_x: " << fi.at(k).at(x).at(y).at(0) << endl;

                            fi.at(k).at(x).at(y).at(1) = 0.5*(g[k][x][y] + g[k][x][y+1]) + 0.25*(dy/dt)*(cs[k][x][y] - cs[k][x][y+1]);
                            //if(k==0)cout << "fi_y: " << fi.at(k).at(x).at(y).at(1) << endl;
                        }
                    }
                }

                break;
            }
    }
    return fi;
}
