#include "FORCE.h"

using namespace std;

FORCE::FORCE(std::string const_in, std::string formel_in, std::string winkel_in, std::string kreis_in, std::string save_in) : numerische_methode("FORCE", const_in, formel_in, winkel_in, kreis_in, save_in)
{

}

vector< vector< vector< vector <double> > > > FORCE::calc_method_flux()
{
    cout << "Berechne FORCE Fluss..." << endl;

    double cs[gs.u.size()][raster.getwidth()][raster.getheight()];
    double fd[gs.u.size()][raster.getwidth()][raster.getheight()];
    double gd[gs.u.size()][raster.getwidth()][raster.getheight()];
    double f_lax[gs.u.size()][raster.getwidth()][raster.getheight()];
    double f_rie[gs.u.size()][raster.getwidth()][raster.getheight()];
    double g_lax[gs.u.size()][raster.getwidth()][raster.getheight()];
    double g_rie[gs.u.size()][raster.getwidth()][raster.getheight()];

    vector< vector < vector< vector<double> > > > f_force (gs.u.size(), vector< vector< vector<double> > > ( raster.getwidth(), vector< vector<double> > ( raster.getheight() , vector<double> (dimension , 0.0) ) ) );

    switch(raster.getdim())
    {
    case(1):
        {
            Raster u_rie(raster.getwidth());
            for(int i = 0 ; i < CELLS[0]+2*ordnung+1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    cs[k][i][0] = gs.solve_u(raster.get_Zelle(i),k);
                }

                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    fd[k][i][0] = gs.solve_f_u(raster.get_Zelle(i),k);
                }
            }
            //LaxFriedrichFluss berechnen
            for(int i = 0 ; i < CELLS[0] +ordnung+1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    f_lax[k][i][0] = 0.5*(fd[k][i][0] + fd[k][i+1][0]) + 0.5*(dx/dt)*(cs[k][i][0] - cs[k][i+1][0]);
                }
            }
            //Richtmyer Fluss berechnen
            for(int i = 0 ; i < CELLS[0] +ordnung +1 ; i++)
            {
                u_rie.set_Zelle_d(0.5*(cs[0][i][0] + cs[0][i+1][0]) + 0.5*(dt/dx)*(fd[0][i][0] - fd[0][i+1][0]),i);
                u_rie.set_Zelle_ux((0.5*(cs[1][i][0] + cs[1][i+1][0]) + 0.5*(dt/dx)*(fd[1][i][0] - fd[1][i+1][0]))/u_rie.get_Zelle(i).d,i);
                u_rie.set_Zelle_uxr(0.5*(cs[2][i][0] + cs[2][i+1][0]) + 0.5*(dt/dx)*(fd[2][i][0] - fd[2][i+1][0]),i);
                u_rie.set_Zelle_p(ct*pow(u_rie.get_Zelle(i).d,g),i);
            }
            for(int i = 0 ; i < CELLS[0]+2*ordnung+1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    f_rie[k][i][0] = gs.solve_f_u(u_rie.get_Zelle(i),k);
                }
            }

            //FORCE Fluss berechnen
            for(int i = 0 ; i < CELLS[0] +ordnung +1 ; i++)
            {
                for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                {
                    f_force.at(k).at(i).at(0).at(0) = 0.5*(f_lax[k][i][0] + f_rie[k][i][0]);
                }
            }
            break;
        }
    case(2):
        {
            //2D
            Raster u_rie_f(raster.getwidth(), raster.getheight());
            Raster u_rie_g(raster.getwidth(), raster.getheight());

            for(int x = 0 ; x < CELLS[0]+2*ordnung+1 ; x++)
            {
                for(int y = 0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                {
                    for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                    {
                        cs[k][x][y] = gs.solve_u(raster.get_Zelle(x,y),k);
                        //cout << "u: " << cs[k][x][y] << endl;
                    }
                    for(unsigned int k = 0 ; k < gs.f_u.size() ; k++)
                    {
                        fd[k][x][y] = gs.solve_f_u(raster.get_Zelle(x,y),k);
                        //cout << "f: " << fd[k][x][y] << endl;
                    }
                    for(unsigned int k = 0 ; k < gs.g_u.size() ; k++)
                    {
                        gd[k][x][y] = gs.solve_g_u(raster.get_Zelle(x,y),k);
                        //cout << "g: " << gd[k][x][y] << endl;
                        //cin >> g_rie[0][0][0];
                    }
                }
            }

            //Lax-Friedrich Fluss berechnen

            for(int x = 0 ; x < CELLS[0]+ordnung+1 ; x++)
            {
                for(int y = 0 ; y < CELLS[1]+ordnung+1 ; y++)
                {
                    for(unsigned int k = 0 ; k < gs.u.size() ; k++)
                    {
                        f_lax[k][x][y] = 0.5*(fd[k][x][y] + fd[k][x+1][y]) + 0.25*(dx/dt)*(cs[k][x][y] - cs[k][x+1][y]);
                        //if(k==0)cout << "fi_x: " << fi.at(k).at(x).at(y).at(0) << endl;

                        g_lax[k][x][y] = 0.5*(gd[k][x][y] + gd[k][x][y+1]) + 0.25*(dy/dt)*(cs[k][x][y] - cs[k][x][y+1]);
                        //if(k==0)cout << "fi_y: " << fi.at(k).at(x).at(y).at(1) << endl;
                    }
                }
            }

            //Richtmyer Fluss berechnen

            for(int x=0 ; x < CELLS[0]+ordnung+1 ; x++)
            {
                for(int y=0 ; y < CELLS[1]+ordnung+1 ; y++)
                {
                    u_rie_f.set_Zelle_d((0.5*(cs[0][x][y]+cs[0][x+1][y])+(dt/2*dx)*(fd[0][x][y]-fd[0][x+1][y])),x,y);
                    u_rie_f.set_Zelle_ux((0.5*(cs[1][x][y]+cs[1][x+1][y])+(dt/2*dx)*(fd[1][x][y]-fd[1][x+1][y]))/u_rie_f.get_Zelle(x,y).d,x,y);
                    u_rie_f.set_Zelle_uy((0.5*(cs[2][x][y]+cs[2][x+1][y])+(dt/2*dx)*(fd[2][x][y]-fd[2][x+1][y]))/u_rie_f.get_Zelle(x,y).d,x,y);
                    u_rie_f.set_Zelle_uxr((0.5*(cs[3][x][y]+cs[3][x+1][y])+(dt/2*dx)*(fd[3][x][y]-fd[3][x+1][y])),x,y);
                    u_rie_f.set_Zelle_uyr((0.5*(cs[4][x][y]+cs[4][x+1][y])+(dt/2*dx)*(fd[4][x][y]-fd[4][x+1][y])),x,y);
                    u_rie_f.set_Zelle_p(ct*pow(u_rie_f.get_Zelle(x,y).d,g),x,y);

                    u_rie_g.set_Zelle_d((0.5*(cs[0][x][y]+cs[0][x][y+1])+(dt/2*dy)*(gd[0][x][y]-gd[0][x][y+1])),x,y);
                    u_rie_g.set_Zelle_ux((0.5*(cs[1][x][y]+cs[1][x][y+1])+(dt/2*dy)*(gd[1][x][y]-gd[1][x][y+1]))/u_rie_g.get_Zelle(x,y).d,x,y);
                    u_rie_g.set_Zelle_uy((0.5*(cs[2][x][y]+cs[2][x][y+1])+(dt/2*dy)*(gd[2][x][y]-gd[2][x][y+1]))/u_rie_g.get_Zelle(x,y).d,x,y);
                    u_rie_g.set_Zelle_uxr((0.5*(cs[3][x][y]+cs[3][x][y+1])+(dt/2*dy)*(gd[3][x][y]-gd[3][x][y+1])),x,y);
                    u_rie_g.set_Zelle_uyr((0.5*(cs[4][x][y]+cs[4][x][y+1])+(dt/2*dy)*(gd[4][x][y]-gd[4][x][y+1])),x,y);
                    u_rie_g.set_Zelle_p(ct*pow(u_rie_g.get_Zelle(x,y).d,g),x,y);
                }
            }

            for(int x=0 ; x < CELLS[0]+2*ordnung+1 ; x++)
            {
                for(int y=0 ; y < CELLS[1]+2*ordnung+1 ; y++)
                {
                    for(unsigned int k=0 ; k < gs.u.size() ; k++)
                    {
                        f_rie[k][x][y] = gs.solve_f_u(u_rie_f.get_Zelle(x,y),k);
                        g_rie[k][x][y] = gs.solve_g_u(u_rie_g.get_Zelle(x,y),k);
                    }
                }
            }

            //FORCE Fluss berechnen

            for(int x=0 ; x < CELLS[0]+ordnung+1 ; x++)
            {
                for(int y=0 ; y < CELLS[1]+ordnung+1 ; y++)
                {
                    for(unsigned int k=0 ; k < gs.u.size() ; k++)
                    {
                        f_force.at(k).at(x).at(y).at(0) = 0.5*(f_lax[k][x][y] + f_rie[k][x][y]);
                        f_force.at(k).at(x).at(y).at(1) = 0.5*(g_lax[k][x][y] + g_rie[k][x][y]);
                    }
                }
            }

            break;
        }
    }
    return f_force;
}
