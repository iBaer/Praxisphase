/*
 * cluster_builder.cpp
 *
 *  Created on: 25.07.2016
 *      Author: pascal
 */

#include "cluster_builder.h"
#include <sstream>
#include <set>
#include "grid.h"

using namespace std;

Cluster_Builder::Cluster_Builder() {

	marked_cells = nullptr;

}

Cluster_Builder::~Cluster_Builder() {

}

void Cluster_Builder::build_clusters(Grid* grid, int* &marked_cells, vector<Cluster_Square>& clusters){

	this->marked_cells = marked_cells;
	this->grid_main = grid;

	// 3. Feineres Gitter erzeugen
	int cluster_amount = binary_clustering(marked_cells, grid_main->grid_size_total[0], grid_main->grid_size_total[1]);
	//vector<Cluster_Square> clusters;	//TODO: Grid anhängen, das feinere Gitter (Cluster) haben kann
	// TODO: Grid Level speichern!?
	int with_adjacencies;
	// Initiales Rechteck Cluster bilden
	square_clustering(clusters, marked_cells, grid_main, cluster_amount);
	// Bilde größere Rechtecke solange sich Cluster berühren
	if (cluster_amount > 1) {
		do {
			with_adjacencies = cluster_adjacency_check(clusters, marked_cells, cluster_amount);
			cout << " ---> Adjencencies found: " << with_adjacencies << endl;
			if (with_adjacencies != 0 && cluster_amount > 1)
				square_cluster_merge(clusters, marked_cells, grid_main, cluster_amount);
			else
				break;
		} while (with_adjacencies != 0 || cluster_amount <= 1);
	}
}

// Start-Algorithmus zur Suche von Clustern
int Cluster_Builder::binary_clustering(int* &marked_cells, int x_max, int y_max){

	int pos = 0;
	int nth_cluster = 0;
	int* grid_size = new int[2];
	grid_size[0] = x_max;
	grid_size[1] = y_max;
	for(int y=0;y<y_max;y++){
		for(int x=0;x<x_max;x++){
			pos = x + y * grid_size[0];
			if(marked_cells[pos]==-1){
				//cout << x << " " << y << " " << marked_cells[pos] << " pos " << pos <<endl;
				nth_cluster = nth_cluster + 1;
				clusterize(marked_cells, grid_size, nth_cluster, x, y);
			}
		}
	}


	//TODO: Remove
	 ofstream myfile ("amr_cluster.txt");
	  if (myfile.is_open())
	  {
		for (int y = 0; y < grid_size[1]; y++) {
			for (int x = 0; x < grid_size[0]; x++) {
				pos = x + y * grid_size[0];
				myfile << x << " " << y << " " << marked_cells[pos] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";

		delete[] grid_size;
		//delete[] clustered_cells;
	return nth_cluster;
}

//TODO: Gefahr wegen Rekursion?
// Rekursiver "Cluster" Algorithmus
void Cluster_Builder::clusterize(int* &clustered_cells, int* grid_size, int nth_cluster, int x, int y){
	int pos = x + y * grid_size[0];
	if(x>= 0 && x< grid_size[0] && y>= 0 && y< grid_size[1] && clustered_cells[pos]==-1){
		//cout << "Marking at " << x << "," << y << " | old: " << clustered_cells[pos] << " | new: " << nth_cluster << endl;
		clustered_cells[pos] = nth_cluster;
		clusterize(clustered_cells, grid_size, nth_cluster, x-1, y);
		clusterize(clustered_cells, grid_size, nth_cluster, x, y-1);
		clusterize(clustered_cells, grid_size, nth_cluster, x+1, y);
		clusterize(clustered_cells, grid_size, nth_cluster, x, y+1);
	}
}

void Cluster_Builder::square_clustering(vector<Cluster_Square> &clusters, int* &marked_cells, Grid* grid, int cluster_amount){
	int pos = 0;
	int x_max = grid->grid_size_total[0];
	int y_max = grid->grid_size_total[1];

	for(int n = 1; n <= cluster_amount; n++){
		int first = 1;
		clusters.push_back(Cluster_Square());
		for (int y = 0; y < y_max; y++) {
			for (int x = 0; x < x_max; x++) {

				pos = x + y * x_max;

				if(marked_cells[pos]==n){
					if(first==1){
						clusters.at(n-1).pos_x_min = x;
						clusters.at(n-1).pos_x_max = x;
						clusters.at(n-1).pos_y_min = y;
						clusters.at(n-1).pos_y_max = y;
						first=0;
					}
					else if(x < clusters.at(n-1).pos_x_min)
						clusters.at(n-1).pos_x_min = x;
					else if(x > clusters.at(n-1).pos_x_max)
						clusters.at(n-1).pos_x_max = x;
					else if(y < clusters.at(n-1).pos_y_min)
						clusters.at(n-1).pos_y_min = y;
					else if(y > clusters.at(n-1).pos_y_max)
						clusters.at(n-1).pos_y_max = y;
				}
			}
		}
	}

	// Marking new clusters
	for(int n = 0; n < cluster_amount; n++){
		cout << "Cluster #"<<n+1<<" from "<< clusters.at(n).pos_x_min<<","<<clusters.at(n).pos_y_min<<" to "<< clusters.at(n).pos_x_max<<","<<clusters.at(n).pos_y_max<<endl;
		for (int y = clusters.at(n).pos_y_min; y <= clusters.at(n).pos_y_max; y++) {
			for (int x = clusters.at(n).pos_x_min; x <= clusters.at(n).pos_x_max; x++) {
				pos = x + y * x_max;
				marked_cells[pos] = n+1;
			}
		}
	}

	//TODO: Remove
	 ofstream myfile ("amr_cluster_squared.txt");
	  if (myfile.is_open())
	  {
		for (int y = 0; y < y_max; y++) {
			for (int x = 0; x < x_max; x++) {
				pos = x + y * x_max;
				myfile << x << " " << y << " " << marked_cells[pos] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";

	  return;

}

int Cluster_Builder::cluster_adjacency_check(vector<Cluster_Square> &clusters, int* &marked_cells, int cluster_amount) {

	int with_adjacencies = 0;

	for (int n = 0; n < cluster_amount; n++) {
		for (int l = n + 1; l < cluster_amount; l++) {
			if (clusters[n].pos_x_min - 1 <= clusters[l].pos_x_max && clusters[n].pos_x_max + 1 >= clusters[l].pos_x_min) {
				if (clusters[n].pos_y_min - 1 <= clusters[l].pos_y_max && clusters[n].pos_y_max + 1 >= clusters[l].pos_y_min) {
					clusters[n].adjacencies.insert(l);
					with_adjacencies++;
				}
			}
		}
	}

	//TODO: Output Optional
	for (int n = 0; n < cluster_amount; n++) {
		cout << "Cluster #" << n + 1 << " adjacent to: ";
		if(clusters[n].adjacencies.begin() == clusters[n].adjacencies.end())
			cout << "NONE";
		for (set<int>::iterator it = clusters[n].adjacencies.begin(); it != clusters[n].adjacencies.end(); ++it){
			if(it != clusters[n].adjacencies.begin())
				cout << ", ";
			cout << "Cluster #" << *it+1;
		}
		cout << endl;
	}

	return with_adjacencies;
}

// Vereint benachbarte Cluster zu einem größeren ! Rechteck !
void Cluster_Builder::square_cluster_merge(vector<Cluster_Square> &clusters, int* &marked_cells, Grid* grid, int &cluster_amount){

	int* cluster_to_destroy = new int[cluster_amount]();
	int x_max = grid->grid_size_total[0];
	int y_max = grid->grid_size_total[1];

	// Fügt Clustern komplette Liste mit Nachbarn zu, damit auch nicht direkt verbundene Nachbarn dabei sind
	for (int n = cluster_amount-1; n >=0; n--) {
		if(clusters[n].adjacencies.begin() != clusters[n].adjacencies.end()){
			for (set<int>::iterator it_target = clusters[n].adjacencies.begin(); it_target != clusters[n].adjacencies.end(); ++it_target){
				if(clusters[*it_target].adjacencies.begin() != clusters[*it_target].adjacencies.end()){
					for (set<int>::iterator it_source = clusters[*it_target].adjacencies.begin(); it_source != clusters[*it_target].adjacencies.end(); ++it_source){
						clusters[n].adjacencies.insert(*it_source);
					}
				}
			}
		}
	}


	for(int n = 0; n < cluster_amount; n++){
		// Zerstöre alle zusammenhängenden Cluster mit höherem Index
		for (set<int>::iterator it = clusters[n].adjacencies.begin(); it != clusters[n].adjacencies.end(); ++it){
			if(*it>n){
				cluster_to_destroy[*it] = 1;
			}
		}
	}

	//TODO: Output Optional
	for (int n = 0; n < cluster_amount; n++) {
		cout << "MergeCluster #" << n + 1 << " adjacent to: ";
		if(clusters[n].adjacencies.begin() == clusters[n].adjacencies.end())
			cout << "NONE";
		for (set<int>::iterator it = clusters[n].adjacencies.begin(); it != clusters[n].adjacencies.end(); ++it){
			if(it != clusters[n].adjacencies.begin())
				cout << ", ";
			cout << "Cluster #" << *it+1;
		}
		cout << endl;
	}


	// Verbindet Cluster in dem es min/max Positionen anpasst
	for (int n = 0; n < cluster_amount; n++) {
		for (set<int>::iterator it = clusters[n].adjacencies.begin(); it != clusters[n].adjacencies.end(); ++it){
			if(clusters[n].pos_x_min > clusters[*it].pos_x_min)
				clusters[n].pos_x_min = clusters[*it].pos_x_min;
			if(clusters[n].pos_x_max < clusters[*it].pos_x_max)
				clusters[n].pos_x_max = clusters[*it].pos_x_max;
			if(clusters[n].pos_y_min > clusters[*it].pos_y_min)
				clusters[n].pos_y_min = clusters[*it].pos_y_min;
			if(clusters[n].pos_y_max < clusters[*it].pos_y_max)
				clusters[n].pos_y_max = clusters[*it].pos_y_max;
		}
		clusters[n].adjacencies.clear();
	}

	// Löscht Cluster, die mit anderen "verbunden" sind
	for (int n = cluster_amount-1; n >=0; n--) {
		if(cluster_to_destroy[n] == 1){
			clusters.erase(clusters.begin()+n);
			cluster_amount--;
		}
	}
	delete[] cluster_to_destroy;

	int pos = 0;
	// Markierierungen mit Indices der kombinierten Cluster
	for(int n = 0; n < cluster_amount; n++){
		cout << "Cluster #"<<n+1<<" from "<< clusters.at(n).pos_x_min<<","<<clusters.at(n).pos_y_min<<" to "<< clusters.at(n).pos_x_max<<","<<clusters.at(n).pos_y_max<<endl;
		for (int y = clusters.at(n).pos_y_min; y <= clusters.at(n).pos_y_max; y++) {
			for (int x = clusters.at(n).pos_x_min; x <= clusters.at(n).pos_x_max; x++) {
				pos = x + y * x_max;
				marked_cells[pos] = n+1;
			}
		}
	}

	//TODO: Remove
	 ofstream myfile ("amr_cluster_combined.txt");
	  if (myfile.is_open())
	  {
		for (int y = 0; y < y_max; y++) {
			for (int x = 0; x < x_max; x++) {
				pos = x + y * x_max;
				myfile << x << " " << y << " " << marked_cells[pos] << endl;
			}
		}
	    myfile.close();
	  }
	  else cout << "Unable to open file";
}
