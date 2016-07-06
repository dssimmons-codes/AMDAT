// Custom output example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include <iostream>
#include <vector>
#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-3,x_max=3;
const double y_min=-3,y_max=3;
const double z_min=0,z_max=6;

// Set up the number of blocks that the container is divided
// into.
const int n_x=3,n_y=3,n_z=3;

int main() {
    
    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block.
    container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                  false,false,false,8);
    
    // Import the monodisperse test packing and output the Voronoi
    // tessellation in gnuplot and POV-Ray formats.
    //con.import("pack_six_cube");
    
    // Randomly add particles into the container
    if (true) {
        static time_t seed=time(NULL);
        srand((unsigned int)seed);
        seed+=3333;
        for(int i=0;i<=5;i++) {
            double x=x_min+((double)std::rand()/(double)RAND_MAX)*(x_max-x_min);
            double y=y_min+((double)std::rand()/(double)RAND_MAX)*(y_max-y_min);
            double z=z_min+((double)std::rand()/(double)RAND_MAX)*(z_max-z_min);
            //std::cout << i << " " << x << " " << y << " " << z << "\n";
            con.put(i,x,y,z);
        }
        for(int i=7;i<=10;i++) {
            double x=x_min+((double)std::rand()/(double)RAND_MAX)*(x_max-x_min);
            double y=y_min+((double)std::rand()/(double)RAND_MAX)*(y_max-y_min);
            double z=z_min+((double)std::rand()/(double)RAND_MAX)*(z_max-z_min);
            //std::cout << i << " " << x << " " << y << " " << z << "\n";
            con.put(i,x,y,z);
        }
    }
    
    
    // Do a custom output routine to store a variety of face-based
    // statistics. Store the particle ID and position, the number of faces
    // the total face area, the order of each face, the areas of each face,
    // the vertices making up each face, and the neighboring particle (or
    // wall) corresponding to each face.
    con.print_custom("%i %n","packing.custom2");
    
    std::vector<std::vector<int> > neighlist=con.get_neighList();
    for (size_t i=0;i<neighlist.size();++i) {
        for (size_t ii=0;ii<neighlist.at(i).size();++ii) {
            std::cout << neighlist.at(i).at(ii) << " ";
        } std::cout << "\n";
    }
}
