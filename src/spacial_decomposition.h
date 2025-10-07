// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Spacial Decomposition class declaration: class to box atoms into spacial cells by position over a range of time*/
/*Written by David S. Simmons*/

#ifndef SPACIAL_DECOMPOSITION
#define SPACIAL_DECOMPOSITION

#include "atom_trajectory.h"
#include "system.h"

namespace std{

class Spacial_Decomposition: public Analysis
{

    int ** atomID;			//list of unique atom IDs corresponding to atoms in each cell at each time
    Coordinate ** spacial_cells;
    
    int n_timesteps;
    int n_cells[3];			//number of cells in x,y,and z dimensions, counting overhand cells
    int n_real_cells[3];		//number of non-overhang cells (cells in box bondaries) in x, y, and z directions
    int * atom_count;			//number of atoms in each cell at each time
    Coordinate cell_size;		//spatial size of cells
    
    /*box boundaries*/
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;
    
    int * cell_capacity;		//currently allocated storage capacity for each cell
    int total_cells;			//total number of cells
    int countmax;			//the most atoms in any cell
    int overhang;			//number of overhang cells on each side

    
    void initialize(System* sys, float min_box_size, int cell_capacity, int over);	//initializer; min_box_size is a user specification for the minimum size of the box; the method will find the box size in each dimension that is an even divisor of the simulation size that is at least equal to this.

    void displacementkernel(int, int, int ,int, int, int, int){};
    void cellkernel(int,Coordinate,Coordinate){};
    void change_capacity(int index, int capacity);
    
  public:
    Spacial_Decomposition(System* sys, float min_cell_size, int cell_capacity, int over);
    ~Spacial_Decomposition();
    void clear_memory();
    void atomkernel(Trajectory * traj);
    
    void atomlist_kernel(int species_index, int molecule_index, int atom_type, int atom_index);	
    void loop_cell(Analysis* analysis, int timegap, int nextii, Coordinate coordinate1, int atom1ID, int deltax, int deltay, int deltaz)const;
    void loop_cell_ID(Analysis * analysis, int thisii, int nextii, int timegapii, Coordinate coordinate1, int atom1ID, int deltax, int deltay, int deltaz)const;
    int convert_index(int timeii, int xii, int yii, int zii)const;
    void fill_overhang();
    int show_overhang()const{return overhang;};
    
    void preprocess_list(Particle_List*){};
};

}

#endif