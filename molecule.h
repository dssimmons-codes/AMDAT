/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to hold table of atom trajectories in a given molecule.  Right now there is not clearly a reason to have this, but it exists because it is the correct scope at which to implement a method to calculate the movement of the center of mass of a molecule, for example.  This may be desirable later.  It also permits a heirarchal storage of sysyem data in which the form of the hierarchy specifies the nature of the system.*/
/*Written by David S. Simmons*/

#ifndef MOLECULE
#define MOLECULE


#include "atom_trajectory.h"
#include "coordinate.h"
#include "trajectory.h"

namespace std {

class Molecule:public Trajectory
{
    int total_atoms;
    int n_atomtypes;
    int * n_atoms;
    
    int moleculeID;
    int species;		//stores species of molecule
    
    bool atoms_unwrapped;
    
    Atom_Trajectory** atoms;
    void clear_memory();			//method to clear memory allocated to arrays so that they may be recreated
    void calculate_mass();
    
  public:
    Molecule(int types=0,int ts=0);
    Molecule(int,int*,int);
    ~Molecule();
    void set(int, int *,int);
    int typecount(int);				//method to return number of atoms of a given type
    int atomcount();				//method to return number of atoms of a given type
    float * msd(int, int, int);		//method to pass on msd of a particular atom
    Atom_Trajectory show_atom(int,int);	//method to return atom_trajectory object in molecule
    int unwrap_atoms(const Coordinate &);			//unwrap all atom trajectories (only works w/ constant volume)
    int wrap_atoms(const Coordinate *, Coordinate **);				//wrap all atom trajectories				
    void calculate_center_of_mass(const Coordinate *, Coordinate**);
    void calculate_center_of_mass();
    
    Coordinate show_unwrapped(int atom_type, int atom_index, int timestep)const;	//show unwrapped coordinate of a given atom at a given time
    
    Atom_Trajectory* show_atom_trajectory(int, int);
    float atom_distance(int,int,int,int);
    
    float gyration_radius();
    
    void set_coordinate(int, int, const Coordinate &, int);
    void set_unwrapped(int, int, const Coordinate &, int);
    void set_velocity(int, int, const Coordinate &, int);
    void set_moleculeID(int ID){moleculeID=ID;};	//set unique molecule ID
    void set_species(int spec){species=spec;};		//set species
    int show_species()const{return species;};		//return species index
    
    void ID_to_atoms();		//Pass unique ID of molecule down to constituent atoms
};

}

#endif