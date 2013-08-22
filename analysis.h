/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Parent analysis class. (Abstract class)*/
/*Written by David S. Simmons*/

/*All analytical classes that wish to use the data structure's built-in looping methods for different data sets should inherit from this class.  The 'kernel' methods should be overloaded in each daughter class to be passed to the data structure loop methods.*/

/*The displacementkernel method is for use with analytical methods that need to consider displacements in time.  Typically its first three arguments will be the displacement time and two timesteps yielding that displacement time.  Its next four will be the same arguments as those for the atomkernel method.*/

/*The cellkernel method is for use with loops over cells in methods that refer to a spacial decomposition of the atom positions.  It receives the timegap and two coordinates from a loop over a cell within a spacial_decomposition object.*/

/*The setkernel method is for use any time some set of atom or molecules must be looped over.  Its four arguments are typically species index, molecule index, atom type, and atom index.  These are enough to identify any individual atom or molecule within the system.*/

/*The postprocess method automatically runs after the system loop is complete, and is typically intended for such tasks as data normalization.*/

/*This class is defined to do analysis on a single system; some analytical techniques performed on multiple systems can be built by using an array of Analysis objects, one for each system included in the analysis*/

#ifndef ANALYSIS
#define ANALYSIS

#include <string>
#include <stdio.h>
#include "trajectory_list.h"

namespace std{
	
class System;
class Particle_List;
	
class Analysis
{
  protected:
    System * system;	//the system on which analysis is to be performed
    
    Trajectory_List * trajectory_list;			//Array of trajectory_lists used by this analysis tool
    
    virtual void preprocess(){};
    virtual void preprocess_list(Particle_List*){};
    virtual void postprocess(){};		//method to automatically run after loop for postprocessing
    virtual void postprocess_list(Particle_List*){};		//method to automatically run after loop over time-varying atom list for postprocessing
    virtual void postprocess_list(){};		//method to automatically run after loop over time-varying trajectory list for postprocessing
    
  public:
    
    /*Constructors and assignment*/
    
    Analysis();				//default constructor
    Analysis(const Analysis &);		//copy constructor
    Analysis operator = (const Analysis &);	//assignment
    
    virtual void analyze(Trajectory_List *){cout<<"Error: Trajectory list targets not implemented for this analysis method.\n";};
    
    /*System loop methods*/
    virtual void atomkernel(Trajectory*){cout<<"Error: System set targets not implemented for this analysis method.\n";};  //new
    virtual void displacementkernel(int,int,int,Trajectory*){}; //new
    
    /*Spatial decomposition loop methods*/
    virtual void cellkernel(int, Coordinate, Coordinate,int atom_type=0){cout<<"Error: Cell loops not implemented for this analysis method.\n";};
    virtual void cellkernelID(int, int, int, int, int, Coordinate, Coordinate){cout<<"Error: Cell loops not implemented for this analysis method.\n";};
    
    /*Trajectory list loop methods*/
    virtual void list_displacementkernel(int, int, int){cout<<"Error: Trajectory list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Trajectory*){cout<<"Error: Trajectory list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Trajectory*, int, int, int){cout<<"Error: Trajectory list targets not fully implemented for this analysis method.\n";};
    
    /*Methods to use with trajectory list bins*/
    virtual void bin_hook(Trajectory_List *,int,int,int){cout<<"Error: Trajectory list bins not implemented for this analysis method.\n";};
    virtual void postprocess_bins(){};
    
    /*Particle list methods - outdated and slated for removal*/
    //virtual void list_displacementkernel(int, int, int, Particle_List*){};
    virtual void atomlist_kernel(int,int,int,int,int,Particle_List*){};
    
    
    virtual void write(string)const{cout<<"Error: No standard write method implemented for this analysis method.\n";};		//generic method for writing results to file
    
    /*Methods for employing the loops implemented in class System over various sets of atoms*/
    virtual void atomlist(Particle_List *);
    void all_moleculecom();
    int species_moleculecom(){return 1;};
    void species_moleculecom(int species_index);
    void species_moleculecom(string species_name);
    int atom_species(){return 3;};
    void atom_species(int species_index, int atom_type, int atom_index);
    void atom_species(string species_name, string atomtype_name, int atom_index);
    int type_molecule(){return 3;};
    void type_molecule(int species_index, int molecule_index, int atom_type);	
    void type_molecule(string species_name, int molecule_index, string atomtype_name);
    int molecule(){return 2;};
    void molecule(int species_index, int molecule_index);
    void molecule(string species_name, int molecule_index);
    int species(){return 1;};
    void species(int species_index);
    void species(string species_name);
    int type_species(){return 2;};
    void type_species(int species_index, int atom_type);
    void type_species(string species_name, string atomtype_name);
    int type_system(){return 1;};
    void type_system(int atom_type);
    void type_system(string atomtype_name);
    void all();

    void analyze(string runline);
    
    System* show_system(){return system;};		//method to return system on which analysis was performed
	virtual bool isThreadSafe() {return false;};
};
}

#endif
