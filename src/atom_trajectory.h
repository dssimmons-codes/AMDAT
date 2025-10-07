// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to store atom data and trajectory*/
/*Written by David S. Simmons*/

#ifndef ATOM_TRAJECTORY
#define ATOM_TRAJECTORY

#include "trajectory.h"


namespace std {



/*-------------------------*/
/*--Atom Trajectory Class--*/
/*-------------------------*/

/*Class to handle an entire atom trajectory*/
class Atom_Trajectory: public Trajectory
{    
    int atomID;							//unique ID of atom
    int moleculeID;						//unique ID of molecule within which this atom resides
    
  public:
    using Trajectory::set;		//no idea what this is but the code doesn't compile properly with it. Who added it and why?
    Atom_Trajectory();		//default constructor
    Atom_Trajectory(int timesteps);				//constructor to instantiate coordinate array without defining it
    Atom_Trajectory(int,Coordinate *);				//constructor to instantiate coordinate array with full definition
    void set(int,int,int m=1);					//change number of timesteps and reallocate memory accordingly
    void set(int,int,Coordinate*,int);				//method to fully define object, including coordinate list
        
    int show_atomID()const{return atomID;};
    int show_moleculeID()const{return moleculeID;};
    void set_atomID(int ID){atomID = ID;};
    
    void set_moleculeID(int ID){moleculeID=ID;};		//set unique ID of molecule within which atome resides
};

}

#endif
