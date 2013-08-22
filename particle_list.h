/*Particle_List class - stores identifying data for a list of atoms*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#ifndef PARTICLE_LIST
#define PARTICLE_LIST

#include "analysis.h"
#include <iostream>
#include "system.h"

namespace std{

class Particle_List
{
  protected:
    /*arrays specifying atoms, each containing data for n_atoms atoms at n_times times*/
    int ** species_index;
    int ** molecule_index;
    int ** type_index;
    int ** atom_index;
    int ** trajectoryIDs;
    int n_times;		//number of list times
    int * n_atoms;		//number of atoms at each t
    int capacity;		//atom number capacity of arrays
    int * loopindex;		//index of next atom in external loop
    
  public:
    Particle_List(int n_times = 0, int capacity = 0);
    void clear();
    virtual void atomloop(Analysis*);
    virtual void atomloop(Analysis*,int);
    virtual int show_n_atoms(int timeindex){return n_atoms[timeindex];};
    virtual void write_count(string);
    int show_n_times(){return n_times;};
    virtual int show_atom(int timeii, int atomii);
    
};


}

#endif