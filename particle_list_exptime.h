/*Particle_List class for exponential time spacing - stores identifying data for a list of atoms at the start time of each exponential block*/
/*Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#ifndef PARTICLE_LIST_EXPTIME
#define PARTICLE_LIST_EXPTIME

#include "particle_list.h"
#include <iostream>

namespace std{
	
class Particle_List_Exptime: public Particle_List, public Analysis
{
  public:
    Particle_List_Exptime(int n_times = 0, int capacity = 0);
    void atomloop(Analysis*);
    void atomloop(Analysis*, int);
    int show_n_atoms(int timeindex);
    void write_count(string);
    //int show_atom(int timeindex, int atomii);
    int show_n_atoms(int,int);
};
	
}

#endif