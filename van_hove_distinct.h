/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Van_Hove_Distinct: a class for distinct part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#ifndef VAN_HOVE_DISTINCT
#define VAN_HOVE_DISTINCT

#include "time_correlation_function.h"
#include "spacial_decomposition.h"
#include "particle_list.h"

namespace std{

class Van_Hove_Distinct: public Time_Correlation_Function
{
    int min_cell_size;
    int max_cell_range;
    Spacial_Decomposition const * cellobject;
    void initialize(System*sys, int bin_count, const Spacial_Decomposition* cells, float min_size, float value_max, float value_min);

    void atom_list(int atomcount, int* species_index, int* molecule_index, int* atom_type, int* atom_index);
    
    /*computational members - just used in particle listcalculations; no useful value later on*/
    int nexttime;
    int timegap;
    
  public:
    Van_Hove_Distinct();
    Van_Hove_Distinct(System*sys, int bin_count, const Spacial_Decomposition* cells, float min_size=0, float value_max=0, float value_min=0);
    void set(System*sys, int bin_count, const Spacial_Decomposition* cells, float min_size=0, float value_max=0, float value_min=0){initialize(sys, bin_count, cells, min_size, value_max, value_min);};
    
    void displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj);
    void atomkernel(Trajectory * traj);
    
    void single_atom(int species_index, int molecule_index, int atom_type, int atom_index);
    void cellkernel(int timegap, Coordinate coordinate1, Coordinate coordinate2, int atom_type=0);
    void postprocess();
    
    void atomlist_kernel(int, int, int, int, int, Particle_List*);
    void list_displacementkernel(int timegapii, int thisii, int nextii, Particle_List* particle_list);
    void atomlist(Particle_List * particle_list){system->displacement_loop_list(this, particle_list);postprocess_list(particle_list);};
    
};

}

#endif