/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class Van_Hove_Self: Class for self part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#ifndef VAN_HOVE_SELF
#define VAN_HOVE_SELF

#include "time_correlation_function.h"
#include "particle_list.h"

namespace std{
	
class Van_Hove_Self: public Time_Correlation_Function
{
    void initialize(System* sys, int bin_count, float value_max, float value_min);
    
    void atom_list(int atomcount, int* species_index, int* molecule_index, int* atom_type, int* atom_index);
    void cellkernel(int, int, int, Atom_Trajectory, Atom_Trajectory){};
    
    /*computational members - just used in trajectory listcalculations; no useful value later on*/
    int timegap, nexttime, currenttime;
    
    void preprocess();
    void preprocess_list(Particle_List* particle_list);
    
    
  public:
    Van_Hove_Self(); 
    Van_Hove_Self(System* sys, int bin_count, float value_max=0, float value_min=0);
    void set(System* sys, int bin_count, float value_max=0, float value_min=0);
    
    void atomlist(Particle_List * particle_list){preprocess_list(particle_list);system->displacement_loop_list(this, particle_list);postprocess_list(particle_list);};
    
    
    void displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj);
    void atomkernel(Trajectory * traj);
    
 
    void list_displacementkernel(int timegapii, int thisii, int nextii, Particle_List* particle_list);
    void postprocess();
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel (Trajectory* current_trajectory);

};

}

#endif