/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class Van_Hove_Self: Class for self part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#ifndef VAN_HOVE_SELF
#define VAN_HOVE_SELF

#include "space-time_correlation_function.h"

namespace std{
	
class Van_Hove_Self: public Space_Time_Correlation_Function
{
    void initialize(System* sys, int bin_count, float value_max);
    
    //void atom_list(int atomcount, int* species_index, int* molecule_index, int* atom_type, int* atom_index);


  public:
    Van_Hove_Self(); 
    Van_Hove_Self(System* sys, int bin_count, float value_max=0);
    void set(System* sys, int bin_count, float value_max=0);
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel (Trajectory* current_trajectory, int timegapii, int thisii, int nextii);

};

}

#endif