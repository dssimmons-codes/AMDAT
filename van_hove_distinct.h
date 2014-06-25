/*Amorphous Molecular Dynamics Analysis Toolkit (aMDAT)*/
/*Van_Hove_Distinct: a class for calculation of the distinct part of the Van Hove correlation function.*/
/*Written by David S. Simmons*/

#ifndef VAN_HOVE_DISTINCT
#define VAN_HOVE_DISTINCT

#include "space-time_correlation_function.h"
#include "trajectory_list_bins.h"

namespace std{

class Van_Hove_Distinct: public Space_Time_Correlation_Function
{ 
    Trajectory_List_Bins * binned_trajectories;
    
    Trajectory_List * currentlist0;
    Trajectory_List * currentlist1;
    
    bool use_binned;
    
    int nx, ny, nz;		//number of bins in the x, y, and z dimensions

  public:
    Van_Hove_Distinct();
    
    Van_Hove_Distinct(System*sys, Trajectory_List_Bins binnedtraj, int bin_count, float value_max=0);
    Van_Hove_Distinct(System*sys, int bin_count, float value_max=0);
    void set(System*sys, int bin_count, float value_max=0);
    
    
    void analyze(Trajectory_List*);
    void analyze(Trajectory_List*, Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*, int, int, int);
    void listkernel2(Trajectory*, Trajectory*, int, int, int);
    
};

}

#endif