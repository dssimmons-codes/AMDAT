/*Molecular Dynamics Analysis Toolkit*/
/*Class to calculate non-gaussian parameter*/
/*Written by David S. Simmons*/

#ifndef NON_GAUSSIAN_PARAMETER
#define NON_GAUSSIAN_PARAMETER

#include "system.h"
#include <string.h>
#include "mean_square_displacement.h"

namespace std{

class Non_Gaussian_Parameter: public Analysis
{
    int n_times;
    float * ngp;
    int * weighting;
    Mean_Square_Displacement const * msd;
    float * timetable;
    int atomcount;
    float * n_atoms;
    void postprocess();
    
    //calculation variables
    int currenttime, nexttime, currenttimegap;
    
  public:
    Non_Gaussian_Parameter();
    Non_Gaussian_Parameter(const Non_Gaussian_Parameter &);
    Non_Gaussian_Parameter(System*, const Mean_Square_Displacement*);
    Non_Gaussian_Parameter operator=(const Non_Gaussian_Parameter &);
    
    void displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj);
    void atomkernel(Trajectory * traj);
    
    void write(string)const;
    int max()const;

    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void postprocess_list();
    
    void bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii);
    void postprocess_bins();
};

}

#endif