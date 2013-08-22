/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*RgTensor_Stats class: calculates time-dependent Rg Tensor statistics for a set of trajectories*/
/*Written by David S. Simmons*/

#include "analysis.h"
#include "rgtensor.h"

#ifndef RGTENSOR_STATS
#define RGTENSOR_STATS

namespace std{

class RgTensor_Stats:public Analysis
{
  float * mean_gyration_radius;
  threefloat * eigenvalues; // stores mean eigenvalues as a function of time
  threefloat ** single_eigenvalues; // stores single-trajectory eigenvalues as a function of time;
  float ** single_relative_asphericity; // stores single-trajectory relative asphericity as a function of time;
  float * mean_relative_asphericity;
  float ** single_gyration_radius; // stores single-trajectory gyration radius as a function of time;
  float ** rel_asphericity_dist;
  float ** gyration_rad_dist;
  
  int n_times;
  int n_trajectories;
  void atomkernel(int,int,int,int){};
  void displacementkernel(int,int,int,int,int,int,int){};
  
  bool rel_asphericity_dist_calc;
  bool gyration_rad_dist_calc;
  int rel_asphericity_bins;
  int gyration_rad_bins;
  float max_gyr_rad_binned;
  
  
  /*dummy variables*/
  int traj_index;	//variable to keep track of which member of the list we are on
  
  
  public:
    RgTensor_Stats(System*);
    
    void analyze (Trajectory_List*);
    void listkernel(Trajectory*);
    
    void write(string);
    void write_rel_asphericity_dist(string);
    void write_gyration_rad_dist(string);
    
    void calc_rel_asphericity_dist(int nbins);
    void calc_gyration_rad_dist(int,float);  
    
};

}

#endif