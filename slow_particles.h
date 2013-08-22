/*Slow_Particles class: identifies and stores list of slow trajectories at start of each exponential block*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#ifndef SLOW_PARTICLES
#define SLOW_PARTICLES

#include "analysis.h"
#include "exptime_trajectory_list.h"
#include "system.h"
#include "gaussian_comparison.h"

namespace std{

class Slow_Particles: public Exptime_Trajectory_List
{
  private:
    Gaussian_Comparison * gaussian_comparison;
    int displacement_time_index;
    float maxdistance;
    int atoms_considered;		//number of atoms tested

  public:
    Slow_Particles();
    Slow_Particles(System * sys, Gaussian_Comparison * gc);
    void set(System *, Gaussian_Comparison *);
    
    void displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj);
    void atomkernel(Trajectory * traj);
};

}

#endif