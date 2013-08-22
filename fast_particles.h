/*Fast_Particles class: identifies and stores list of fast particles at start of each exponential block*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#ifndef FAST_PARTICLES
#define FAST_PARTICLES

#include "analysis.h"
#include "exptime_trajectory_list.h"
#include "system.h"
#include "gaussian_comparison.h"

namespace std{

class Fast_Particles: public Exptime_Trajectory_List
{
  private:
    Gaussian_Comparison * gaussian_comparison;
    int displacement_time_index;
    float mindistance;
    int atoms_considered;		//number of atoms tested

  public:
    Fast_Particles();
    Fast_Particles(System * sys, Gaussian_Comparison * gc);
    void set(System *, Gaussian_Comparison *);
    void atomkernel(Trajectory * traj);
    void displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj);
    
};

}

#endif