/*Fast_Particles class: identifies and stores list of fast particles at start of each exponential block*/
/*Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#ifndef AVERAGE_PARTICLES
#define AVERAGE_PARTICLES

#include "analysis.h"
#include "exptime_trajectory_list.h"
#include "system.h"
#include "gaussian_comparison.h"

namespace std{

class Average_Particles: public Exptime_Trajectory_List
{
    Gaussian_Comparison * gaussian_comparison;
    int displacement_time_index;
    float mindistance;
    float maxdistance;
    int atoms_considered;		//number of atoms tested

  public:
    Average_Particles();
    Average_Particles(System * sys, Gaussian_Comparison * gc);
    void set(System *, Gaussian_Comparison *);
    void atomkernel(Trajectory * traj);
    void displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj);
};

}

#endif