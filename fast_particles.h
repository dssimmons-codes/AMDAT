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

  public:
    Fast_Particles();
    Fast_Particles(System * sys, Gaussian_Comparison * gc);
    Fast_Particles(System * sys, int timeindex, float distancethreshold);	//construct to find fast particles based on user-specificed distance cutoff
    void set(System *, Gaussian_Comparison *);		//construct to find fast_particles based on distance cutoff specified by crossover between actual and Gaussian Van Hove computed by Gaussian_Comparison
    
    Analysis_Type what_are_you(){Analysis_Type type = fast_particles; return type;};		//virtual method to report the type of analysis
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *, int, int, int);
};

}

#endif