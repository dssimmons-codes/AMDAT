/*Particles_Between class - stores a list of trajectories of particles that are geometrically between pairs of particles in a second trajectory list*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons and Pierre Kawak*/

/**/

#include <iostream>

#ifndef PARTICLES_BETWEEN
#define PARTICLES_BETWEEN

#include "trajectory_list.h"
#include "analysis_onetime.h"
#include "system.h"


namespace std{

class Particles_Between: public Trajectory_List, public Analysis_Onetime//, public Value_List <int>
{
//  float maxdistance;    //I'm guessing we need some variable like this to store the max total distance between the two secondary particles
  float dist_cutoff;         //I'm guessing we need a second variable for that lateral distance.
  float theta_cutoff;

  public:
    Particles_Between();
    Particles_Between(System*);
    Particles_Between(System*,float,float);  //put other required arguments in this one
    Particles_Between(const Particles_Between &);
    ~Particles_Between();

    Particles_Between operator = (const Particles_Between &);	//assignment

    Analysis_Type what_are_you(){Analysis_Type type = particles_between; return type;};   //this analysis type needs to be added at top of ananlysis parent class

    void preprocess();    //this is where any preprocessing prior to loops should live
    void preprocess2();    //this is where any preprocessing prior to loops should live
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void bin_hook(Trajectory_List *,int,int,int);
    void postprocess_bins();
    void set(System * syst, float d_cutoff, float t_cutoff);

//    Value_List<int> * vlist;

};

}
#endif
