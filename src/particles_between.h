/*Particles_Between class - stores a list of trajectories of particles that are geometrically between pairs of particles in a second trajectory list*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons and Pierre Kawak*/

/**/

#ifndef PARTICLES_BETWEEN
#define PARTICLES_BETWEEN

#include <iostream>
#include "trajectory_list.h"
#include "analysis_onetime.h"
#include "system.h"

namespace std{

class Particles_Between: public Trajectory_List, public Analysis_Onetime
{
  float maxdistance;    //I'm guessing we need some variable like this to store the max total distance between the two secondary particles
  float radius;         //I'm guessing we need a second variable for that lateral distance. 
  
  public:
  Particles_Between();
  Particles_Between(System*);
  Particles_Between(System*,float,float);  //put other required arguments in this one
  Particles_Between(const Particles_Between &);
  ~Particles_Between();
  
  Particles_Between operator = (const Particles_Between &);	//assignment
  
  Analysis_Type what_are_you(){Analysis_Type type = particles_between; return type;};   //this analysis type needs to be added at top of ananlysis parent class
  
  void preprocess();    //this is where any preprocessing prior to loops should live 
  void timekernel(int timeii){timekernel2(timeii);};
  void timekernel2(int timeii);
  void listkernel(Trajectory *, int, int, int);
  void listkernel2(Trajectory *, Trajectory *, int, int, int);
  void postprocess_list();
  
};

}
#endif
