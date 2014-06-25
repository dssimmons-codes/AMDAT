/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to allow definition of multi-trajectory bodies other than molecules
/*Written by David S. Simmons*/

#ifndef MULTIBODY
#define MULTIBODY

#include "trajectory.h"
#include coordinate.h"

namespace std
{

  class Multibody:public Trajectory
  {
    int n_trajectories;
    int * trajectories
    
    int multibody_ID;
    
    void calculate_mass;
    void calculate_center_of_mass();
    
    public:
      Multibody();
      Multibody(const Multibody &);
      ~Multibody();
      void operator =(const Multibody &);
      
      
      Coordinate shortvector(int trajii1,int trajii2,int timeii);		//Returns shortest vector between two trajectories at a given time (first trajectory index, second trajectory index, time)

      void show_coordinates(int timeii,Coordinate* list);			//returns list of coordinates of trajectories at time ii. Coordinates are returned via Coordinate* list, which must be of length equal to n_trajectories;
    
  }
  
}

#endif