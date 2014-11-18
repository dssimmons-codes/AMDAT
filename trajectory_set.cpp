/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Trajectory_set class - Class to store sets of user-defined trajectories*/
/*Written by David S. Simmons*/

#include "trajectory_set.h"
#include "trajectory.h"

using namespace std;

Trajectory_Set::Trajectory_Set()
{
    n_trajectories=0;
    trajectories=new Trajectory [n_trajectories];
}

Trajectory_Set::Trajectory_Set(const Trajectory_Set & copy)
{
    int trajii;

    n_trajectories = copy.n_trajectories;
    trajectories=new Trajectory [n_trajectories];

    for(trajii=0;trajii<n_trajectories;trajii++)
    {
        trajectories[trajii]=copy.trajectories[trajii];
    }
}

Trajectory_Set::~Trajectory_Set()
{
    delete [] trajectories;
}


Trajectory_Set Trajectory_Set::operator=(const Trajectory_Set & copy)
{
    int trajii;

    if(this!=&copy)
    {
        delete [] trajectories;

        n_trajectories = copy.n_trajectories;
        trajectories=new Trajectory [n_trajectories];

        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=copy.trajectories[trajii];
        }
    }
    
    return *this;
}



Trajectory_Set::Trajectory_Set(Multibody_Set * multibodies, bool centertype)
{
    int trajii;
    n_trajectories = multibodies->show_n_multibodies();

    trajectories = new Trajectory [n_trajectories];

    if(centertype) //compute and save centroid trajectory
    {
        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=multibodies->show_multibody(trajii)->centroid_trajectory();
        }
    }
    else    //compute and save COM trajectory
    {
        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=multibodies->show_multibody(trajii)->center_of_mass_trajectory();
        }
    }
}



void Trajectory_Set::trajectories_from_multibodies(Multibody_Set * multibodies, bool centertype)
{
  if(centertype)
  {
    trajectories_from_centroid(multibodies);
  }
  else
  {
    trajectories_from_com(multibodies);
  }
}


void Trajectory_Set::trajectories_from_centroid(Multibody_Set * multibodies)
{
  
  int trajii;
 
  delete [] trajectories;
  
  n_trajectories = multibodies->show_n_multibodies();

  trajectories = new Trajectory [n_trajectories];
  
  for(trajii=0;trajii<n_trajectories;trajii++)
  {
    trajectories[trajii]=multibodies->show_multibody(trajii)->centroid_trajectory();
  }
}




void Trajectory_Set::trajectories_from_com(Multibody_Set * multibodies)
{
  
  int trajii;
 
  delete [] trajectories;
  
  n_trajectories = multibodies->show_n_multibodies();

  trajectories = new Trajectory [n_trajectories];
  
  for(trajii=0;trajii<n_trajectories;trajii++)
  {
    trajectories[trajii]=multibodies->show_multibody(trajii)->centroid_trajectory();
  }
}
