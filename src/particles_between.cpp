//Add standard methods up here. Contructors, equality operator, etc.
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

#include "particles_between.h"

using namespace std;

/*Default constructor*/
Particles_Between::Particles_Between()
{
  int timeii;

  n_times = 0;

  capacity = 0;
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
  }
}

Particles_Between::Particles_Between(System * syst)
{
  int timeii;

  sys=syst;
  system = const_cast<System*>(sys);
   
  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.
  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
  }
}

Particles_Between::Particles_Between(System * syst, float maxdist, float rad)
{
  int timeii;

  sys=syst;
  system = const_cast<System*>(sys);
   
  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.
  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
  }
  maxdistance=maxdist;
  radius=rad;
}


Particles_Between::Particles_Between(const Particles_Between & copy):Trajectory_List(copy), Analysis_Onetime(copy)
{
//  Trajectory_List::Trajectory_List(copy);
//  Analysis::Analysis(copy);
  maxdistance=copy.maxdistance;
  radius=copy.radius;
  
}

Particles_Between::~Particles_Between()
{
////  Trajectory_List::~Trajectory_List();
////  Analysis::~Analysis();
}

Particles_Between Particles_Between::operator = (const Particles_Between & copy)
{
  sys=copy.system;
  system = const_cast<System*>(sys);
  trajectory_list = copy.trajectory_list;
  n_trajectories = copy.n_trajectories;
  included = copy.included;
//  Trajectory_List::operator=(copy);
//  Analysis::operator=(copy);

  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.

  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  maxdistance=copy.maxdistance;
  radius=copy.radius;

}

void Particles_Between::timekernel2(int timeii)
{
   trajectory_list->listloop(this,0, timeii, 0);  
}


void Particles_Between::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Particles_Between::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii) //need to check back on which of these time indices is the right one.
{
   Trajectory* traj3;   
   bool isbetween;    //boolean to store whether the particle is in the 'between' set at this time
   
   int n_trajs=trajectory_list2->show_n_trajectories(thisii);  //record how many particles are in the second trajectory list at this time
   
   for(int trajii=0;trajii<n_trajs;trajii++)  //here is the second loop over particles
   {
      traj3 = (*trajectory_list2)(thisii,trajii);  //this is the syntax to pull the second particle in set 2
      //we now have pointers to all 3 particles here.
      //traj1 is the central particle you are trying to sort
      //traj2 is the first particle in the group you are looking between
      //traj3 is the second particle in the group you are looking between
         
      //do some math here to determine whether particle is in between
      //if it is, set isbetween to true.
      //if it is not, set isbetween to false
      
         if(isbetween)
         {
            addtrajectory(thisii,traj1);        //this line will add the trajectory to the trajectory list
         }
   }
      
  float distance;
  if(traj2!=traj3)
  {
    distance=(traj3->show_coordinate(thisii)-(traj2->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
//    bin(thisii,distance);
  }
   
}

void Particles_Between::preprocess()
{
}

void Particles_Between::preprocess2()
{
}

void Particles_Between::analyze (Trajectory_List *)
{
}

void Particles_Between::analyze (Trajectory_List *, Trajectory_List*)
{
}

void Particles_Between::bin_hook(Trajectory_List *,int,int,int)
{
}

void Particles_Between::postprocess_bins()
{
  //I don't think any postprocessing is needed in this case
}

