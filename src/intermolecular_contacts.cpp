/*Methods for Intermolecular_Contacts class: identifies and stores list of particles that are contacting (within specified distance threshold) a particle in another molecule*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include "intermolecular_contacts.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;

Intermolecular_Contacts::Intermolecular_Contacts()
{
  int timeii;

  threshold = -1;
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



Intermolecular_Contacts::Intermolecular_Contacts(System * sys)
{
  int timeii;

  system = sys;

  threshold=-1;
  capacity=system->show_n_atoms()+system->show_n_molecules();
  n_times = system->show_n_timesteps();

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


Intermolecular_Contacts::Intermolecular_Contacts(System * sys, float distance_threshold)
{
  int timeii;

  system = sys;
  threshold = distance_threshold;
  capacity=system->show_n_atoms()+system->show_n_molecules();
  n_times = system->show_n_timesteps();

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


void Intermolecular_Contacts::set(System * syst)
{

  int cap = syst->show_n_trajectories();
  reset(syst,cap);

  threshold=-1;

}



void Intermolecular_Contacts::set(System * sys, float distance_threshold)
{

  int cap = sys->show_n_trajectories();
  reset(sys,cap);

  threshold = distance_threshold;

}


//copy constructor
Intermolecular_Contacts::Intermolecular_Contacts(const Intermolecular_Contacts & copy):Analysis_Onetime(copy),Everytime_Trajectory_List(copy)
{
    threshold=copy.threshold;
}


//assignment operator
Intermolecular_Contacts Intermolecular_Contacts::operator = (const Intermolecular_Contacts & copy)
{
    threshold=copy.threshold;
    system=copy.system;
    int cap = system->show_n_trajectories();
    reset(system, cap);
    Analysis_Onetime::operator=(copy);
    Everytime_Trajectory_List::operator=(copy);
}


void Intermolecular_Contacts::timekernel2(int timeii)
{
   trajectory_list->listloop(this,0, timeii, 0);
}


void Intermolecular_Contacts::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Intermolecular_Contacts::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;

  if (is_included(thisii, traj1->show_trajectory_ID())) //skip if traj1 is already included
    return;

  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    if (distance<threshold)
    {
        if(n_trajectories[timegapii]==capacity){cout<<"Error: particle list memory allocation full.\n";exit(1);}
        addtrajectory(thisii,traj1);

    }
  }
}
