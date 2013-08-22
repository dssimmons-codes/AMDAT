/*Methods for Fast_Particles class: identifies and stores list of fast particles at start of each exponential block*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#include "fast_particles.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;

Fast_Particles::Fast_Particles()
{
  int timeii;

  displacement_time_index=-1;
  mindistance = -1;
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

  atoms_considered = 0;
}



Fast_Particles::Fast_Particles(System * sys, Gaussian_Comparison * gc)
{
  int timeii;

  system = sys;
  gaussian_comparison = gc;
  displacement_time_index = gaussian_comparison->show_time_index();
  mindistance = gaussian_comparison->show_fastboundary();
  capacity=system->show_n_atoms()+system->show_n_molecules();
  n_times = system->show_n_exponentials();

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
  }
  atoms_considered = 0;
}



void Fast_Particles::set(System * syst, Gaussian_Comparison * gc)
{

  int cap = syst->show_n_trajectories();
  reset(syst,cap);

  gaussian_comparison = gc;
  displacement_time_index = gaussian_comparison->show_time_index();
  mindistance = gaussian_comparison->show_fastboundary();



  atoms_considered = 0;
}


void Fast_Particles::atomkernel(Trajectory * traj)
{
  system->displacement_loop(this, traj, displacement_time_index,bool(0));
  atoms_considered++;
}



void Fast_Particles::displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj)
{
  float square_displacement;		//particle displacement
  int expindex;			//index of exponential block

  square_displacement = pow(traj->distance(thisii,nextii),2);  //calculate particle displacement

  expindex = int((float(thisii)-float(system->show_frt()))/float(system->show_n_exponential_steps()));	//calculate which exponential block this corresponds to
  if(square_displacement > mindistance)
  {
    if(n_trajectories[expindex]==capacity){cout<<"Error: particle list memory allocation full.\n";exit(1);}
    addtrajectory(expindex,traj);
  }
}
