// Copyright (C) 2010-2025 David S. Simmons

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
}


Fast_Particles::Fast_Particles(System * sys, int timeindex, float distance_threshold)
{
  int timeii;

  system = sys;
  displacement_time_index = timeindex;
  mindistance = distance_threshold;
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
}


void Fast_Particles::set(System * syst, Gaussian_Comparison * gc)
{

  int cap = syst->show_n_trajectories();
  reset(syst,cap);

  gaussian_comparison = gc;
  displacement_time_index = gaussian_comparison->show_time_index();
  mindistance = gaussian_comparison->show_fastboundary();

}



void Fast_Particles::set(System * sys, int timeindex, float distance_threshold)
{

  int cap = sys->show_n_trajectories();
  reset(sys,cap);

  gaussian_comparison = 0;
  displacement_time_index = timeindex;
  mindistance = distance_threshold;

}


void Fast_Particles::analyze(Trajectory_List * t_list)
{
  trajectory_list=t_list;
  system->displacement_list(this,displacement_time_index,bool(0));
  postprocess_list();
}


void Fast_Particles::list_displacementkernel(int timegapii,int thisii, int nextii)
{

//  weighting[timegapii]+=trajectory_list->show_n_trajectories(currenttime);
//  //weighting[timegapii]+=(trajectory_list[0]).show_n_trajectories(currenttime);
//  (trajectory_list[0]).listloop(this,currenttime);
  (trajectory_list[0]).listloop(this,timegapii, thisii, nextii);
}


void Fast_Particles::listkernel(Trajectory* current_trajectory, int timegapii, int currenttime, int nexttime)
{
  float square_displacement;		//particle displacement
  int expindex;			//index of exponential block

  square_displacement = pow(current_trajectory->distance(currenttime,nexttime),2);  //calculate particle displacement

  expindex = int((float(currenttime)-float(system->show_frt()))/float(system->show_n_exponential_steps()));	//calculate which exponential block this corresponds to
  if(square_displacement > mindistance)
  {
    if(n_trajectories[expindex]==capacity){cout<<"Error: particle list memory allocation full.\n";exit(1);}
    addtrajectory(expindex,current_trajectory);
  }
}