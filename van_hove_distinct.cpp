/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Van_Hove_Distinct: a class for distinct part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#include "van_hove_distinct.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;


/*----------------------------------------------------------------------------------------------*/

/*Constructor*/

Van_Hove_Distinct::Van_Hove_Distinct()
{
  system = 0;
  n_bins = 0;
  max_value = 0;
  n_times = 0;
  bin_size = 0;
  timetable = 0;
  correlation = new float * [1];
  weighting = new int [0];
  correlation[0]=new float [0];
}



/*----------------------------------------------------------------------------------------------*/

Van_Hove_Distinct::Van_Hove_Distinct(System*sys, Trajectory_List_Bins binnedtraj, int bin_count, float value_max)
{
  int timeii, binii;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = true;
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
  
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
  
}


Van_Hove_Distinct::Van_Hove_Distinct(System*sys, int bin_count, float value_max)
{
  int timeii, binii;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = false;
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
 
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
  
}



void Van_Hove_Distinct::set(System*sys, int bin_count, float value_max)
{
  int timeii, binii;
  
  delete [] trajectory_list;
//  delete [] currentlists;
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  
  delete [] correlation;
  delete [] weighting;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = false;
  
  
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
  
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
 
  
}


void Van_Hove_Distinct::analyze(Trajectory_List * t_list1)
{

  
  trajectory_list[0]=t_list1[0];
  trajectory_list[1]=t_list1[0];
 
  system->displacement_list(this);
  postprocess_list();
}


void Van_Hove_Distinct::analyze(Trajectory_List * t_list1, Trajectory_List * t_list2)
{
  trajectory_list[0]=t_list1[0];
  trajectory_list[1]=t_list2[0];
  
  if(!use_binned)
  {
    currentlist0=&trajectory_list[0];
    currentlist1=&trajectory_list[1];
  }
  
  system->displacement_list(this);
  postprocess_list();
  
}




void Van_Hove_Distinct::list_displacementkernel(int timegapii, int thisii, int nextii)
{
  int xii, yii, zii;	//bin of first trajectory_list

  if(use_binned)
  {
    for(xii=0;xii<nx;xii++)
    {
      for(yii=0;yii<ny;yii++)
      {
	for(zii=0;zii<nz;zii++)
	{
	
	
	}
      }
    }
  }
  else
  {
    weighting[timegapii]+=trajectory_list[0].show_n_trajectories(thisii)*trajectory_list[1].show_n_trajectories(nextii);
//    current_list1->listloop(this,timegapii, thisii, nextii);
  }
}


void Van_Hove_Distinct::listkernel(Trajectory* traj1, int timegapii, int thisii, int nextii)
{
//  current_list2->listloop(this, traj1, timegapii, thisii, nextii);
}




void Van_Hove_Distinct::listkernel2(Trajectory* traj1, Trajectory* traj2, int timegapii, int thisii, int nextii)
{
  float distance;
  if(traj1!=traj2)
  {
    (traj2->show_coordinate(nextii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    bin(timegapii,distance);
  }
  else
  {
    weighting[timegapii]--;
  }
  
  
}

