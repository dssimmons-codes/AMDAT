// Copyright (C) 2010-2025 David S. Simmons

/*Amorphou Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_square_displacement.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "static_trajectory_list.h"
#include <omp.h>
using namespace std;

Mean_Square_Displacement::Mean_Square_Displacement()
{
  n_times = 0;
  atomcount = 0;
}

Mean_Square_Displacement::~Mean_Square_Displacement() = default;

Mean_Square_Displacement::Mean_Square_Displacement(const Mean_Square_Displacement & copy)
{
  system = copy.system;
  trajectory_list = copy.trajectory_list;
  n_times = copy.n_times;
  atomcount = copy.atomcount;
  msd = copy.msd;
  weighting = copy.weighting;
  timetable = copy.timetable;
}



/** **/
Mean_Square_Displacement::Mean_Square_Displacement(System*sys):Mean_Square_Displacement()
{
  initialize(sys);
}




Mean_Square_Displacement& Mean_Square_Displacement::operator = (const Mean_Square_Displacement & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    trajectory_list = copy.trajectory_list;
    n_times = copy.n_times;
    atomcount = copy.atomcount;
    msd = copy.msd;
    weighting = copy.weighting;
    timetable = copy.timetable;
  }
  return *this;
}


void Mean_Square_Displacement::initialize(System* sys)
{
  system = sys;
  n_times = system->show_n_timegaps();
  msd.assign(n_times, {});
  weighting.assign(n_times, {});

  float * times = system->displacement_times();
  timetable.assign(times, times+n_times);
  delete [] times;

  atomcount = 0;
}



/*Methods to do analysis using trajectory list*/

void Mean_Square_Displacement::analyze(Trajectory_List * t_list)
{
  trajectory_list=t_list;
  system->displacement_list(this,false);
  postprocess_list();
}

void Mean_Square_Displacement::list_displacementkernel(int timegapii,int thisii, int nextii)
{

  int currenttime=thisii;
  int nexttime=nextii;
  int currenttimegap=timegapii;

//  weighting[timegapii]+=trajectory_list->show_n_trajectories(currenttime);
//  //weighting[timegapii]+=(trajectory_list[0]).show_n_trajectories(currenttime);
//  (trajectory_list[0]).listloop(this,currenttime);
  weighting[timegapii][0]+=trajectory_list->show_n_trajectories(thisii);
  (trajectory_list[0]).listloop(this,timegapii, thisii, nextii);
}



void Mean_Square_Displacement::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  msd[timegapii][0]+=current_trajectory->distance(thisii,nextii)*current_trajectory->distance(thisii,nextii);
}


void Mean_Square_Displacement::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        msd[timeii][0] /= float(weighting[timeii][0]);

  }
}



/*Method to write MSD data to file*/

void Mean_Square_Displacement::write(string filename)const
{
  int timeii;

  cout << "\nWriting msd to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean square displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<msd[timeii][0]<<"\n";
  }
}


void Mean_Square_Displacement::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting msd to file.";

  output << "Mean square displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<msd[timeii][0]<<"\n";
  }
}

void Mean_Square_Displacement::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Mean_Square_Displacement::postprocess_bins()
{
  postprocess_list();
}
