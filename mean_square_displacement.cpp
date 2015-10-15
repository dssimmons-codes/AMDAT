/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_square_displacement.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "static_trajectory_list.h"
#include <omp.h>
using namespace std;


Mean_Square_Displacement::Mean_Square_Displacement()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  msd = new float [n_times];
  weighting = new float [n_times];

  atomcount = 0;
}


Mean_Square_Displacement::Mean_Square_Displacement(const Mean_Square_Displacement & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  msd = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=copy.msd[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}



/** **/
Mean_Square_Displacement::Mean_Square_Displacement(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  msd = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;

}




Mean_Square_Displacement Mean_Square_Displacement::operator = (const Mean_Square_Displacement & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] msd;
  delete [] weighting;

  msd = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=copy.msd[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Mean_Square_Displacement::initialize(System* sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] msd;
  delete [] weighting;

  msd = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=0;
    weighting[timeii]=0;
  }
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

  currenttime=thisii;
  nexttime=nextii;
  currenttimegap=timegapii;

//  weighting[timegapii]+=trajectory_list->show_n_trajectories(currenttime);
//  //weighting[timegapii]+=(trajectory_list[0]).show_n_trajectories(currenttime);
//  (trajectory_list[0]).listloop(this,currenttime);
  #pragma omp atomic
  weighting[timegapii]+=trajectory_list->show_n_trajectories(thisii);
  (trajectory_list[0]).listloop(this,timegapii, thisii, nextii);
}



void Mean_Square_Displacement::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  #pragma omp atomic
  msd[timegapii]+=pow(current_trajectory->distance(thisii,nextii),2);
}


void Mean_Square_Displacement::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        msd[timeii] /= float(weighting[timeii]);

  }
}



/*Method to write MSD data to file*/

void Mean_Square_Displacement::write(string filename)const
{
  int timeii;

  cout << "\nWriting msd to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean square displacement data created by AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<msd[timeii]<<"\n";
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
