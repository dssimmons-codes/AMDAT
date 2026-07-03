// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_displacement.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "static_trajectory_list.h"
using namespace std;

const int PAD = 16;

Mean_Displacement::Mean_Displacement()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  md = new Coordinate * [n_times];
  weighting = new int * [n_times];

  atomcount = 0;
}


Mean_Displacement::Mean_Displacement(const Mean_Displacement & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  md = new Coordinate * [n_times];
  weighting = new int * [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii] = new Coordinate[PAD];
    weighting[timeii] = new int [PAD];
  }

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii][0]=copy.md[timeii][0];
    weighting[timeii][0]=copy.weighting[timeii][0];
  }
}



/** **/
Mean_Displacement::Mean_Displacement(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  md = new Coordinate * [n_times];
  weighting = new int * [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii] = new Coordinate[PAD];
    weighting[timeii] = new int [PAD];
  }

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii][0]=0;
  }
  atomcount = 0;

}




Mean_Displacement Mean_Displacement::operator = (const Mean_Displacement & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] md[timeii];
    delete [] weighting[timeii];
  }

  delete [] md;
  delete [] weighting;

  md = new Coordinate * [n_times];
  weighting = new int * [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii] = new Coordinate[PAD];
    weighting[timeii] = new int [PAD];
  }

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii][0]=copy.md[timeii][0];
    weighting[timeii][0]=copy.weighting[timeii][0];
  }

  }

  return *this;

}


void Mean_Displacement::initialize(System* sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] md[timeii];
    delete [] weighting[timeii];
  }

  delete [] md;
  delete [] weighting;

  md = new Coordinate * [n_times];
  weighting = new int * [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii] = new Coordinate[PAD];
    weighting[timeii] = new int [PAD];
  }

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii][0]=0;
  }
  atomcount = 0;
}

void Mean_Displacement::preprocess()
{
  int timeii;
  weighting_temp = system->timegap_weighting();
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii][0]=weighting_temp[timeii];
  }
}


/*Methods to do analysis using trajectory list*/

void Mean_Displacement::analyze(Trajectory_List * t_list)
{
  trajectory_list=t_list;

  system->displacement_list(this);
  postprocess_list();
}

void Mean_Displacement::list_displacementkernel(int timegapii,int thisii, int nextii)
{

  int currenttimegap=timegapii;
  int currenttime=thisii;
  int nexttime=nextii;

  weighting[timegapii][0]+=trajectory_list->show_n_trajectories(currenttime);
  //weighting[timegapii]+=(trajectory_list[0]).show_n_trajectories(currenttime);
  (trajectory_list[0]).listloop(this,currenttimegap,currenttime,nexttime);
}



void Mean_Displacement::listkernel(Trajectory* current_trajectory, int currenttimegap, int currenttime, int nexttime)
{
  Coordinate c_next = current_trajectory->show_unwrapped(nexttime);
  Coordinate c_this = current_trajectory->show_unwrapped(currenttime);
  md[currenttimegap][0]+=c_next-c_this;
}



void Mean_Displacement::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        md[timeii][0] /= float(weighting[timeii][0]);

  }
}



/*Method to write MD data to file*/

void Mean_Displacement::write(string filename)
{
  int timeii;

  cout << "\nWriting mean displacement to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<md[timeii][0].show_x()<<"\t"<<md[timeii][0].show_y()<<"\t"<<md[timeii][0].show_z()<<"\n";
  }
}

void Mean_Displacement::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting mean displacement to file.";

  output << "Mean displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<md[timeii][0].show_x()<<"\t"<<md[timeii][0].show_y()<<"\t"<<md[timeii][0].show_z()<<"\n";
  }
}



void Mean_Displacement::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Mean_Displacement::postprocess_bins()
{
  postprocess_list();
}
