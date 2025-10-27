// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_unsteady_displacement.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "static_trajectory_list.h"
#include <omp.h>
using namespace std;


Mean_Unsteady_Displacement::Mean_Unsteady_Displacement()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  md = new Coordinate [n_times];
  weighting = new int [n_times];
}


Mean_Unsteady_Displacement::Mean_Unsteady_Displacement(const Mean_Unsteady_Displacement & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii]=copy.md[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}



Mean_Unsteady_Displacement::Mean_Unsteady_Displacement(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timesteps();

   //allocate memory for mean square displacement data
  md = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
  }

}




Mean_Unsteady_Displacement Mean_Unsteady_Displacement::operator = (const Mean_Unsteady_Displacement & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;

  delete [] md;
  delete [] weighting;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii]=copy.md[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Mean_Unsteady_Displacement::initialize(System* sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timesteps();

   //allocate memory for mean square displacement data

  delete [] md;
  delete [] weighting;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
  }
}



/*Methods to do analysis using trajectory list*/

void Mean_Unsteady_Displacement::analyze(Trajectory_List * t_list)
{
  int timeii;
  trajectory_list=t_list;
  
  weighting[0] = 1;
  
  for(timeii=0;timeii<n_times-1;timeii++)
  {
    (trajectory_list[0]).listloop(this,0, timeii, 0);
  }
  
  //system->displacement_list(this);
  postprocess_list();
}



void Mean_Unsteady_Displacement::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  #pragma omp atomic
  weighting[thisii+1]++;
  md[thisii+1]+=current_trajectory->show_unwrapped(thisii+1)-current_trajectory->show_unwrapped(thisii);
}


void Mean_Unsteady_Displacement::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        md[timeii] /= float(weighting[timeii]);

  }
}



/*Method to mean velocity data to file*/

void Mean_Unsteady_Displacement::write(string filename)
{
  int timeii;

  cout << "\nWriting mean incremental displacement to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean incremental displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << system->show_time(timeii)<<"\t"<<md[timeii].show_x()<<"\t"<<md[timeii].show_y()<<"\t"<<md[timeii].show_z()<<"\n";
  }
}


void Mean_Unsteady_Displacement::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting mean incremental displacement to file.";

  output << "Mean incremental displacement data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << system->show_time(timeii)<<"\t"<<md[timeii].show_x()<<"\t"<<md[timeii].show_y()<<"\t"<<md[timeii].show_z()<<"\n";
  }
}