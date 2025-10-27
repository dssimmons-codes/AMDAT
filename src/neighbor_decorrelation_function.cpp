// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate neighbor decorrelation function*/
/*Written by David S. Simmons*/

#include "neighbor_decorrelation_function.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "static_trajectory_list.h"
#include <omp.h>
using namespace std;


Neighbor_Decorrelation_Function::Neighbor_Decorrelation_Function()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  ndf = new float [n_times];
  weighting = new float [n_times];

  atomcount = 0;
  
  n_list=0;
}


Neighbor_Decorrelation_Function::Neighbor_Decorrelation_Function(const Neighbor_Decorrelation_Function & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  ndf = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    ndf[timeii]=copy.ndf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  
  n_list=copy.n_list;
}



/** **/
Neighbor_Decorrelation_Function::Neighbor_Decorrelation_Function(System*sys, Neighbor_List* nlist)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  ndf = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    ndf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  
  n_list=nlist;

}




Neighbor_Decorrelation_Function Neighbor_Decorrelation_Function::operator = (const Neighbor_Decorrelation_Function & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] ndf;
  delete [] weighting;

  ndf = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    ndf[timeii]=copy.ndf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  
  n_list=copy.n_list;

  }

  return *this;

}


void Neighbor_Decorrelation_Function::initialize(System* sys, Neighbor_List* nlist)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] ndf;
  delete [] weighting;

  ndf = new float [n_times];
  weighting = new float [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    ndf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  
  n_list=nlist;
}






void Neighbor_Decorrelation_Function::analyze(Trajectory_List * t_list)
{
  trajectory_list=t_list;
  system->displacement_list(this,false);
  postprocess_list();
}

void Neighbor_Decorrelation_Function::list_displacementkernel(int timegapii,int thisii, int nextii)
{
  weighting[timegapii]+=trajectory_list->show_n_trajectories(thisii);
  (trajectory_list[0]).listloop(this,timegapii, thisii, nextii);
}



void Neighbor_Decorrelation_Function::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  int trajID;
  
  trajID=current_trajectory->show_trajectory_ID();
  
  ndf[timegapii]+=n_list->n_persistent_neighbors(trajID,thisii,nextii);
}


void Neighbor_Decorrelation_Function::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        ndf[timeii] /= float(weighting[timeii]);

  }
}



/*Method to write MSD data to file*/

void Neighbor_Decorrelation_Function::write(string filename)const
{
  int timeii;

  cout << "\nWriting neighbor decorrelation function to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Neighbor decorrelation function data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ndf[timeii]<<"\n";
  }
}


void Neighbor_Decorrelation_Function::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting neighbor decorrelation function to file.";

  output << "Neighbor decorrelation function data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ndf[timeii]<<"\n";
  }
}

void Neighbor_Decorrelation_Function::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Neighbor_Decorrelation_Function::postprocess_bins()
{
  postprocess_list();
}
