/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_velocity_unsteady.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "static_trajectory_list.h"
#include <omp.h>
using namespace std;


Mean_Velocity_Unsteady::Mean_Velocity_Unsteady()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  mean_velocity = new Coordinate [n_times];
  weighting = new int [n_times];
}


Mean_Velocity_Unsteady::Mean_Velocity_Unsteady(const Mean_Velocity_Unsteady & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;

  mean_velocity = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    mean_velocity[timeii]=copy.mean_velocity[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}



Mean_Velocity_Unsteady::Mean_Velocity_Unsteady(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timesteps();

   //allocate memory for mean square displacement data
  mean_velocity = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
  }

}




Mean_Velocity_Unsteady Mean_Velocity_Unsteady::operator = (const Mean_Velocity_Unsteady & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;

  delete [] mean_velocity;
  delete [] weighting;

  mean_velocity = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    mean_velocity[timeii]=copy.mean_velocity[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Mean_Velocity_Unsteady::initialize(System* sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timesteps();

   //allocate memory for mean square displacement data

  delete [] mean_velocity;
  delete [] weighting;

  mean_velocity = new Coordinate [n_times];
  weighting = new int [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
  }
}



/*Methods to do analysis using trajectory list*/

void Mean_Velocity_Unsteady::analyze(Trajectory_List * t_list)
{
  int timeii;
  trajectory_list=t_list;
  for(timeii=0;timeii<n_times;timeii++)
  {
    (trajectory_list[0]).listloop(this,0, timeii, 0);
  }
  
  //system->displacement_list(this);
  postprocess_list();
}




void Mean_Velocity_Unsteady::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  #pragma omp atomic
  weighting[thisii]++;
  mean_velocity[thisii]+=current_trajectory->show_velocity(thisii);
}


void Mean_Velocity_Unsteady::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        mean_velocity[timeii] /= float(weighting[timeii]);

  }
}



/*Method to mean velocity data to file*/

void Mean_Velocity_Unsteady::write(string filename)
{
  int timeii;

  cout << "\nWriting mean velocity to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean velocity data created by AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << system->show_time(timeii)<<"\t"<<mean_velocity[timeii].show_x()<<"\t"<<mean_velocity[timeii].show_y()<<"\t"<<mean_velocity[timeii].show_z()<<"\n";
  }
}