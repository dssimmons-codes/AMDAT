// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate radial_distribution_function from n^2 method*/
/*Written by David S. Simmons*/

#include <math.h>
#include "mean_closest_distance.h"
#include "version.h"

#define PI 3.1459265

using namespace std;


Mean_Closest_Distance::Mean_Closest_Distance():Analysis_Onetime()
{
  n_times=0;
  mean_mcd=0;
  
  time_mcd= new float [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  total_center_atoms=0;
  
}



Mean_Closest_Distance::Mean_Closest_Distance(System*sys, int timescheme)
{
  int timeii;
  
  system=sys;
  
  total_center_atoms=0;
  mean_mcd=0;
  
  time_scheme = timescheme;
  n_times=determine_n_times();
  
  time_mcd= new float [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
  }
}


Mean_Closest_Distance::Mean_Closest_Distance(const Mean_Closest_Distance & copy):Analysis_Onetime(copy)
{
  int timeii;
  
  system=copy.system;
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  total_center_atoms=copy.total_center_atoms;
    mean_mcd=copy.mean_mcd;
  
  
  time_mcd= new float [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
  }
}

Mean_Closest_Distance Mean_Closest_Distance::operator=(const Mean_Closest_Distance & copy)
{
  int timeii, binii;
  
  if(this!=&copy)
  {
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  delete [] time_mcd;
    
  system=copy.system;
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  total_center_atoms=copy.total_center_atoms;
  mean_mcd=copy.mean_mcd;
  
  time_mcd= new float [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
  }
  }
  return *this;
}

void Mean_Closest_Distance::set(System*sys, int timescheme)
{
  int timeii;
  float minboxsize;
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  delete [] time_mcd;
  
  system=sys;
  
  time_scheme = timescheme;
  n_times=determine_n_times();

  
  time_mcd= new float [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  total_center_atoms=0;
    mean_mcd=0;
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
  }
}


void Mean_Closest_Distance::timekernel2(int timeii)
{
   n_atoms_i[timeii]=trajectory_list->show_n_trajectories(system_time(timeii));
   n_atoms_j[timeii]=trajectory_list2->show_n_trajectories(system_time(timeii));
   
   total_center_atoms+=n_atoms_i[timeii];
    float * time_mcd;
   trajectory_list->listloop(this,0, timeii, 0);

}


void Mean_Closest_Distance::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  tempdistance = (system->min_box_dimensions()).max();
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
  mean_mcd+=tempdistance;
}


void Mean_Closest_Distance::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    if(distance<tempdistance)
    {
        tempdistance=distance;
    }
  }
}


 
 void Mean_Closest_Distance::postprocess_list()
 {
   
   mean_mcd=mean_mcd/float(total_center_atoms);
   
 }
 
 
 void Mean_Closest_Distance::write(string filename)
 {
  int binii;

  ofstream output (filename.c_str());

  cout << "\nWriting mcd to file " <<filename<<".";cout.flush();

  output << "Mean closest distance data created by AMDAT v." << VERSION << "\n";
 
  output << mean_mcd;

  output.close();
 }
 
 
  void Mean_Closest_Distance::write(ofstream& output)const
 {
  int binii;

  cout << "\nWriting mcd to file.";

  output << "Mean closest distance data created by AMDAT v." << VERSION << "\n";
 
  output << mean_mcd;

 }

 
 
