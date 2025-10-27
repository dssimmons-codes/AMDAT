// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate radial_distribution_function from n^2 method*/
/*Written by David S. Simmons*/

#include <math.h>
#include "radial_distribution_function.h"
#include "generated/version.h"

#define PI 3.1459265

using namespace std;


Radial_Distribution_Function::Radial_Distribution_Function():Analysis_Onetime()
{
  max_distance=0;
  n_bins=0;
  bin_size=0;
  n_times=0;
  
  time_rdf= new float * [n_times];
  mean_rdf= new float [n_bins];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
}



Radial_Distribution_Function::Radial_Distribution_Function(System*sys, int nbins, int timescheme, float maxdistance)
{
  int timeii, binii;
  float minboxsize;
  
  system=sys;
  n_bins=nbins;
  minboxsize=(system->min_box_dimensions()).min();

  if(std::fabs(minboxsize - 1.0f) <= 1e-7f)
  {
    minboxsize=(system->min_box_dimensions()).min2D();
  }

  if(maxdistance>0&&maxdistance<minboxsize/2)
  {
    max_distance=maxdistance;
  }
  else
  {
    max_distance = minboxsize/2;
    cout << "\nMax length scale binned selected by user for rdf is not between 0 and half the smallest box dimension. Defaulting to half smallest box dimension. Max_distance: " << max_distance << "\n";
  }
  
  time_scheme = timescheme;
  n_times=determine_n_times();
  bin_size = max_distance/float(n_bins);
  
  time_rdf= new float * [n_times];
  mean_rdf= new float [n_bins];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_rdf[timeii]=new float [n_bins];
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
    for(binii=0;binii<n_bins;binii++)
    {
      time_rdf[timeii][binii]=0;
    }
    for(binii=0;binii<n_bins;binii++)
    {
      mean_rdf[binii]=0;
    }
  }
}


Radial_Distribution_Function::Radial_Distribution_Function(const Radial_Distribution_Function & copy):Analysis_Onetime(copy)
{
  int timeii, binii;
  
  system=copy.system;
  n_bins=copy.n_bins;
  max_distance=copy.max_distance;
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  bin_size = copy.bin_size;
  
  time_rdf= new float * [n_times];
  mean_rdf= new float [n_bins];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_rdf[timeii]=new float [n_bins];
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
    for(binii=0;binii<n_bins;binii++)
    {
      time_rdf[timeii][binii]=copy.time_rdf[timeii][binii];
    }
    for(binii=0;binii<n_bins;binii++)
    {
      mean_rdf[binii]=copy.mean_rdf[binii];
    }
  }
}

Radial_Distribution_Function Radial_Distribution_Function::operator=(const Radial_Distribution_Function & copy)
{
  int timeii, binii;
  
  if(this!=&copy)
  {
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  delete [] mean_rdf;
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] time_rdf[timeii];
  }
  delete [] time_rdf;
    
  system=copy.system;
  n_bins=copy.n_bins;
  max_distance=copy.max_distance;
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  bin_size = copy.bin_size;
  
  time_rdf= new float * [n_times];
  mean_rdf= new float [n_bins];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_rdf[timeii]=new float [n_bins];
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
    for(binii=0;binii<n_bins;binii++)
    {
      time_rdf[timeii][binii]=copy.time_rdf[timeii][binii];
    }
    for(binii=0;binii<n_bins;binii++)
    {
      mean_rdf[binii]=copy.mean_rdf[binii];
    }
  }
  }
  return *this;
}

void Radial_Distribution_Function::set(System*sys, int nbins, int timescheme, float maxdistance)
{
  int timeii, binii;
  float minboxsize;
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  delete [] mean_rdf;
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] time_rdf[timeii];
  }
  delete [] time_rdf;
  
  system=sys;
  n_bins=nbins;
  minboxsize=(system->min_box_dimensions()).min();
  
  if(std::fabs(minboxsize - 1.0f) <= 1e-7f)
  {
    minboxsize=(system->min_box_dimensions()).min2D();
  }

  if(maxdistance>0&&maxdistance<minboxsize/2)
  {
    max_distance=maxdistance;
  }
  else
  {
    max_distance = minboxsize/2;
    cout << "\nMax length scale binned selected by user for rdf is not between 0 and half the smallest box dimension. Defaulting to half smallest box dimension.\n";
  }
  
  time_scheme = timescheme;
  n_times=determine_n_times();
  bin_size = max_distance/float(n_bins);
  
  time_rdf= new float * [n_times];
  mean_rdf= new float [n_bins];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_rdf[timeii]=new float [n_bins];
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
    for(binii=0;binii<n_bins;binii++)
    {
      time_rdf[timeii][binii]=0;
    }
  }
  for(binii=0;binii<n_bins;binii++)
    {
      mean_rdf[binii]=0;
    }
}


void Radial_Distribution_Function::timekernel2(int timeii)
{
   n_atoms_i[timeii]=trajectory_list->show_n_trajectories(system_time(timeii));
   n_atoms_j[timeii]=trajectory_list2->show_n_trajectories(system_time(timeii));
   trajectory_list->listloop(this,0, timeii, 0);
}


void Radial_Distribution_Function::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Radial_Distribution_Function::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    bin(thisii,distance);
  }
}


void Radial_Distribution_Function::bin(int timestep, float distance)
 {
  int binindex;
  binindex = int((distance)/bin_size);
  
  if(binindex>=0 && binindex<n_bins)
  {
    (time_rdf[timestep][binindex])++;
  }
 }
 
 
 void Radial_Distribution_Function::postprocess_list()
 {
   int timeii,binii;
   float boxvolume;
   Coordinate boxsize;
   float rhoj;
   float rshell,shellvolume;
   int n_i_total=0;
   
   for(timeii=0;timeii<n_times;timeii++)
   {
     boxsize = system->size(system_time(timeii));
     boxvolume=boxsize.show_x()*boxsize.show_y()*boxsize.show_z();
     rhoj=float(n_atoms_j[timeii])/boxvolume;
     n_i_total+=n_atoms_i[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       rshell = float(binii)*bin_size;						//determine inner radius of bin
       shellvolume = (4.0/3.0)*PI*(pow(rshell+bin_size,3.0)-pow(rshell,3.0));		//calculate volume of bin
       time_rdf[timeii][binii]/=(float(n_atoms_i[timeii])*rhoj*shellvolume);
       mean_rdf[binii]+=time_rdf[timeii][binii]*float(n_atoms_i[timeii]);
     }
   }
   
   for(binii=0;binii<n_bins;binii++)
   {
     
     mean_rdf[binii]/=float(n_i_total);
   }
 }
 
 #ifdef NEVER
 
 Radial_Distribution_Function operator+ (const Radial_Distribution_Function & increment)
 {
   
   Radial_Distribution_Function() temp;
   
   
   //compute time_rdf
   for(timeii=0;timeii<n_times;timeii++)
   {
     temp.n_atoms_i[timeii]=n_atoms_i+increment[timeii].i_atoms_j[timeii];
     temp.n_atoms_j[timeii]=n_atoms_j+increment[timeii].n_atoms_j[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       time_rdf
     }
   }
 
   
   //compute and normalize mean_rdf
   int n_i_total = 0;
   for(timeii=0;timeii<n_times;timeii++)
   {
     n_i_total+=n_atoms_i[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       temp.mean_rdf[binii]+=time_rdf[timeii]*float(n_atoms_i[timeii]);
     }
   }
   
   for(binii=0;binii<n_bins;binii++)
   {
     temp/ mean_rdf[binii]/=float(n_i_total);
   }
 }
 
 #endif
 
 void Radial_Distribution_Function::write(string filename)
 {
  int binii;

  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Radial distribution function data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(binii=0;binii<n_bins;binii++)
  {
    output << float(binii)*bin_size+bin_size/2 << "\t";
  }
  output << "\n";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << mean_rdf[binii]  <<  "\t";
  }
  
  output << "\n";

  output.close();
 }
 
 
  void Radial_Distribution_Function::write(ofstream& output)const
 {
  int binii;

  cout << "\nWriting rdf to file.";

  /*Write first row - list of bin numbers*/
  output << "Radial distribution function data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(binii=0;binii<n_bins;binii++)
  {
    output << float(binii)*bin_size+bin_size/2 << "\t";
  }
  output << "\n";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << mean_rdf[binii]  <<  "\t";
  }
  
  output << "\n";

 }

 
 
 
 void Radial_Distribution_Function::structure_factor(string filename, int n_ks)
 {
   int timeii, kii, binii;
   float * klist;
   float * structure_fac;
   float min_k;
   float mean_rho=0;
   Coordinate boxsize;
   float boxvolume, rshell;
   
   for(timeii=0;timeii<n_times;timeii++)
   {
     boxsize = system->size(system_time(timeii));
     boxvolume=boxsize.show_x()*boxsize.show_y()*boxsize.show_z();
     mean_rho+=float(n_atoms_j[timeii])/boxvolume/n_times;
   }
   
   min_k = 2.*PI/max_distance;
   
   klist = new float [n_ks];
   structure_fac = new float [n_ks];
   
   for(kii=0;kii<n_ks;kii++)
   {
     klist[kii]=min_k*(kii+1);
     structure_fac[kii]=0;
     for(binii=0;binii<n_bins;binii++)
     {
       rshell = (float(binii)+0.5)*bin_size;						//determine mean radius of bin
       structure_fac[kii]+=(mean_rdf[binii]-1)*rshell*sin(klist[kii]*rshell)/klist[kii]*bin_size;
     }
     structure_fac[kii]=1.+4.*PI*mean_rho*structure_fac[kii];
   }
   
   ofstream output (filename.c_str());

   cout << "\nWriting structure factor to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Structure factor data calculated from g(r) created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(kii=0;kii<n_ks;kii++)
  {
    output << klist[kii] << "\t";
  }
  output << "\n";
  
  for(kii=0;kii<n_ks;kii++)
  {
    output << structure_fac[kii]  <<  "\t";
  }
  
  output << "\n";

  output.close();
   
 }
 
