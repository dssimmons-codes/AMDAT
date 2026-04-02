// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit*/
/*Methods for class to calculate non-gaussian parameter*/
/*Written by David S. Simmons*/

#include "non_gaussian_parameter.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"

using namespace std;

Non_Gaussian_Parameter::Non_Gaussian_Parameter()
{
  system = 0;
  n_times = 0;
  msd = 0;
  
  ngp = new float* [n_times];
  weighting = new long int * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]= new float [16];
    weighting[timeii]= new long int [8];
  }
  timetable = 0;
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii][0]=0;
    weighting[timeii][0]=0;
  }
  atomcount = 0;
}



Non_Gaussian_Parameter::Non_Gaussian_Parameter(const Non_Gaussian_Parameter & copy)
{
  system = copy.system;
  n_times = copy.n_times;
  msd = copy.msd;
  
  ngp = new float * [n_times];
  weighting = new long int * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]= new float [16];
    weighting[timeii]= new long int [8];
  }
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii][0]=copy.ngp[timeii][0];
    weighting[timeii][0]=copy.weighting[timeii][0];
  }
  atomcount = copy.atomcount;
}


Non_Gaussian_Parameter Non_Gaussian_Parameter::operator =(const Non_Gaussian_Parameter & copy)
{
 if(this!=&copy)
 {
  system = copy.system;
  n_times = copy.n_times;
  msd = copy.msd;
  
  ngp = new float * [n_times];
  weighting = new long int * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]= new float [16];
    weighting[timeii]= new long int [8];
  }
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii][0]=copy.ngp[timeii][0];
    weighting[timeii][0]=0;;
  }
  atomcount = copy.atomcount;
 }
 return *this;
}


Non_Gaussian_Parameter::Non_Gaussian_Parameter(System* sys, const Mean_Square_Displacement * m)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();
  msd = m;

  //allocate memory for ngp data and msd data
  ngp = new float * [n_times];
  weighting = new long int * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]= new float [16];
    weighting[timeii]= new long int [8];
  }
  timetable = system -> displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii][0]=0;
    weighting[timeii][0]=0;;
  }
  atomcount = 0;
}


void Non_Gaussian_Parameter::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;
	system->displacement_list(this,false);
	postprocess_list();
}



void Non_Gaussian_Parameter::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	int currenttime=thisii;
	int nexttime=nextii;
	int currenttimegap=timegapii;
	weighting[timegapii][0]+=trajectory_list->show_n_trajectories(thisii);
	(trajectory_list[0]).listloop(this,timegapii,thisii,nextii);
}



void Non_Gaussian_Parameter::listkernel(Trajectory* current_trajectory)
{
//	#pragma omp atomic
//	ngp[currenttimegap]+=pow(current_trajectory->distance(currenttime,nexttime),4);
}


void Non_Gaussian_Parameter::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  double dist = current_trajectory->distance(thisii,nextii);
	ngp[timegapii][0]+=dist*dist*dist*dist;
}


void Non_Gaussian_Parameter::postprocess_list()
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
    double timestepmsd = msd->show(timeii);
		ngp[timeii][0] *= (3.0/(float(weighting[timeii][0])))/(5*timestepmsd*timestepmsd);
		ngp[timeii][0] -= 1.0;
	}
}



void Non_Gaussian_Parameter::write(string filename)const
{
  int timeii;

  cout << "\nWriting non-Gaussian parameter to file " << filename <<".";

  ofstream output(filename.c_str());

  output << "Non-Gaussian parameter data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ngp[timeii][0]<<"\n";
  }
}

void Non_Gaussian_Parameter::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting non-Gaussian parameter to file.";

  output << "Non-Gaussian parameter data created by AMDAT v." << amdat::build::SEMVER << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ngp[timeii][0]<<"\n";
  }
}


int Non_Gaussian_Parameter::max()const
{
  int maxtime = -1; 
  int timeii;
  float maxvalue=0;

  for(timeii=1;timeii<n_times;timeii++)
  {
    if(ngp[timeii][0]>maxvalue)
    {
      maxvalue = ngp[timeii][0];
      maxtime = timeii;
    }
  }

  return maxtime;
}


void Non_Gaussian_Parameter::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);
  
}


void Non_Gaussian_Parameter::postprocess_bins()
{
  postprocess_list();
}

