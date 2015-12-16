/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_square_displacement_2d.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;

Mean_Square_Displacement_2D::Mean_Square_Displacement_2D()
{
  system=0;
  plane="";
  n_times=0;
  atomcount=0;
  msd = new float [n_times];
  weighting = new int [n_times];
  timetable = new float [n_times];
  distancefun = &Trajectory::distance;
}

Mean_Square_Displacement_2D::~Mean_Square_Displacement_2D()
{
  delete [] msd;
  delete [] weighting;
  delete [] timetable;
}

Mean_Square_Displacement_2D::Mean_Square_Displacement_2D(const Mean_Square_Displacement_2D & copy)
{
  system = copy.system;
  plane = copy.plane;
  n_times = copy.n_times;
  atomcount = copy.atomcount;
  distancefun = copy.distancefun;
  trajectory_list=copy.trajectory_list;
  
  msd = new float [n_times];
  weighting = new int [n_times];
  timetable = system->displacement_times();
  
  for (int timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=copy.msd[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}

Mean_Square_Displacement_2D Mean_Square_Displacement_2D::operator=(const Mean_Square_Displacement_2D & copy)
{
  if(this!=&copy)
  {
  
  system = copy.system;
  plane = copy.plane;
  n_times = copy.n_times;
  atomcount = copy.atomcount;
  distancefun = copy.distancefun;
  trajectory_list = trajectory_list;
  
  delete [] msd;
  delete [] weighting;
  delete [] timetable;
    
  msd = new float [n_times];
  weighting = new int [n_times];
  timetable = system->displacement_times();
  
  for (int timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=copy.msd[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  
  }
  return *this;
}

Mean_Square_Displacement_2D::Mean_Square_Displacement_2D(System*sys,string orientation)
{
  initialize(sys,orientation);
}


void Mean_Square_Displacement_2D::initialize(System* sys,string orientation)
{
  int timeii;
	
  system = sys;
  plane = orientation;
  
  if(plane=="xy") distancefun = &Trajectory::distance_xy;
  else if(plane=="xz") distancefun = &Trajectory::distance_xz;
  else if(plane=="yz") distancefun = &Trajectory::distance_yz;
  else
  {
    cout<<"Error: plane command "<<orientation<<" not understood.\n";
    exit(1);
  }
  
  n_times = system->show_n_timegaps();
	
   //allocate memory for mean square displacement data
  msd = new float [n_times];
  weighting =system->timegap_weighting();
  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    msd[timeii]=0;
  }
  atomcount = 0;
}



/*--------Methods to do analysis using trajectory lists------------*/

void Mean_Square_Displacement_2D::analyze(Trajectory_List * t_list)
{
  int timeii;
  trajectory_list=t_list;
  for(timeii=0;timeii<n_times;timeii++)
  {
	  weighting[timeii]=0;
  }
  system->displacement_list(this);
  postprocess_list();
}

void Mean_Square_Displacement_2D::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	currenttimegap=timegapii;
	weighting[timegapii]+=trajectory_list[0].show_n_trajectories(currenttime);
	(trajectory_list[0]).listloop(this,currenttime);
}

void Mean_Square_Displacement_2D::listkernel(Trajectory* current_trajectory)
{
	msd[currenttimegap]+=pow((current_trajectory->*distancefun)(currenttime,nexttime),2);
}

void Mean_Square_Displacement_2D::postprocess_list()
{
	int timeii;
	for(timeii=0;timeii<n_times;timeii++)
	{
		msd[timeii] /= float(weighting[timeii]);
	}
}

/*--------Methods to do analysis using binned trajectory lists------------*/

void Mean_Square_Displacement_2D::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;
  list_displacementkernel(timegapii, thisii, nextii);
}

void Mean_Square_Displacement_2D::postprocess_bins()
{
  postprocess_list();
}

/*----------Method to write MSD data to file-------------*/

void Mean_Square_Displacement_2D::write(string filename)
{
  int timeii;
  
  cout << "\nWriting msd to file.";
  
  ofstream output(filename.c_str());
  
  output << "2-D mean square displacement data for "<< plane <<" plane created by MDAT v." << VERSION << "\n"; 
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<msd[timeii]<<"\n";
  }
}

void Mean_Square_Displacement_2D::write(ofstream& output)const
{
  int timeii;
  
  cout << "\nWriting 2D msd to file.";
  
  output << "2-D mean square displacement data for "<< plane <<" plane created by MDAT v." << VERSION << "\n"; 
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<msd[timeii]<<"\n";
  }
}