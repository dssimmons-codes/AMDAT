/*Molecular Dynamics Analysis Toolkit*/
/*Methods for class to calculate non-gaussian parameter*/
/*Written by David S. Simmons*/

#include "non_gaussian_parameter.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;


Non_Gaussian_Parameter::Non_Gaussian_Parameter()
{
  system = 0;
  n_times = 0;
  msd = 0;
  
  ngp = new float [n_times];
  n_atoms = new float [n_times];
  weighting = 0;
  timetable = 0;
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=0;
    n_atoms[timeii]=0;
  }
  atomcount = 0;
}



Non_Gaussian_Parameter::Non_Gaussian_Parameter(const Non_Gaussian_Parameter & copy)
{
  system = copy.system;
  n_times = copy.n_times;
  msd = copy.msd;
  
  ngp = new float [n_times];
  n_atoms = new float [n_times];
  weighting = system -> timegap_weighting();
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=copy.ngp[timeii];
    n_atoms[timeii]=copy.ngp[timeii];
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
  
  ngp = new float [n_times];
  n_atoms = new float [n_times];
  weighting = system -> timegap_weighting();
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=copy.ngp[timeii];
    n_atoms[timeii]=copy.ngp[timeii];
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
  ngp = new float [n_times];
  n_atoms = new float [n_times];
  weighting = system -> timegap_weighting();
  timetable = system -> displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=0;
    n_atoms[timeii]=0;
  }
  atomcount = 0;
}


void Non_Gaussian_Parameter::displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj)
{
  ngp[timegap]+=pow(traj->distance(thisii,nextii),4);
}



void Non_Gaussian_Parameter::atomkernel(Trajectory * traj)
{
  system->displacement_loop(this, traj);
  atomcount ++;
}


void Non_Gaussian_Parameter::postprocess()
{
  int timeii;

  for(timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii] *= (3.0/(float(weighting[timeii])*float(atomcount)))/(5.0*pow((msd->show(timeii)),2.0));
    ngp[timeii] -= 1.0;
  }
}



void Non_Gaussian_Parameter::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;

	system->displacement_list(this);

	postprocess_list();
}



void Non_Gaussian_Parameter::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	currenttimegap=timegapii;
	n_atoms[timegapii]+=trajectory_list[0].show_n_trajectories(currenttime)/float(weighting[timegapii]);
	(trajectory_list[0]).listloop(this,currenttime);
}



void Non_Gaussian_Parameter::listkernel(Trajectory* current_trajectory)
{
	ngp[currenttimegap]+=pow(current_trajectory->distance(currenttime,nexttime),4);

}



void Non_Gaussian_Parameter::postprocess_list()
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
		ngp[timeii] *= (3.0/(float(weighting[timeii])*n_atoms[timeii]))/(5*pow((msd->show(timeii)),2.0));
		ngp[timeii] -= 1.0;
	}
}



void Non_Gaussian_Parameter::write(string filename)const
{
  int timeii;

  cout << "\nWriting non-Gaussian parameter to file " << filename <<".";

  ofstream output(filename.c_str());

  output << "Non-Gaussian parameter data created by MDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ngp[timeii]<<"\n";
  }
}



int Non_Gaussian_Parameter::max()const
{
  int maxtime, timeii;
  float maxvalue=0;

  for(timeii=1;timeii<n_times;timeii++)
  {
    if(ngp[timeii]>maxvalue)
    {
      maxvalue = ngp[timeii];
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

