/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "mean_displacement.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "static_trajectory_list.h"
using namespace std;


Mean_Displacement::Mean_Displacement()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  md = new Coordinate [n_times];
  weighting = new int [n_times];

  atomcount = 0;
}


Mean_Displacement::Mean_Displacement(const Mean_Displacement & copy)
{
  int timeii;

  system = copy.system;
  trajectory_list = copy.trajectory_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii]=copy.md[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}



/** **/
Mean_Displacement::Mean_Displacement(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  md = new Coordinate [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
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

  delete [] md;
  delete [] weighting;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    md[timeii]=copy.md[timeii];
    weighting[timeii]=copy.weighting[timeii];
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

  delete [] md;
  delete [] weighting;

  md = new Coordinate [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
  }
  atomcount = 0;
}

void Mean_Displacement::preprocess()
{
  weighting = system->timegap_weighting();
}

/*Methods to do analysis using system loops*/



void Mean_Displacement::displacementkernel(int timegap,int thisii, int nextii, Trajectory* traj)
{
  md[timegap]+=traj->show_unwrapped(nextii)-traj->show_unwrapped(thisii);
}



void Mean_Displacement::atomkernel(Trajectory * traj)
{
	system->displacement_loop(this, traj);
	atomcount ++;
}


void Mean_Displacement::postprocess()
{

  for(int timeii=0;timeii<n_times;timeii++)
  {

        md[timeii] /= (float(weighting[timeii])*float(atomcount));

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

  currenttime=thisii;
  nexttime=nextii;
  currenttimegap=timegapii;

  weighting[timegapii]+=trajectory_list->show_n_trajectories(currenttime);
  //weighting[timegapii]+=(trajectory_list[0]).show_n_trajectories(currenttime);
  (trajectory_list[0]).listloop(this,currenttime);
}



void Mean_Displacement::listkernel(Trajectory* current_trajectory)
{
  md[currenttimegap]+=current_trajectory->show_unwrapped(nexttime)-current_trajectory->show_unwrapped(currenttime);
}



void Mean_Displacement::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        md[timeii] /= float(weighting[timeii]);

  }
}



/*Method to write MD data to file*/

void Mean_Displacement::write(string filename)
{
  int timeii;

  cout << "\nWriting mean displacement to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Mean displacement data created by AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<md[timeii].show_x()<<"\t"<<md[timeii].show_y()<<"\t"<<md[timeii].show_z()<<"\n";
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
