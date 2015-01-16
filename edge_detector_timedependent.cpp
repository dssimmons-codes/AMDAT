/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class containing information about system composition*/
/*Written by Daniel A. Hunsicker*/


#include <stdlib.h>
#include "edge_detector_timedependent.h"
#include "version.h"
#include <omp.h>

using namespace std;

Edge_Detector_Timedependent::Edge_Detector_Timedependent()
{
    system=0;
    n_times=0;
    max_magnitude=0;
    time_dependent_edge = new Coordinate [n_times];
}

Edge_Detector_Timedependent::Edge_Detector_Timedependent(System* sys, Coordinate vector)
{
  system=sys;
  n_times=system->show_n_timesteps();
  max_magnitude=0;
  unit_vector=vector.unit_vector();
  time_dependent_edge = new Coordinate [n_times];
}



Edge_Detector_Timedependent::Edge_Detector_Timedependent(const Edge_Detector_Timedependent & copy):Analysis(copy)
{
  n_times=copy.n_times;
  unit_vector=copy.unit_vector;
  time_dependent_edge = new Coordinate [n_times];
  max_magnitude=copy.max_magnitude;
  for (int timeii=0; timeii<n_times;timeii++)
  {
      time_dependent_edge[timeii]=copy.time_dependent_edge[timeii];
  }
}


Edge_Detector_Timedependent::~Edge_Detector_Timedependent()
{
  delete [] time_dependent_edge;
}



Edge_Detector_Timedependent Edge_Detector_Timedependent::operator=(const Edge_Detector_Timedependent & copy)
{
  if(this!=&copy)
  {
    system=copy.system;
    tokenize=copy.tokenize;
    trajectory_list=copy.trajectory_list;
    n_times=copy.n_times;
    unit_vector=copy.unit_vector;
    time_dependent_edge = new Coordinate [n_times];
    for (int timeii=0; timeii<n_times;timeii++)
    {
        time_dependent_edge[timeii]=copy.time_dependent_edge[timeii];
    }
    max_magnitude=copy.max_magnitude;
  }
  return *this;
}


void Edge_Detector_Timedependent::analyze(Trajectory_List * t_list)
{
    trajectory_list=t_list;

    for (int timeii=0; timeii<n_times;timeii++)
    {
      if(trajectory_list->show_n_trajectories(timeii)>0)
      {
	max_magnitude=((*trajectory_list)(timeii,0)->show_coordinate(timeii))&unit_vector;
	time_dependent_edge[timeii]=(*trajectory_list)(timeii,0)->show_coordinate(timeii);
      }
      trajectory_list->listloop(this,0,timeii,0);
    }
}


void Edge_Detector_Timedependent::listkernel(Trajectory* traj,int timegapii, int thisii, int nextii)
{
    float current_magnitude = (traj->show_coordinate(thisii))&unit_vector;
    if(current_magnitude>max_magnitude)
    {
      max_magnitude=current_magnitude;
      time_dependent_edge[thisii]=traj->show_coordinate(thisii);
    }
}


void Edge_Detector_Timedependent::write(string filename)
{

    cout << "\nWriting edge to file " << filename << "." << endl;

    ofstream output(filename.c_str());

    output << "Edge data created by AMDAT v." << VERSION << "\n";
    output << "Unit vector is " <<  unit_vector.show_x() << "\t" <<unit_vector.show_y() << "\t" << unit_vector.show_z() << "\n";
    output << "time\tx\ty\tz\n";
    for (int timeii=0; timeii<n_times;timeii++)
    {
      output << system->show_time(timeii) << "\t" << time_dependent_edge[timeii].show_x() << "\t"<<time_dependent_edge[timeii].show_y()<<"\t"<<time_dependent_edge[timeii].show_z()<<"\n";
    }
}

