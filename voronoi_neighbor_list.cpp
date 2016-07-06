/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Voronoi_Neighbor_List class - builds neighbor list from distance thresholding*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "voronoi_neighbor_list.h"
#include "voro++.hh"
#include "container.hh"

using namespace std;



Voronoi_Neighbor_List::Voronoi_Neighbor_List():Neighbor_List()
{
  system=0;
  n_times=1;
}



Voronoi_Neighbor_List::Voronoi_Neighbor_List(const Voronoi_Neighbor_List& copy):Neighbor_List(copy)
{
  int typeii;
 
  system=copy.system;
    
}




Voronoi_Neighbor_List Voronoi_Neighbor_List::operator=(const Voronoi_Neighbor_List& copy)
{
  int typeii;
  
  if(this!=&copy)
  {

    Neighbor_List::operator=(copy);
    
    system=syst;
  }

  return *this;
}


Voronoi_Neighbor_List::~Voronoi_Neighbor_List()
{

}


Voronoi_Neighbor_List::Voronoi_Neighbor_List(System* sys, int firsttime, int lasttime):Neighbor_List(sys)
{
  int startoffset=firsttime;
  int endoffset;
 
  system=sys;
  computed_times.clear();
  
  if(firsttime==-1)
  {
    computed_times.clear();
    computed_times.resize(n_times,true);
  }
  else
  {
    if(lasttime==-1)
    {
      endoffset=firsttime;
    }
    else
    {
      endoffset=lasttime;
    }
  
    int n_blocks = system->show_n_exponentials();
    int blocksize = system->show_n_exponential_steps();
    int blockstart;
    
    computed_times.resize(n_times,false);
  
    if(firsttime>blocksize||lasttime>blocksize)
    {
      cout<<"\nError: time offsets selected for neighbor list construction must fit within the block.\n";
      exit(0);
    }
  
    for(int blockii=0;blockii<n_blocks;blockii++)
    {
      blockstart = blockii*blocksize;
      computed_times[blockstart]=true;
      for(int offsetii=startoffset; offsetii<=endoffset;offsetii++)
      {
	computed_times[blockstart+offsetii]=true;
      }
    }
  }
  
}



void Voronoi_Neighbor_List::timekernel2(int timeii)
{
  float x_min,x_max,y_min,y_max,z_min,z_max;
  int trajii, neighii;
  int n_x,n_y,n_z;
  int n_trajectories;
  vector<vector<int>> voronoi_neighbors;
  float n_traj;
  
  if(computed_times[timeii])
  {
    n_traj=trajectory_list->show_n_trajectories(timeii);
    
    n_x=n_y=n_z=int(pow(n_traj/8.0,1.0/3.0));
 
    n_trajectories=system->show_n_trajectories();
    
    x_min=(system->boundaries())[0].show_x();
    y_min=(system->boundaries())[0].show_y();
    z_min=(system->boundaries())[0].show_z();
    x_max=(system->boundaries())[1].show_x();
    y_max=(system->boundaries())[1].show_y();
    z_max=(system->boundaries())[1].show_z();
    
    voronoi = new voro::container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,true,true,true,8);
		  
    
   
    trajectory_list->listloop(this,0, timeii, 0);
    
    voronoi->print_custom("%i %n","packing.custom2");
    
    voronoi_neighbors=voronoi->get_neighList();
    
    for(trajii=0;trajii<n_trajectories;trajii++)
    {
      if((included[timeii])(trajii))
      {
	for(neighii=0;neighii<voronoi_neighbors[trajii].size();neighii++)
	{
	  neighbors[timeii][trajii].push_back(system->show_trajectory(voronoi_neighbors[trajii][neighii]));
	  
	}
      }
    }
    
    delete voronoi;
  }
}

void Voronoi_Neighbor_List::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  int trajectory1ID;
  Coordinate position;
  
  trajectory1ID=current_trajectory->show_trajectory_ID();
  (included[thisii])(trajectory1ID,1);
  position=current_trajectory->show_coordinate(thisii);
  
  voronoi->put(trajectory1ID,position.show_x(),position.show_y(),position.show_z());
 
}

void Voronoi_Neighbor_List::postprocess_list()
{
  values.resize(neighbors.size());
  int timeii,trajii;
  
  for(timeii=0;timeii<neighbors.size();timeii++)
  {
    if(computed_times[timeii])
    {
      values[timeii].resize(neighbors[timeii].size());
      for(trajii=0;trajii<neighbors[timeii].size();trajii++)
      {
	if(included[timeii](trajii))
	{
	  values[timeii][trajii]=neighbors[timeii][trajii].size();
	}
      }
    }
  }
}