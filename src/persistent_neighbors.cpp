// Copyright (C) 2010-2025 David S. Simmons

/*Methods for Persistent_Neighbors class - Identifies particles participating in stringlike cooperative rearrangements and converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibody.h"
#include "generated/version.h"
#include "persistent_neighbors.h"
#include "system.h"


using namespace std;



Persistent_Neighbors::Persistent_Neighbors():Dynamic_Cluster_Multibodies()
{
  neighbor_list=0;
  
}


Persistent_Neighbors::Persistent_Neighbors(const Persistent_Neighbors& copy):Dynamic_Cluster_Multibodies(copy)
{
  neighbor_list=copy.neighbor_list;
}



Persistent_Neighbors Persistent_Neighbors::operator=(const Persistent_Neighbors& copy)
{
  if(this!=&copy)
  {
    Dynamic_Cluster_Multibodies::operator=(copy);
    neighbor_list=copy.neighbor_list;
  }
  return *this;
}


Persistent_Neighbors::Persistent_Neighbors(System * syst, int tgap, Neighbor_List* nlist):Dynamic_Cluster_Multibodies(syst,tgap)
{
  neighbor_list=nlist;
}



//checks for persistent neighbors - this is a very inefficient way to do this - it is n^2 when anything based on a neighbor list should be n. If I recoded from the level of dynamic_cluster_multibodies I could make this vastly more efficient by just getting a list of persistent neighbors from the neighbor_list and checking it against the included trajectories at the same time, since the neighbor list is likely to be shorter. Can consider reworking later.
bool Persistent_Neighbors::clustered_check(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  bool check;
  int trajectory1ID;
  
  trajectory1ID=trajectory1->show_trajectory_ID();
  
  check=(neighbor_list->is_neighbor(thisii,trajectory1ID, trajectory2))&&(neighbor_list->is_neighbor(nextii,trajectory1ID, trajectory2));

  return check;
}



Coordinate Persistent_Neighbors::get_imageoffset(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  return (trajectory1->show_coordinate(thisii)).closest_image(trajectory2->show_coordinate(thisii),system->size(thisii));
}

