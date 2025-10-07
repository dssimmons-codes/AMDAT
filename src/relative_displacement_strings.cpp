// Copyright (C) 2010-2025 David S. Simmons

/*Methods for Relative_Displacement_Strings class - Identifies particles participating in stringlike cooperative rearrangements and converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibody.h"
#include "version.h"
#include "relative_displacement_strings.h"
#include "system.h"


using namespace std;




Relative_Displacement_Strings::Relative_Displacement_Strings():Dynamic_Cluster_Multibodies()
{
  neighbor_list=0;
  steps_for_averaging=0;
  threshold=0;
  
}


Relative_Displacement_Strings::Relative_Displacement_Strings(const Relative_Displacement_Strings& copy):Dynamic_Cluster_Multibodies(copy)
{
  neighbor_list=copy.neighbor_list;
  steps_for_averaging=copy.steps_for_averaging;
  threshold=copy.threshold;
}



Relative_Displacement_Strings Relative_Displacement_Strings::operator=(const Relative_Displacement_Strings& copy)
{
  if(this!=&copy)
  {
    Dynamic_Cluster_Multibodies::operator=(copy);
    neighbor_list=copy.neighbor_list;
    steps_for_averaging=copy.steps_for_averaging;
    threshold=copy.threshold;
  }
  return *this;
}


Relative_Displacement_Strings::Relative_Displacement_Strings(System * syst, int tgap, Neighbor_List* nlist, float thresh, int avgsteps):Dynamic_Cluster_Multibodies(syst,tgap)
{
  threshold=thresh;
  neighbor_list=nlist;
  steps_for_averaging=avgsteps;
}




bool Relative_Displacement_Strings::clustered_check(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  bool check;
  float initial_separation,distance;
  int trajectory1ID;
  int timeii;
  
  trajectory1ID=trajectory1->show_trajectory_ID();
  
  check = (neighbor_list->is_neighbor(thisii,trajectory1ID, trajectory2));
  
  if(check)
  {
    initial_separation=0;
    //take average initial distance over range of time
    for(timeii=0;timeii<steps_for_averaging;timeii++)
    {
      initial_separation+=(trajectory2->show_coordinate(thisii+timeii)-trajectory1->show_coordinate(thisii+timeii)).length_unwrapped(system->size(thisii+timeii));
      //cout<<"\t"<<initial_separation;
    }
    initial_separation/=float(steps_for_averaging);
    
    
    distance = (trajectory2->show_coordinate(thisii)-trajectory1->show_coordinate(nextii)).length_unwrapped(system->size(thisii));
    check=check&&(distance<(threshold*initial_separation));
  }
  
  return check;
}



Coordinate Relative_Displacement_Strings::get_imageoffset(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  return (trajectory1->show_coordinate(thisii)).closest_image(trajectory2->show_coordinate(thisii),system->size(thisii));
}

