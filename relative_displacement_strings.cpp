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
  
}


Relative_Displacement_Strings::Relative_Displacement_Strings(const Relative_Displacement_Strings& copy):Dynamic_Cluster_Multibodies(copy)
{
  neighbor_list=copy.neighbor_list;
}



Relative_Displacement_Strings Relative_Displacement_Strings::operator=(const Relative_Displacement_Strings& copy)
{
  if(this!=&copy)
  {
    Dynamic_Cluster_Multibodies::operator=(copy);
    neighbor_list=copy.neighbor_list;
  }
  return *this;
}


Relative_Displacement_Strings::Relative_Displacement_Strings(System * syst, int tgap, Neighbor_List* nlist, float thresh):Dynamic_Cluster_Multibodies(syst,tgap)
{
  threshold=thresh;
  neighbor_list=nlist;
}




bool Relative_Displacement_Strings::clustered_check(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  bool check;
  float initial_separation,distance;
  int trajectory1ID;
  
  trajectory1ID=trajectory1->show_trajectory_ID();
  
  check = (neighbor_list->is_neighbor(thisii,trajectory1ID, trajectory2));
  
  if(check)
  {
    initial_separation=(trajectory2->show_coordinate(thisii)-trajectory1->show_coordinate(thisii)).length_unwrapped(system->size(thisii));
    distance = (trajectory2->show_coordinate(thisii)-trajectory1->show_coordinate(nextii)).length_unwrapped(system->size(thisii));
    check=check&&(distance<(threshold*initial_separation));
  }
  
  return check;
}



Coordinate Relative_Displacement_Strings::get_imageoffset(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  return (trajectory1->show_coordinate(thisii)).closest_image(trajectory2->show_coordinate(thisii),system->size(thisii));
}

