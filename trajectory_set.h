/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Header for Trajectory_set class - Class to store sets of user-defined trajectories
/*Written by David S. Simmons*/

#ifndef TRAJECTORY_SET
#define TRAJECTORY_SET

#include "trajectory.h"

namespace std
{
  
class Trajectory_Set
{
    int n_trajectories;
    Trajectory * trajectories;
  
  public:
    Trajectory_Set();
};

}