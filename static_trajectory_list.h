/*Static_Trajectory_List class - stores a list of trajectories that is the same for all times*/
/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Written by David S. Simmons*/



#ifndef STATIC_TRAJECTORY_LIST
#define STATIC_TRAJECTORY_LIST

#include <iostream>
#include "trajectory_list.h"
#include "analysis.h"
#include "system.h"
#include "trajectory_set.h"

namespace std{

class Static_Trajectory_List: public Trajectory_List, public Analysis
{

  public:
    Static_Trajectory_List();
    Static_Trajectory_List(System*syst, int capacity=0);
    void reset(System*syst, int capacity=0);

    void atomkernel(Trajectory * traj);
    
    void set(System * syst, Trajectory_Set * trajectory_set);		//initialize trajectory list based on trajectory set
};

}
#endif
