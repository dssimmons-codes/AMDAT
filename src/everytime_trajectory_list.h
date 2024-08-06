/*Everytime_Trajectory_List class - stores a list of trajectories that is the same for all times*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons*/

/**/

#ifndef EVERYTIME_TRAJECTORY_LIST
#define EVERYTIME_TRAJECTORY_LIST

#include <iostream>
#include "trajectory_list.h"
#include "analysis_onetime.h"
#include "system.h"

namespace std{

class Everytime_Trajectory_List: public Trajectory_List, public Analysis_Onetime
{

  public:
    Everytime_Trajectory_List();
    Everytime_Trajectory_List(System*syst, int capacity=0);
    void reset(System*syst, int capacity=0);

    void write_count(string)const;

    Analysis_Type what_are_you(){Analysis_Type type = everytime_trajectory_list; return type;};		//virtual method to report the type of analysis

    virtual void atomkernel(int species_index, int moleculeii, int atomtype, int atomindex){};

    virtual void displacementkernel(int,int,int,int,int,int,int){};
};

}
#endif
