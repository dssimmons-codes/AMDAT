/*Exptime_Trajectory_List class - stores a list of trajectories that is the same for all times*/
/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Written by David S. Simmons*/

/**/

#ifndef EXPTIME_TRAJECTORY_LIST
#define EXPTIME_TRAJECTORY_LIST

#include <iostream>
#include "trajectory_list.h"
#include "analysis.h"
#include "system.h"

namespace std{

class Exptime_Trajectory_List: public Trajectory_List, public Analysis
{

  public:
    Exptime_Trajectory_List();
    Exptime_Trajectory_List(System*syst, int capacity=0);
    void reset(System*syst, int capacity=0);

    void write_count(string)const;
    
    Analysis_Type what_are_you(){Analysis_Type type = exptime_trajectory_list; return type;};		//virtual method to report the type of analysis
    
    virtual void atomkernel(int species_index, int moleculeii, int atomtype, int atomindex){};

    virtual void displacementkernel(int,int,int,int,int,int,int){};
};

}
#endif
