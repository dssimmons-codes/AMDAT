/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to determines how a trajectory list decays*/
/*Written by Daniel Hunsicker*/

#ifndef TRAJECTORY_LIST_DECAY
#define TRAJECTORY_LIST_DACAY

#include <iostream>
#include "trajectory.h"
#include "boolean_list.h"
#include "trajectory_list.h"
#include "system.h"

namespace std{


class Trajectory_List_Decay: public Analysis
{
//    System* system;
    Boolean_List currently_included;
    Boolean_List initially_included;
    float * avg_decay;
    int current_num;
    int initial_num;
    int total_initial;
    int * total_current;
    int current_time;
    int n_blocks;
    int blocksize;

  public:
    Trajectory_List_Decay();
    Trajectory_List_Decay(System*);
    Trajectory_List_Decay(const Trajectory_List_Decay &);
    Trajectory_List_Decay operator = (const Trajectory_List_Decay &);
    ~Trajectory_List_Decay();

    void analyze(Trajectory_List * t_list);
    void listkernel(Trajectory *);
    void write(string);
    void postprocess();
};


}
#endif
