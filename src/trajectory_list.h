// Copyright (C) 2010-2025 David S. Simmons

/*Trajectory_List class - stores a list of trajectories*/
/*Amorphous Molecular Dynamics Analysis Toolkit (aMDAT)*/
/*Written by David S. Simmons*/

#ifndef TRAJECTORY_LIST
#define TRAJECTORY_LIST

#include <iostream>
#include "trajectory.h"
#include "boolean_list.h"
#include "trajectory_set.h"
#include "multibody_list.h"

namespace std{

class System;
class Analysis;

class Trajectory_List
{
  protected:
    System * sys;
    Trajectory*** trajectories;		//list of pointers to Trajectory objects, with a separate list at each time
    int * time_conversion;	//table of system frame index to internal time index conversions
    int * n_trajectories;

    int n_times;		//number of internal times in the Trajectory list

    int capacity;		//maximum number of trajectories that can be stored at each time

    mutable Boolean_List * included;	//stores boolean list specifying included trajectories at each time;

    int convert_time(int timeii)const{return time_conversion[timeii];};	//convert requested time (Where the index is the time index from the system object) to internal time index


    int n_atomtypes;

    void trajlist_from_boollist();


    int n_system_trajectories()const;

    void addtrajectory(int, Trajectory*);

  public:
    Trajectory_List(int timecount = 1, int cap = 0);
    Trajectory_List(System* sys, int timecount, int cap, Boolean_List * boollist, int*time_conversion);
    ~Trajectory_List();
    Trajectory_List(const Trajectory_List &); // MEM - copy constructor

    void set(System* sys, int timecount, int cap, Boolean_List * boollist, int*time_conv);
    void set(System * syst, vector<Trajectory_Set*> trajectory_sets, int*time_conv);		//initialize trajectory list based on vector of trajectory sets
    
    void flatten_multibodies(const Multibody_List& mblist);		//set up trajectory list by combining all the trajectories in all multibodies
    
    int show_n_trajectories(int timeii) {return n_trajectories[convert_time(timeii)];}
  
    Trajectory* operator () (int trajii);				//return a requested trajectory at the first time stored
    Trajectory* operator () (int timeii, int trajii);		//return a requested trajectory at a given time
    bool is_included(int timeii,int trajii);                              //returns 1 if trajectory is included at that time
    Boolean_List show_included(int timeii){return included[timeii];}
    void listloop(Analysis* analysis, int time);			//loop over trajectories at a given time
    void listloop(Analysis* analysis, int timegap, int curTime, int nextTime);			//loop over trajectories at a given time
    void listloop2(Analysis* analysis, Trajectory* traj, int timegap, int curTime, int nextTime);			//loop over trajectories at a given time
    int show_n_trajectories(int timeii)const{return n_trajectories[convert_time(timeii)];};	//return number of trajectories at a given time

    virtual void write_count(string)const;
    void write_xyz(string)const;
    void write_full_xyz(string)const;

    Trajectory_List time_intersection(int, int)const;
    Trajectory_List time_union(int,int)const;

    Trajectory_List operator&& (const Trajectory_List &)const;
    Trajectory_List same_timetable_intersection(const Trajectory_List &)const;
    Trajectory_List operator|| (const Trajectory_List &)const;
    Trajectory_List operator = (const Trajectory_List &);

    void inversion(Trajectory_List*,Trajectory_List*);
};

}
#endif
