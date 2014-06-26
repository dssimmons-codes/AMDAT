/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to allow definition of multi-trajectory bodies other than molecules*/
/*Use of this class currently relies upon the assumption that the unwrapped coordinates reflect the true relative position of the trajectories in the multibody. For multibodies consisting of trajectories within a single molecule (all current multibodies), this essentially requires that no molecules crossed the box boundaries at the beginning of the simulation, which may have been well before the start of the simulation trajectory under consideration by AMDAT*/
/*Written by David S. Simmons*/

#ifndef MULTIBODY
#define MULTIBODY

#include "trajectory.h"


namespace std
{

  class System;

  class Multibody
  {
    System * system;
    
    int n_trajectories;
    Trajectory ** trajectories;
    
    //int multibody_ID;
    
    //void calculate_mass;
    //void calculate_center_of_mass();
    //void calculate centroid
    
    public:
      Multibody();
      Multibody(const Multibody &);
      ~Multibody();
      Multibody operator =(const Multibody &);
      Multibody(int n_bodies);		//construct multibody with n_bodies initially corresponding to null pointers
      Multibody(int n_bodies,Trajectory** bodies);
      
      void set(int n_bodies, Trajectory** bodies);
      void set(int body_index, Trajectory * body){trajectories[body_index]=body;};	//set pointer to one of the bodies
      //Coordinate shortvector(int trajii1,int trajii2,int timeii);		//Returns shortest vector between two trajectories at a given time (first trajectory index, second trajectory index, time)

      //void show_coordinates(int timeii,Coordinate* list);			//returns list of coordinates of trajectories at time ii. Coordinates are returned via Coordinate* list, which must be of length equal to n_trajectories;
      
      Trajectory center_of_mass()const;
      Trajectory centroid()const;
      
      //float gyration_radius(int timeii);
      //threefloat principle_axes(int timeii);

  };
  
}

#endif