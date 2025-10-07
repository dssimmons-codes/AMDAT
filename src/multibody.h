// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to allow definition of multi-trajectory bodies other than molecules*/
/*Written by David S. Simmons*/

#ifndef MULTIBODY
#define MULTIBODY

#include <vector>

#include "trajectory.h"
#include "coordinate.h"
#include "trajmath.h"
#include "rgtensor.h"


namespace std
{

  class System;

  class Multibody: public Trajectory
  {
    System * system;

    //int n_trajectories;
    //Trajectory ** trajectories;
    vector<Trajectory*> trajectories;
    int tref;				//reference time at which multibody is defined within consistent image flags
    vector<Coordinate> relative_image_index; //stores relative image flag for each body in the multibody at its reference time of creation so that any calculation involving the relative internal relationships of the body within the multibody correctly compute their relative separation vectors at all times (should always consist of three ints);

    //int multibody_ID;

    //void calculate_mass;
    //void calculate_center_of_mass();
    Coordinate calculate_centroid(int timeii)const;
    
    Coordinate consistent_position(int trajii,int timeii)const;	//returns position of specified trajectory at specified time, using image of each particle such that the multibody does not wrap a boundary. Any computation involving internal degrees of freedom of the multibody should use this method to obtain particle coordinates.
    
    

    public:
      Multibody();
      Multibody(const Multibody &);
      //~Multibody();
      Multibody operator =(const Multibody &);
      Trajectory * operator()(int bodyindex);	//return pointer to trajectory indicated by bodyindex or null pointer if invalid index

      Multibody(System * sys,int reftime=0);
      Multibody(int n_bodies);		//effectively the same as the default constructor
      Multibody(int n_bodies,Trajectory** bodies);

      void set(System * sys);
      void set(int n_bodies, Trajectory** bodies);
      void set(int body_index, Trajectory * body){trajectories[body_index]=body;};	//set pointer to one of the bodies
      //void set(int body_index, Trajectory * body, Coordinate imageoffset)	//set pointer to one of the bodies while setting its image offsets
      
      //Coordinate shortvector(int trajii1,int trajii2,int timeii);		//Returns shortest vector between two trajectories at a given time (first trajectory index, second trajectory index, time)

      //void show_coordinates(int timeii,Coordinate* list);			//returns list of coordinates of trajectories at time ii. Coordinates are returned via Coordinate* list, which must be of length equal to n_trajectories;

      void center_of_mass_trajectory(int trajtype=0);
      void centroid_trajectory(int trajtype=0);

      int show_n_bodies()const{return trajectories.size();};

      float square_gyration_radius(int timeii);
      void gyr_tensor(int timeii, sixfloat*);
      //threefloat principle_axes(int timeii);
      
      bool trajectory_check(Trajectory*);	//returns true if this multibody contains the specified trajectory; returns false otherwise
      void add_body(Trajectory*);
      void set_reftime(int reftime){tref=reftime;};			//set reference time index for correct wrapping/unwrapping of internal degrees of freedom
      void add_body(Trajectory*, Coordinate);
      void absorb_multibody(const Multibody &, Coordinate);
      int show_body_ID(int bodyii){return trajectories[bodyii]->show_trajectory_ID();};	//gives the trajectory ID of specified trajectory in multibody
      void clear();
      
      vector<Trajectory*> show_bodies()const{return trajectories;};

  };

}

#endif
