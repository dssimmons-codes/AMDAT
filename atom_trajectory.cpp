/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for atom_trajectory class*/
/*Written by David S. Simmons*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atom_trajectory.h"
#include <iostream>


using namespace std;



Atom_Trajectory::Atom_Trajectory()
{
  trajtype=atom;
  is_unwrapped=0;
  is_wrapped=0;
  is_velocity=0;
  n_timesteps=0;
  type=-1;
  trajectory_ID=-1;
  mass=1;
  coordinates = new Coordinate[1];
  unwrapped = new Coordinate[1];
  velocity=new Coordinate[1];
  atomID=-1;
  moleculeID=-1;
}

/*Atom_Trajectory Constructor to initialize timestep count and allocate memory accordingly; defaults timesteps to one*/
Atom_Trajectory::Atom_Trajectory(int timesteps)
{
  trajtype=atom;
  n_timesteps = timesteps;					//set number of timesteps (size of trajectory coordinate array)
  coordinates = new Coordinate[n_timesteps];			//allocate memory for coordinate array
  unwrapped = new Coordinate[1];				//allocate memory for displacement array
  velocity = new Coordinate[1];
  is_unwrapped = 0;						//not unwrapped
  atomID = NULL;
  trajectory_ID = atomID;
  mass = 1;							//set default mass to 1
}



/*Atom_Trajectory Constructor to instantiate coordinate array*/
Atom_Trajectory::Atom_Trajectory(int timesteps, Coordinate * coordinatelist)
{
	trajtype=atom;
  n_timesteps = timesteps;					//set number of timesteps (size of trajectory coordinate array)
  coordinates = new Coordinate[n_timesteps];			//allocate memory for coordinate array
  coordinates = coordinatelist;					//define coordinate array
  unwrapped = new Coordinate[1];
  velocity = new Coordinate[1];
  is_unwrapped = 0;						//not unwrapped
  atomID = NULL;
  trajectory_ID = atomID;
  mass = 1;							//set default mass to 1
}

/*--------------------------------------------------------------------------------------*/

/*Method to define trajectory size and allocate memory appropriately*/
void Atom_Trajectory::set(int atom_type,int timesteps, int m)
{
	clear_memory();					//clear memory assigned to existing arrays so they may be reassigned

	type = atom_type;					//set parameters
	n_timesteps=timesteps;
	trajectory_ID = NULL;						
	atomID = trajectory_ID;
	
	coordinates = new Coordinate[n_timesteps];		//allocate new memory to coordinate array

	unwrapped = new Coordinate[1];			//dummy memory allocation for unwrapped coordinate array
	velocity = new Coordinate[1];
	mass = m;
}

/*--------------------------------------------------------------------------------------*/

/*Method to define coordinates and associated paramters*/
void Atom_Trajectory::set(int atom_type,int timesteps, Coordinate* coordinatelist, int m =1)
{
	clear_memory();					//clear memory assigned to existing arrays so they may be reassigned

	type = atom_type;					//set parameters
	n_timesteps=timesteps;
	trajectory_ID = NULL;	
	atomID = trajectory_ID;
  
	coordinates = new Coordinate[n_timesteps];		//allocate new memory to coordinate array
	coordinates = coordinatelist;				//copy coordinate array

	unwrapped = new Coordinate[1];			//dummy memory allocation for unwrapped coordinate array
	velocity = new Coordinate[1];
}
	
/*--------------------------------------------------------------------------------------*/





