/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for Class Trajectory: Parent class for classes that store trajectory information, such as atoms and molecules*/
/*Written by David S. Simmons*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "trajectory.h"

using namespace std;

/*--------------------------------------------------------------------------------------*/

Trajectory::Trajectory()
{
  trajtype=general;
  is_unwrapped=0;
  is_wrapped=0;
  is_velocity=0;
  n_timesteps=0;
  type=-1;
  trajectory_ID=-1;
  mass=-1;
  coordinates = new Coordinate[1];
  unwrapped = new Coordinate[1];
  velocity=new Coordinate[1];
}


/*Trajectory Constructor to initialize timestep count and allocate memory accordingly; defaults timesteps to one*/
Trajectory::Trajectory(int timesteps)
{
	trajtype=general;
	n_timesteps = timesteps;					//set number of timesteps (size of trajectory coordinate array)
	coordinates = new Coordinate[n_timesteps];			//allocate memory for coordinate array
	unwrapped = new Coordinate[1];					//allocate memory for displacement array
	velocity = new Coordinate[1];
	
	for(int timeii=0;timeii<n_timesteps;timeii++)
	{
	  coordinates[timeii].set(0,0,0);
	}
	
	is_unwrapped = 0;						//not unwrapped
	is_wrapped = 1;
	is_velocity = 0;
	trajectory_ID = NULL;
	mass = 1;							//set default mass to 1
}

/*--------------------------------------------------------------------------------------*/

/*Trajectory Constructor to instantiate coordinate array*/
Trajectory::Trajectory(int timesteps, Coordinate * coordinatelist)
{
	trajtype=general;
	n_timesteps = timesteps;					//set number of timesteps (size of trajectory coordinate array)
	coordinates = new Coordinate[n_timesteps];			//allocate memory for coordinate array
	
	for(int timeii=0;timeii<n_timesteps;timeii++)
	{
	  coordinates[timeii] = coordinatelist[timeii];					//define coordinate array
	}
	unwrapped = new Coordinate[1];
	velocity = new Coordinate[1];
	is_unwrapped = 0;						//not unwrapped
	is_wrapped = 1;
	is_velocity = 0;
	trajectory_ID = NULL;
	mass = 1;							//set default mass to 1
}

/*--------------------------------------------------------------------------------------*/


Trajectory::~Trajectory()
{
	delete [] coordinates;
	delete [] unwrapped;
	delete [] velocity;
}

/*--------------------------------------------------------------------------------------*/

/*Method to clear memory allocated to arrays*/
void Trajectory::clear_memory()
{
	delete [] coordinates;
	delete [] unwrapped;
	delete [] velocity;
	is_unwrapped=0;
	is_wrapped=0;
	is_velocity=0;
}

/*--------------------------------------------------------------------------------------*/

/*Method to define trajectory size and allocate memory appropriately*/
void Trajectory::set(int atom_type,int timesteps, float m)
{
	clear_memory();					//clear memory assigned to existing arrays so they may be reassigned

	n_timesteps=timesteps;
	trajectory_ID = NULL;						

	coordinates = new Coordinate[n_timesteps];			//allocate memory for coordinate array
	unwrapped = new Coordinate[1];				//allocate memory for displacement array
	velocity = new Coordinate[1];
	
	for(int timeii=0;timeii<n_timesteps;timeii++)
	{
	  coordinates[timeii].set(0,0,0);
	}
	
	is_unwrapped = 0;						//not unwrapped
	is_wrapped = 1;
	is_velocity = 0;

	mass = m;
}

/*--------------------------------------------------------------------------------------*/

/*Method to define coordinates and associated paramters*/
void Trajectory::set(int atom_type,int timesteps, Coordinate* coordinatelist, int m =1)
{
	clear_memory();					//clear memory assigned to existing arrays so they may be reassigned

	n_timesteps=timesteps;
	
	coordinates = new Coordinate [n_timesteps];
	for(int timeii=0;timeii<n_timesteps;timeii++)
	{
	  coordinates[timeii] = coordinatelist[timeii];					//define coordinate array
	}					//define coordinate array
	unwrapped = new Coordinate[1];
	velocity = new Coordinate[1];
	is_unwrapped = 0;						//not unwrapped
	is_wrapped = 1;
	is_velocity = 0;
	trajectory_ID = NULL;
	mass = 1;	
	
}
	
/*--------------------------------------------------------------------------------------*/
	
/*Method to set single coordinate in coordinate list*/
void Trajectory::set(const Coordinate & coordinate, int timestep)
{
	if(timestep>=n_timesteps)
	{
		cout << "\nTimestep ("<<timestep<<") for atom trajectory greater than number of configurations ("<<n_timesteps<<")!\n";
		exit(1);
	}
	coordinates[timestep]=coordinate;
  
}


/*--------------------------------------------------------------------------------------*/
	
/*Method to set single unwrapped coordinate in coordinate list*/
void Trajectory::set_unwrapped(const Coordinate & coordinate, int timestep)
{
	if(!is_unwrapped)
	{
	  delete [] unwrapped;
	  unwrapped = new Coordinate [n_timesteps];
	  is_unwrapped=true;
	}	
	if(timestep>=n_timesteps)
	{
		cout << "\nTimestep ("<<timestep<<") for atom trajectory greater than number of configurations ("<<n_timesteps<<")!\n";
		exit(1);
	}
	unwrapped[timestep]=coordinate;
  
}


/*--------------------------------------------------------------------------------------*/
	
/*Method to set single velocity in velocity list*/
void Trajectory::set_velocity(const Coordinate & coordinate, int timestep)
{
	if(!is_velocity)
	{
	  delete [] velocity;
	  velocity = new Coordinate [n_timesteps];
	  is_velocity=true;
	}
	if(timestep>=n_timesteps)
	{
		cout << "\nTimestep ("<<timestep<<") for atom trajectory greater than number of configurations ("<<n_timesteps<<")!\n";
		exit(1);
	}
	velocity[timestep]=coordinate;
  
}



/*--------------------------------------------------------------------------------------*/

/*Method to calculate an atom displacement between two time steps.  The time coordinate calculated by this operation is the time between timesteps*/
Coordinate Trajectory::unwrap_displacement(int timestep1, int timestep2,const Coordinate & box_size)
{
	Coordinate displacement;					//instantiate displacement variable

  //Naively calculate displacement, ignoring boundaries
	displacement = coordinates[timestep2] - coordinates[timestep1];
  
  //now adjust for moves across periodic boundaries
	if((displacement.show_x()<box_size.show_x()*.6&&displacement.show_x()>box_size.show_x()*.4)||(displacement.show_y()<box_size.show_y()*.6&&displacement.show_y()>box_size.show_y()*.4)||(displacement.show_z()<box_size.show_z()*.6&&displacement.show_z()>box_size.show_z()*.4))  //check for excessively large real displacements
	{
		cout << "\nWarning: ambiguous particle displacement direction found in unwrapping of trajectory " << trajectory_ID << " due to large displacement magnitude.  Timestep between configurations"<< timestep1 <<" and " << timestep2 << " may be too large.";
	}
	displacement -= box_size*((displacement/(box_size*.55)).integer()); //divisor is .55 rather than .5 in order to prevent it from moving the particle by two box lengths in the case where LAMMPS has allowed a particle to slip outside of the box boundary, potentially yielding a displacement slightly larger than the box size.
	return displacement ;
}

/*--------------------------------------------------------------------------------------*/

/*Method to generate unwrapped coordinate array.*/
void Trajectory::unwrap(const Coordinate & box_size)
{
	int ii;					//index over timestepss
	if(!is_unwrapped)				//only unwrap if unwrapped coordinates do not yet exist
	{
		delete [] unwrapped;			//delete null memory assignment
		unwrapped = new Coordinate[n_timesteps];	//and allocate necessary memory for unwrapped coordinates
		unwrapped[0]=coordinates[0];    //first unwrapped coordinate is same as wrapped coordinate

		for(ii=1;ii<n_timesteps;ii++)		//build unwrapped coordinate set
		{
			unwrapped[ii] = unwrapped[ii-1] + unwrap_displacement(ii-1,ii,box_size);
			//cout <<"\n"<<trajectory_ID<<"\t"<< unwrapped[ii].show_x()<<"\t"<< unwrapped[ii].show_y()<<"\t"<< unwrapped[ii].show_z();
		}
		is_unwrapped = 1;			//note unwrapped state
	}
}





/*--------------------------------------------------------------------------------------*/



/*Method to generate wrapped coordinate array.*/
void Trajectory::wrap(const Coordinate * box_size,Coordinate ** box_boundaries)
{
	int timeii;				//index over timesteps
	if(!is_wrapped)				//only wrap if wrapped coordinates do not yet exist
	{
		delete [] coordinates;			//delete null memory assignment
		coordinates = new Coordinate[n_timesteps];	//and allocate necessary memory for wrapped coordinates

		for(timeii=0;timeii<n_timesteps;timeii++)		//loop over time to build wrapped coordinate set
		{
		  coordinates[timeii]-=box_size[timeii]*((coordinates[timeii]-box_boundaries[timeii][0])/box_size[timeii]).coord_floor();
		}
		is_wrapped = 1;			//note wrapped state
	}
}





/*--------------------------------------------------------------------------------------*/


/*Method to return list of coordinates at requested times.*/
Coordinate * Trajectory::show_coordinates(int*timelist, int listsize)const
{
	int timeii;
	Coordinate* coordinatelist;
	coordinatelist = new Coordinate [listsize];
	
	for(timeii=0;timeii<listsize;timeii++)
	{
		coordinatelist[timeii]=unwrapped[timelist[timeii]];
	}
	
	return coordinatelist;
}



/*--------------------------------------------------------------------------------------*/

/*Method to return distance between position at two timesteps; requires unwrapped coordinates*/
float Trajectory::distance(int first_time, int second_time)const
{
	float distance = (unwrapped[second_time]-unwrapped[first_time]).length();
	return distance;
}

/*--------------------------------------------------------------------------------------*/

/*Method to return distance in xy-plane between position at two timesteps; requires unwrapped coordinates*/
float Trajectory::distance_xy(int first_time, int second_time)const
{
	float distance = (unwrapped[second_time]-unwrapped[first_time]).length_xy();
	return distance;
}

/*--------------------------------------------------------------------------------------*/

/*Method to return distance in xz-plane between position at two timesteps; requires unwrapped coordinates*/
float Trajectory::distance_xz(int first_time, int second_time)const
{
	float distance = (unwrapped[second_time]-unwrapped[first_time]).length_xz();
	return distance;
}

/*--------------------------------------------------------------------------------------*/

/*Method to return distance in yz-plane between position at two timesteps; requires unwrapped coordinates*/
float Trajectory::distance_yz(int first_time, int second_time)const
{
	float distance = (unwrapped[second_time]-unwrapped[first_time]).length_yz();
	return distance;
}

/*--------------------------------------------------------------------------------------*/

/*Method to write trajectory to file*/
void Trajectory::write(string filename)
{
	int timeii;
	
	ofstream output(filename.c_str());
	for(timeii=0;timeii<n_timesteps;timeii++)
	{
		output <<trajectory_ID<<"\t"<< coordinates[timeii].show_x() << "\t" << coordinates[timeii].show_y() << "\t" << coordinates[timeii].show_z() << "\n";
	}
}
