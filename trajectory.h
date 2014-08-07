/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Parent class for classes that store trajectory information, such as atoms and molecules*/
/*Written by David S. Simmons*/

#ifndef TRAJECTORY
#define TRAJECTORY

#include "coordinate.h"
#include <string>

namespace std {

enum trajectory_type {general,atom,molecule};

/*-------------------------*/
/*--Trajectory Class--*/
/*-------------------------*/

class Trajectory
{    
  protected:
    trajectory_type trajtype;	//enum variable to store type of trajectory for later pointer typecasting
    
    bool is_unwrapped;		//boolean variable to determine whether unwrapped coordinates are defined
    bool is_wrapped; 		//boolean variable to determine whether wrapped coordinates are defined
    bool is_velocity;
    
    int n_timesteps;		//number of timesteps in trajectories (same for all trajectories)

    int type;			//type of trajectory
    
    int trajectory_ID;		//unique ID of trajectory
    float mass;			//mass associated with trajectory

    Coordinate *coordinates;		//array of wrapped coordinates
    Coordinate *unwrapped;		//array to hold unwrapped coordinates
    Coordinate *velocity;		//array of velocities

    void clear_memory();					//method to clear memory allocated to arrays
    Coordinate unwrap_displacement(int,int,const Coordinate&);		//method to return atom displacement vector between two selected timesteps

  public:
    Trajectory(int timesteps=0);	//constructor to instantiate coordinate array without defining it
    Trajectory(int,Coordinate *);	//constructor to instantiate coordinate array with full definition
   ~Trajectory();			//destructor
    
    virtual void set(int,int,float m=1);		//change number of timesteps and reallocate memory accordingly
    virtual void set(int,int,Coordinate*,int);	//method to fully define object, including coordinate list
    
    void set(const Coordinate &, int);		//method to set single coordinate in coordinate list
    void set_unwrapped(const Coordinate &, int);
    void set_velocity(const Coordinate &, int);
    void set_mass(float m){mass = m;};	//method to set mass of atom
    void set_type(int t){type=t;};		//set type
    
    void unwrap(const Coordinate &);
    void wrap (const Coordinate *,Coordinate **);
    
    float distance(int,int)const;		//calculate spatial distance between two timesteps
    float distance_xy(int,int)const;		//calculate spatial distance between two timesteps in xy plane
    float distance_xz(int,int)const;		//calculate spatial distance between two timesteps in xz plane
    float distance_yz(int,int)const;		//calculate spatial distance between two timesteps in yz plane
    
    Coordinate show_unwrapped(int T)const{return unwrapped[T];};			
    Coordinate * show_unwrapped()const{return unwrapped;};					//method to return array of displacement
    Coordinate show_coordinate(int T)const{return coordinates[T];};				//method to return single coordinate in array		
    Coordinate * show_coordinates()const{return coordinates;};
    Coordinate * show_coordinates(int*timelist, int listsize)const;		//returns list of coordinates at requested times
    const Coordinate & show_velocity(int T)const{return velocity[T];}			//method to return velocity

    
    float show_mass()const{return mass;};				//method to show mass
    bool check_unwrapped()const{return is_unwrapped;};			//check whether unwrapped coordinates are present
    bool check_wrapped()const {return is_wrapped;};			//check whether wrapped coordinates are present
    void set_wrapped_state(bool val){is_wrapped=val;};			//set whether wrapped coordinates exist
    void set_unwrapped_state(bool val){is_unwrapped=val;};		//set whether unwrapped coordinates exist
    void set_velocity_state(bool val){is_velocity=val;};		//set whether velocity data exists
    int show_type()const{return type;};					//show trajectory type
    
    
    int show_trajectory_ID()const{return trajectory_ID;};
    void set_trajectory_ID(int ID){trajectory_ID = ID;};

    void write(string filename);				//write trajectory to file
	};

}

#endif