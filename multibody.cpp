/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody class - Class to allow definition of multi-trajectory bodies other than molecules*/
/*Written by David S. Simmons*/


#include "multibody.h"
#include "system.h"
#include "coordinate.h"


using namespace std;

/*need to set trajectory_id on creation*/


/*Default constructor*/
Multibody::Multibody()
{
  n_trajectories = 0;
  trajectories = new Trajectory* [n_trajectories];
}



/*Copy constructor*/
Multibody::Multibody(const Multibody & copy)
{
  int trajectoryii;
  
  n_trajectories = copy.n_trajectories;
  trajectories = new Trajectory* [n_trajectories];
  
  /*Copy over pointers to trajectories in multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    trajectories[trajectoryii]=copy.trajectories[trajectoryii];
  }  
}



/*Destructor*/
Multibody::~Multibody()
{
  delete [] trajectories;
}



/*method to return pointer to specified trajectory in multibody*/
Trajectory * Multibody::operator()(int bodyindex)
{
  return trajectories[bodyindex];
}



/*equality operator*/
Multibody Multibody::operator=(const Multibody & copy)
{
  int trajectoryii;
  
  if(this!=&copy)
  {
    n_trajectories = copy.n_trajectories;
    trajectories = new Trajectory* [n_trajectories];
  
    /*Copy over pointers to trajectories in multibody*/
    for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    {
      trajectories[trajectoryii]=copy.trajectories[trajectoryii];
    }  
  }
  
  return *this;
}




/*construct multibody with n_bodies initially corresponding to null pointers*/
Multibody::Multibody(int n_bodies)
{
  int trajectoryii;
  
  n_trajectories=n_bodies;
  trajectories = new Trajectory* [n_trajectories];
  
  /*Copy over pointers to trajectories in multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    trajectories[trajectoryii]=0;
  }  
}



/*Construct multibody based on an array of trajectories*/
Multibody::Multibody(int n_bodies, Trajectory ** bodies)
{
  int trajectoryii;
  
  n_trajectories=n_bodies;
  trajectories = new Trajectory* [n_trajectories];
  
  /*Copy over pointers to trajectories in multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    trajectories[trajectoryii]=bodies[trajectoryii];
  }  
}


void Multibody::set(System*sys)
{
  system=sys;
}


/*reset multibody based on an array of trajectories*/
void Multibody::set(int n_bodies, Trajectory ** bodies)
{
  int trajectoryii;
  delete [] trajectories;
  
  n_trajectories=n_bodies;
  trajectories = new Trajectory* [n_trajectories];
  
  /*Copy over pointers to trajectories in multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    trajectories[trajectoryii]=bodies[trajectoryii];
  }  
  
}



/*return center of mass trajectory of multibody*/
void Multibody::center_of_mass_trajectory()
{
  Coordinate com;
  mass=0;
  int trajectoryii,timeii;
  n_timesteps = system->show_n_timesteps();
  
  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }
  
  type=0;
  
  /*Determine center of mass of multibody at each time based upon unwrapped coordinates and write to */
  for(timeii=0;timeii<n_timesteps;timeii++)
  {
    
    /*calculate center of mass at current time*/
    com.set(0,0,0);
    for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    {
      com+=((trajectories[trajectoryii]->show_unwrapped(timeii))*trajectories[trajectoryii]->show_mass());
    }
    com/=mass;
    
    set_unwrapped(com,timeii);		//set center of mass trajectory coordinate at current time
  }
  
  /*determine wrapped coordinates from unwrapped coordinates*/
  wrap(system->time_dependent_size(),system->time_dependent_boundaries());
}


/*return centroid trajectory of multibody*/
void Multibody::centroid_trajectory()
{
  Coordinate cen;
  mass=0;
  int trajectoryii,timeii;
  n_timesteps = system->show_n_timesteps();
  
  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }
  
  type=0;
  
  /*Determine centroid of multibody at each time based upon unwrapped coordinates and write to */
  for(timeii=0;timeii<n_timesteps;timeii++)
  {
    
    /*calculate centroid at current time*/
    cen=calculate_centroid(timeii);
    set_unwrapped(cen,timeii);		//set center of mass trajectory coordinate at current time
  }
  
  /*determine wrapped coordinates from unwrapped coordinates*/
  wrap(system->time_dependent_size(),system->time_dependent_boundaries());
 
}


Coordinate Multibody::calculate_centroid(int timeii)const
{
  int bodyii;
  Coordinate centroid(0,0,0);
  
  for(bodyii=0;bodyii<n_trajectories;bodyii++)
  {
    centroid += (trajectories[bodyii]->show_unwrapped(timeii));
  }
  centroid/=n_trajectories;
  
  return centroid;
}


/*Method to calculate multibody gyration radis at a given time*/
float Multibody::square_gyration_radius(int timeii)
{
  int bodyii;
  Coordinate centroid(0,0,0);
  float gyration_radius=0;
  
  centroid=calculate_centroid(timeii);
  
  for(bodyii=0;bodyii<n_trajectories;bodyii++)
  {
    gyration_radius+=(trajectories[bodyii]->show_unwrapped(timeii)-centroid).length_sq();
  }
  gyration_radius/=n_trajectories;
  
  return gyration_radius;
  
}
