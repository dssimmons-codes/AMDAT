/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody class - Class to allow definition of multi-trajectory bodies other than molecules
/*Written by David S. Simmons*/

#include "multibody.h"
#include "system.h"
#include "coordinate.h"

using namespace std;


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
Trajectory Multibody::center_of_mass()const
{
  Trajectory comtraj;
  Coordinate com;
  float mass=0;
  int trajectoryii,timeii;
  int trajlength = system->show_n_timesteps();
  
  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }
  
  /*set trajectory type, length, and mass*/
  comtraj.set(0,trajlength,mass);
  
  /*Determine center of mass of multibody at each time based upon unwrapped coordinates and write to */
  for(timeii=0;timeii<trajlength;timeii++)
  {
    
    /*calculate center of mass at current time*/
    com.set(0,0,0);
    for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    {
      com+=((trajectories[trajectoryii]->show_unwrapped(timeii))*trajectories[trajectoryii]->show_mass());
    }
    com/=mass;
    
    comtraj.set_unwrapped(com,timeii);		//set center of mass trajectory coordinate at current time
  }
  
  /*determine wrapped coordinates from unwrapped coordinates*/
  comtraj.wrap(system->time_dependent_size(),system->time_dependent_boundaries());
  
  return comtraj;
}


/*return centroid trajectory of multibody*/
Trajectory Multibody::centroid()const
{
  Trajectory centraj;
  Coordinate cen;
  float mass=0;
  int trajectoryii,timeii;
  int trajlength = system->show_n_timesteps();
  
  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }
  
  /*set trajectory type, length, and mass*/
  centraj.set(0,trajlength,mass);
  
  /*Determine center of mass of multibody at each time based upon unwrapped coordinates and write to */
  for(timeii=0;timeii<trajlength;timeii++)
  {
    
    /*calculate center of mass at current time*/
    cen.set(0,0,0);
    for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    {
      cen+=((trajectories[trajectoryii]->show_unwrapped(timeii)));
    }
    cen/=n_trajectories;
    
    centraj.set_unwrapped(cen,timeii);		//set center of mass trajectory coordinate at current time
  }
  
  /*determine wrapped coordinates from unwrapped coordinates*/
  centraj.wrap(system->time_dependent_size(),system->time_dependent_boundaries());
  
  return centraj;
}