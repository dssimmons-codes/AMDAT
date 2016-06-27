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
  //trajectories = new Trajectory* [n_trajectories];
}

Multibody::Multibody(System * sys)
{
  system=sys;
}



/*Copy constructor*/
Multibody::Multibody(const Multibody & copy)
//:Trajectory(copy)
{
  int trajectoryii;

  trajectories=copy.trajectories;
  tref=copy.tref;
  relative_image_index=copy.relative_image_index;
  //trajectories = new Trajectory* [n_trajectories];

  /*Copy over pointers to trajectories in multibody*/
  //for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  //{
  //  trajectories[trajectoryii]=copy.trajectories[trajectoryii];
  //}
}


#ifdef NEVER
/*Destructor*/
Multibody::~Multibody()
{
}
#endif


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
    //Trajectory::operator=(copy);
    trajectories=copy.trajectories;
    tref=copy.tref;
    relative_image_index=copy.relative_image_index;
    //trajectories = new Trajectory* [n_trajectories];

    /*Copy over pointers to trajectories in multibody*/
    //for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    //{
    //  trajectories[trajectoryii]=copy.trajectories[trajectoryii];
    //}
  }

  return *this;
}




/*construct multibody with n_bodies initially corresponding to null pointers*/
Multibody::Multibody(int n_bodies)
{
  int trajectoryii;

  int n_trajectories=n_bodies;
  //trajectories.resize(n_trajectories);
  tref=-1;
  //trajectories = new Trajectory* [n_trajectories];

  /*Copy over pointers to trajectories in multibody*/
  //for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  //{
   // trajectories[trajectoryii]=0;
  //}
}



/*Construct multibody based on an array of trajectories*/
Multibody::Multibody(int n_bodies, Trajectory ** bodies)
{
  int trajectoryii;

  int n_trajectories=n_bodies;
  //trajectories.resize(n_trajectories);
  //trajectories = new Trajectory* [n_trajectories];

  /*Copy over pointers to trajectories in multibody*/
  //for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  //{
 //   trajectories[trajectoryii]=bodies[trajectoryii];
 // }
}


void Multibody::set(System*sys)
{
  system=sys;
}


/*reset multibody based on an array of trajectories*/
void Multibody::set(int n_bodies, Trajectory ** bodies)
{
  int trajectoryii;
  trajectories.clear();

  int n_trajectories=n_bodies;
  trajectories.resize(n_trajectories);

  /*Copy over pointers to trajectories in multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    trajectories[trajectoryii]=bodies[trajectoryii];
  }

}



/*return center of mass trajectory of multibody*/
void Multibody::center_of_mass_trajectory(int trajtype)
{
  Coordinate com;
  mass=0;
  int trajectoryii,timeii;
  n_timesteps = system->show_n_timesteps();
  int n_trajectories=trajectories.size();

  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }

  type=trajtype;

  /*Determine center of mass of multibody at each time based upon unwrapped coordinates and write to */
  for(timeii=0;timeii<n_timesteps;timeii++)
  {

    /*calculate center of mass at current time*/
    com.set(0,0,0);
    for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
    {
      //com+=((trajectories[trajectoryii]->show_unwrapped(timeii))*trajectories[trajectoryii]->show_mass());
      com+=consistent_position(trajectoryii,timeii)*trajectories[trajectoryii]->show_mass();
    }
    com/=mass;

    set_unwrapped(com,timeii);		//set center of mass trajectory coordinate at current time
  }

  /*determine wrapped coordinates from unwrapped coordinates*/
  wrap(system->time_dependent_size(),system->time_dependent_boundaries());
}


/*return centroid trajectory of multibody*/
void Multibody::centroid_trajectory(int trajtype)
{
  Coordinate cen;
  mass=0;
  int trajectoryii,timeii;
  n_timesteps = system->show_n_timesteps();
  int n_trajectories=trajectories.size();

  /*Calculate total mass of multibody*/
  for(trajectoryii=0;trajectoryii<n_trajectories;trajectoryii++)
  {
    mass+=trajectories[trajectoryii]->show_mass();
  }

  type=trajtype;

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
  int n_trajectories=trajectories.size();
  for(bodyii=0;bodyii<n_trajectories;bodyii++)
  {
    //centroid += (trajectories[bodyii]->show_unwrapped(timeii));
    centroid+=consistent_position(bodyii,timeii);
  }
  centroid/=float(n_trajectories);
  //centroid=trajectories[bodyii-1]->show_unwrapped(timeii);

  return centroid;
}


/*Method to calculate multibody gyration radis at a given time*/
float Multibody::square_gyration_radius(int timeii)
{
  int bodyii;
  Coordinate centroid(0,0,0);
  float gyration_radius=0;
  int n_trajectories=trajectories.size();

  centroid=calculate_centroid(timeii);

  for(bodyii=0;bodyii<n_trajectories;bodyii++)
  {
    //gyration_radius+=(trajectories[bodyii]->show_unwrapped(timeii)-centroid).length_sq();
    gyration_radius+=(consistent_position(bodyii,timeii)-centroid).length_sq();
  }
  gyration_radius/=float(n_trajectories);

  return gyration_radius;

}

void Multibody::gyr_tensor(int timeii, sixfloat* g_tensor)
{
    Coordinate * coordinates;
    int trajii;
    int n_trajectories=trajectories.size();

    coordinates = new Coordinate [n_trajectories];

    for(trajii=0;trajii<n_trajectories;trajii++)
    {
        //coordinates[trajii]=trajectories[trajii]->show_unwrapped(timeii);
	coordinates[trajii]=consistent_position(trajii,timeii);
    }

    gyration_tensor(coordinates,  n_trajectories, g_tensor);

}

#ifdef NEVER
threefloat Multibody::principle_axes(int timeii)
{
  TNT::Array1D<float> eigvals;
  TNT::Array2D<float> rgarray(3,3,0.0);
  threefloat axes;

  sixfloat rgtensor;
  rgeensor = gyrtensor(timeii);

  rgarray[0][0]=rgtensor[blockii][timeii][0];
  rgarray[1][0]=rgtensor[blockii][timeii][1];
  rgarray[2][0]=rgtensor[blockii][timeii][2];
  rgarray[0][1]=rgtensor[blockii][timeii][1];
  rgarray[1][1]=rgtensor[blockii][timeii][3];
  rgarray[2][1]=rgtensor[blockii][timeii][4];
  rgarray[0][2]=rgtensor[blockii][timeii][2];
  rgarray[1][2]=rgtensor[blockii][timeii][4];
  rgarray[2][2]=rgtensor[blockii][timeii][5];

  JAMA::Eigenvalue <float> eigmat(rgarray);

  eigmat.getRealEigenvalues(eigvals);

  axes[0]=pow(eigvals[0],0.5);
  axes[1]=pow(eigvals[1],0.5);
  axes[2]=pow(eigvals[2],0.5);

  return axes;
}
#endif


bool Multibody::trajectory_check(Trajectory* check)
{
  
  for(int bodyii=0;bodyii<trajectories.size();bodyii++)
  {
    if(trajectories[bodyii]==check)
    {
      return true;
    }
  }
  return false;
}



Coordinate Multibody::consistent_position(int trajii, int timeii)const
{
  return ((trajectories[trajii])->show_coordinate(tref))+((system->size(tref))*relative_image_index[trajii])+(trajectories[trajii])->displacement_vector(tref,timeii);
}

void Multibody::add_body(Trajectory* new_trajectory)	//add body, assuming that relative image index is 0 and tref is 0 (usually will not produce robustly useful results)
{
  Coordinate def(0,0,0);
  trajectories.push_back(new_trajectory);
  relative_image_index.push_back(def);
}

void Multibody::add_body(Trajectory* new_trajectory, Coordinate imageoffset)
{
  trajectories.push_back(new_trajectory);
  relative_image_index.push_back(imageoffset);
  
}


void Multibody::absorb_multibody(const Multibody& target)
{
  int bodyii;
  
  for(bodyii=0;bodyii<target.trajectories.size();bodyii++)	//loop over bodies in other multibody
  {
    if(!trajectory_check(target.trajectories[bodyii]))	//only proceed if the body being absorbed is not already in this multibody
    {
      add_body(target.trajectories[bodyii]);	//add body in other multibody to this one
    }
  }
}