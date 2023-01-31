//Some mathematical tools for arrays of coordinates

#include "trajmath.h"
#include <iostream>

//namespace std{
using namespace std;
/*Function to calculate center of mass of coordinates  when all masses are the same*/
Coordinate std::center_of_mass(const Coordinate* coordinates, int n_coordinates)
{
  int coordinateii;
  Coordinate com;	//create center of mass coordinate
  
  //sum up all coordinates
  for(coordinateii=0;coordinateii<n_coordinates;coordinateii++)
  {
    com+=coordinates[coordinateii];
  }
  com/=n_coordinates;	//divide by number of coordinates
  
  return com;
}

/*-------------------------------------------------------------------------------------------------*/


/*Function to calculate a running value of the center of mass of coordinates when all masses are the same*/
Coordinate* std::center_of_mass_running(const Coordinate* coordinates, int n_coordinates)
{
  int coordinateii;
  Coordinate * com;			//create array of center of mass coordinates
  com = new Coordinate [n_coordinates];	//allocate memory for it
  
  com[0]=coordinates[0];	//first center of mass is just first coordinate
  
  
  //unnormalized center of mass up to each coordinate is just that coordinate plus the previous center of masses
  for(coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]=com[coordinateii-1]+coordinates[coordinateii];
  }
  
  //normalize each center of mass by the number of coordinates summed to obtain it
  for(coordinateii=0;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]/=(coordinateii+1);
  }
  
  return com;
}


/*-------------------------------------------------------------------------------------------------*/


/*Function to calculate a running value of the centroid of a trajectory, approximating the integral over the trajectory by a sum over the available time points using the trapezoidal rule.*/
Coordinate* std::center_of_mass_running(const Coordinate* coordinates, float * times, int n_coordinates)
{
  int coordinateii;
  Coordinate * com;			//create array of center of mass coordinates
  com = new Coordinate [n_coordinates];	//allocate memory for it
  
  com[0]=coordinates[0];	//first center of mass is just first coordinate
  
  com[1]=(coordinates[1]+coordinates[0])*(times[1]-times[0]);
  
  //unnormalized center of mass up to each coordinate is just that coordinate plus the previous center of masses
  for(coordinateii=2;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]=com[coordinateii-1]+(coordinates[coordinateii]+coordinates[coordinateii-1])*(times[coordinateii]-times[coordinateii-1]);
  }
  
  //normalize each center of mass by the time span included in it
  for(coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]/=(2*(times[coordinateii]-times[0]));
  }
  
  return com;
}



/*-------------------------------------------------------------------------------------------------*/


/*Function to calculate a running value of the centroid of a trajectory, approximating the integral over the trajectory by a sum over the available time points using the trapezoidal rule.*/
Coordinate* std::center_of_mass_running(const Coordinate* coordinates, int * times, int n_coordinates)
{
  int coordinateii;
  Coordinate * com;			//create array of center of mass coordinates
  com = new Coordinate [n_coordinates];	//allocate memory for it
  
  com[0]=coordinates[0];	//first center of mass is just first coordinate
  
  com[1]=(coordinates[1]+coordinates[0])*float(times[1]-times[0]);
  
  //unnormalized center of mass up to each coordinate is just that coordinate plus the previous center of masses
  for(coordinateii=2;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]=com[coordinateii-1]+(coordinates[coordinateii]+coordinates[coordinateii-1])*float(times[coordinateii]-times[coordinateii-1]);
  }
  
  //normalize each center of mass by the time span included in it
  for(coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]/=(2*float(times[coordinateii]-times[0]));
  }
  
  return com;
}


/*--------------------------------------------------------------------------------------------*/


//calculate gyration tensor for the overall trajectory. Output is to float gyration_tensor [], which must be of size 6
void std::gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat gyration_tensor)
{
  int coordinateii;
  Coordinate com;
  Coordinate rel_coord;
  com=center_of_mass(coordinates,n_coordinates);	//determine center of mass of coordinates
  
  //initialize elements of gyration tensor
  gyration_tensor[0]=0;
  gyration_tensor[1]=0;
  gyration_tensor[2]=0;
  gyration_tensor[3]=0;
  gyration_tensor[4]=0;
  gyration_tensor[5]=0;
  
  //sum tensor elementsover coordinates
  for(coordinateii=0;coordinateii<n_coordinates;coordinateii++)
  {
    rel_coord = coordinates[coordinateii] - com;	//subtract out the center of mass
    
    gyration_tensor[0]+=rel_coord.show_x()*rel_coord.show_x();
    gyration_tensor[1]+=rel_coord.show_x()*rel_coord.show_y();
    gyration_tensor[2]+=rel_coord.show_x()*rel_coord.show_z();
    gyration_tensor[3]+=rel_coord.show_y()*rel_coord.show_y();
    gyration_tensor[4]+=rel_coord.show_y()*rel_coord.show_z();
    gyration_tensor[5]+=rel_coord.show_z()*rel_coord.show_z();
  }
  
  //normalize gyration tensor by number of coordinates
  gyration_tensor[0]/=n_coordinates;
  gyration_tensor[1]/=n_coordinates;
  gyration_tensor[2]/=n_coordinates;
  gyration_tensor[3]/=n_coordinates;
  gyration_tensor[4]/=n_coordinates;
  gyration_tensor[5]/=n_coordinates;
}


/*--------------------------------------------------------------------------------------------*/


//calculate running value of the gyration tensor for the trajectory. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void std::gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat * gyration_tensor)
{
  int coordinateii, innercoordinateii;
  Coordinate * com;
  com = new Coordinate [n_coordinates];
  Coordinate rel_coord;
  com=center_of_mass_running(coordinates,n_coordinates);	//determine center of mass of coordinates

  //loop over subtrajectory ending in each coordinate
  for(coordinateii=0;coordinateii<n_coordinates;coordinateii++)
  {
    //initialize elements of gyration tensor
    gyration_tensor[coordinateii][0]=0;
    gyration_tensor[coordinateii][1]=0;
    gyration_tensor[coordinateii][2]=0;
    gyration_tensor[coordinateii][3]=0;
    gyration_tensor[coordinateii][4]=0;
    gyration_tensor[coordinateii][5]=0;
    
    //sum tensor over coordinates in trajectory
    for(innercoordinateii=0;innercoordinateii<=coordinateii;innercoordinateii++)
    {
      rel_coord = coordinates[innercoordinateii] - com[coordinateii];	//subtract out center of mass
      
      gyration_tensor[coordinateii][0]+=rel_coord.show_x()*rel_coord.show_x();
      gyration_tensor[coordinateii][1]+=rel_coord.show_x()*rel_coord.show_y();
      gyration_tensor[coordinateii][2]+=rel_coord.show_x()*rel_coord.show_z();
      gyration_tensor[coordinateii][3]+=rel_coord.show_y()*rel_coord.show_y();
      gyration_tensor[coordinateii][4]+=rel_coord.show_y()*rel_coord.show_z();
      gyration_tensor[coordinateii][5]+=rel_coord.show_z()*rel_coord.show_z();
    }
    
    //normalize gyration tensor by number of coordinates up to each endpoint
    gyration_tensor[coordinateii][0]/=(coordinateii+1);
    gyration_tensor[coordinateii][1]/=(coordinateii+1);
    gyration_tensor[coordinateii][2]/=(coordinateii+1);
    gyration_tensor[coordinateii][3]/=(coordinateii+1);
    gyration_tensor[coordinateii][4]/=(coordinateii+1);
    gyration_tensor[coordinateii][5]/=(coordinateii+1);
  }
}


/*------------------------------------------------------------------------------------------------*/


//calculate running value of the gyration tensor for trajectories beginning at the start of the trajectory and ending between first_endcoord and last_endcoord.  Output is to float * gyration_tensor [], which must be a pointer to an array of size (last_endcoord-first_endcoord+1) by 6
void std::gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat * gyration_tensor, int first_endcoord, int last_endcoord)
{
  int coordinateii, innercoordinateii;
  int n_trajectories = last_endcoord-first_endcoord+1;
  Coordinate * com;
  com = new Coordinate [n_trajectories];
  Coordinate rel_coord;
  com=center_of_mass_running(coordinates,n_coordinates);	//determine center of mass of coordinates
  
  //loop over subtrajectory ending in each coordinate
  for(coordinateii=first_endcoord;coordinateii<last_endcoord;coordinateii++)
  {
     //initialize elements of gyration tensor
    gyration_tensor[coordinateii-first_endcoord][0]=0;
    gyration_tensor[coordinateii-first_endcoord][1]=0;
    gyration_tensor[coordinateii-first_endcoord][2]=0;
    gyration_tensor[coordinateii-first_endcoord][3]=0;
    gyration_tensor[coordinateii-first_endcoord][4]=0;
    gyration_tensor[coordinateii-first_endcoord][5]=0;
    
    //sum tensor over coordinates in trajectory
    for(innercoordinateii=0;innercoordinateii<=coordinateii;innercoordinateii++)
    {
      rel_coord = coordinates[innercoordinateii] - com[coordinateii];	//subtract out center of mass
      
      gyration_tensor[coordinateii-first_endcoord][0]+=rel_coord.show_x()*rel_coord.show_x();
      gyration_tensor[coordinateii-first_endcoord][1]+=rel_coord.show_x()*rel_coord.show_y();
      gyration_tensor[coordinateii-first_endcoord][2]+=rel_coord.show_x()*rel_coord.show_z();
      gyration_tensor[coordinateii-first_endcoord][3]+=rel_coord.show_y()*rel_coord.show_y();
      gyration_tensor[coordinateii-first_endcoord][4]+=rel_coord.show_y()*rel_coord.show_z();
      gyration_tensor[coordinateii-first_endcoord][5]+=rel_coord.show_z()*rel_coord.show_z();
    }
    
    //normalize gyration tensor by number of coordinates up to each endpoint
    gyration_tensor[coordinateii-first_endcoord][0]/=(coordinateii+1);
    gyration_tensor[coordinateii-first_endcoord][1]/=(coordinateii+1);
    gyration_tensor[coordinateii-first_endcoord][2]/=(coordinateii+1);
    gyration_tensor[coordinateii-first_endcoord][3]/=(coordinateii+1);
    gyration_tensor[coordinateii-first_endcoord][4]/=(coordinateii+1);
    gyration_tensor[coordinateii-first_endcoord][5]/=(coordinateii+1);
  }
}

//}



/*--------------------------------------------------------------------------------------------*/


//calculate running value of the gyration tensor for a trajectory with uneven time spacing, using trapezoidal rule. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void std::gyration_tensor(const Coordinate * coordinates, float * times, int n_coordinates, sixfloat * gyration_tensor)
{
  int coordinateii, innercoordinateii;
  Coordinate * com;
  com = new Coordinate [n_coordinates];
  Coordinate rel_coord1, rel_coord2;
  com=center_of_mass_running(coordinates,times,n_coordinates);	//determine center of mass of coordinates

  
    gyration_tensor[0][0]=0;
    gyration_tensor[0][1]=0;
    gyration_tensor[0][2]=0;
    gyration_tensor[0][3]=0;
    gyration_tensor[0][4]=0;
    gyration_tensor[0][5]=0;
  
  
  //loop over subtrajectory ending in each coordinate
  for(coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    //initialize elements of gyration tensor
    gyration_tensor[coordinateii][0]=0;
    gyration_tensor[coordinateii][1]=0;
    gyration_tensor[coordinateii][2]=0;
    gyration_tensor[coordinateii][3]=0;
    gyration_tensor[coordinateii][4]=0;
    gyration_tensor[coordinateii][5]=0;
    
    //sum tensor over coordinates in trajectory
    for(innercoordinateii=1;innercoordinateii<=coordinateii;innercoordinateii++)
    {
      rel_coord1 = coordinates[innercoordinateii] - com[coordinateii];	//subtract out center of mass
      rel_coord2 = coordinates[innercoordinateii-1] - com[coordinateii];	//subtract out center of mass
      
      gyration_tensor[coordinateii][0]+=(rel_coord1.show_x()*rel_coord1.show_x()+rel_coord2.show_x()*rel_coord2.show_x())*(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][1]+=(rel_coord1.show_x()*rel_coord1.show_y()+rel_coord2.show_x()*rel_coord2.show_y())*(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][2]+=(rel_coord1.show_x()*rel_coord1.show_z()+rel_coord2.show_x()*rel_coord2.show_z())*(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][3]+=(rel_coord1.show_y()*rel_coord1.show_y()+rel_coord2.show_y()*rel_coord2.show_y())*(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][4]+=(rel_coord1.show_y()*rel_coord1.show_z()+rel_coord2.show_y()*rel_coord2.show_z())*(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][5]+=(rel_coord1.show_z()*rel_coord1.show_z()+rel_coord2.show_z()*rel_coord2.show_z())*(times[innercoordinateii]-times[innercoordinateii-1]);
    }
    
    //normalize gyration tensor by number of coordinates up to each endpoint
    gyration_tensor[coordinateii][0]/=(2*(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][1]/=(2*(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][2]/=(2*(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][3]/=(2*(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][4]/=(2*(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][5]/=(2*(times[coordinateii]-times[0]));
  }
}

/*--------------------------------------------------------------------------------------------*/


//calculate running value of the gyration tensor for a trajectory with uneven time spacing, using trapezoidal rule. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void std::gyration_tensor(const Coordinate * coordinates, int * times, int n_coordinates, sixfloat * gyration_tensor)
{
  int coordinateii, innercoordinateii;
  Coordinate * com;
  com = new Coordinate [n_coordinates];
  Coordinate rel_coord1, rel_coord2;
  com=center_of_mass_running(coordinates,times,n_coordinates);	//determine center of mass of coordinates

  
    gyration_tensor[0][0]=0;
    gyration_tensor[0][1]=0;
    gyration_tensor[0][2]=0;
    gyration_tensor[0][3]=0;
    gyration_tensor[0][4]=0;
    gyration_tensor[0][5]=0;
  
  
  //loop over subtrajectory ending in each coordinate
  for(coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    //initialize elements of gyration tensor
    gyration_tensor[coordinateii][0]=0;
    gyration_tensor[coordinateii][1]=0;
    gyration_tensor[coordinateii][2]=0;
    gyration_tensor[coordinateii][3]=0;
    gyration_tensor[coordinateii][4]=0;
    gyration_tensor[coordinateii][5]=0;
    
    //sum tensor over coordinates in trajectory
    for(innercoordinateii=1;innercoordinateii<=coordinateii;innercoordinateii++)
    {
      rel_coord1 = coordinates[innercoordinateii] - com[coordinateii];	//subtract out center of mass
      rel_coord2 = coordinates[innercoordinateii-1] - com[coordinateii];	//subtract out center of mass
      
      gyration_tensor[coordinateii][0]+=(rel_coord1.show_x()*rel_coord1.show_x()+rel_coord2.show_x()*rel_coord2.show_x())*float(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][1]+=(rel_coord1.show_x()*rel_coord1.show_y()+rel_coord2.show_x()*rel_coord2.show_y())*float(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][2]+=(rel_coord1.show_x()*rel_coord1.show_z()+rel_coord2.show_x()*rel_coord2.show_z())*float(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][3]+=(rel_coord1.show_y()*rel_coord1.show_y()+rel_coord2.show_y()*rel_coord2.show_y())*float(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][4]+=(rel_coord1.show_y()*rel_coord1.show_z()+rel_coord2.show_y()*rel_coord2.show_z())*float(times[innercoordinateii]-times[innercoordinateii-1]);
      gyration_tensor[coordinateii][5]+=(rel_coord1.show_z()*rel_coord1.show_z()+rel_coord2.show_z()*rel_coord2.show_z())*float(times[innercoordinateii]-times[innercoordinateii-1]);
    }
    
    //normalize gyration tensor by number of coordinates up to each endpoint
    gyration_tensor[coordinateii][0]/=(2*float(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][1]/=(2*float(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][2]/=(2*float(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][3]/=(2*float(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][4]/=(2*float(times[coordinateii]-times[0]));
    gyration_tensor[coordinateii][5]/=(2*float(times[coordinateii]-times[0]));
  }
}
