/*Calculates center of mass of a set of coordinates*/

#include "center_of_mass.h"

using namespace std;

/*Function to calculate center of mass of coordinates when all masses are the same*/
Coordinate center_of_mass(Coordinate* coordinates, int n_coordinates)
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

/*--------------------------------------------------------------------------------*/


/*Function to calculate a running value of the center of mass of coordinates when all masses are the same*/
Coordinate* center_of_mass_running(Coordinate* coordinates, int n_coordinates)
{

  Coordinate * com;			//create array of center of mass coordinates
  com = new Coordinate [n_coordinates];	//allocate memory for it
  int coordinateii=0;
  com[0]=coordinates[coordinateii];	//first center of mass is just first coordinate


  //unnormalized center of mass up to each coordinate is just that coordinate plus the previous center of masses
  for(int coordinateii=1;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii].set(0,0,0);
    com[coordinateii]=com[coordinateii-1]+coordinates[coordinateii];
  }

  //normalize each center of mass by the number of coordinates summed to obtain it
  for(int coordinateii=0;coordinateii<n_coordinates;coordinateii++)
  {
    com[coordinateii]/=(coordinateii+1);
  }

  return com;
}
