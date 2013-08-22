//Header file for a number of mathematical tools for arrays of coordinates

//problem seems to be that function declarations in .h and .cpp files are somehow being treated as distinct

#include "coordinate.h"
#ifndef TRAJMATH
#define TRAJMATH

namespace std {
  
typedef float sixfloat [6];
 
Coordinate center_of_mass(const Coordinate* coordinates, int n_coordinates);
  
Coordinate * center_of_mass_running(const Coordinate* coordinates, int n_coordinates);
 
/*Function to calculate a running value of the centroid of a trajectory, approximating the integral over the trajectory by a sum over the available time points using the trapezoidal rule.*/
Coordinate* center_of_mass_running(const Coordinate* coordinates, float * times, int n_coordinates);

/*Function to calculate a running value of the centroid of a trajectory, approximating the integral over the trajectory by a sum over the available time points using the trapezoidal rule.*/
Coordinate* center_of_mass_running(const Coordinate* coordinates, int * times, int n_coordinates);

//calculate gyration tensor for the overall trajectory. Output is to float gyration_tensor [], which must be of size 6
void gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat gyration_tensor);

//calculate running value of the gyration tensor for the trajectory. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat * gyration_tensor);

//calculate running value of the gyration tensor for trajectories beginning at the start of the trajectory and ending between first_endcoord and last_endcoord.  Output is to float * gyration_tensor [], which must be a pointer to an array of size (last_endcoord-first_endcoord+1) by 6
void gyration_tensor(const Coordinate * coordinates, int n_coordinates, sixfloat * gyration_tensor, int first_endcoord, int last_endcoord);

//calculate running value of the gyration tensor for a trajectory with uneven time spacing, using trapezoidal rule. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void gyration_tensor(const Coordinate * coordinates, float * times, int n_coordinates, sixfloat * gyration_tensor);

//calculate running value of the gyration tensor for a trajectory with uneven time spacing, using trapezoidal rule. Output is to float * gyration_tensor [], which must be a pointer to an array of size n_coordinates by 6
void gyration_tensor(const Coordinate * coordinates, int * times, int n_coordinates, sixfloat * gyration_tensor);


}

#endif