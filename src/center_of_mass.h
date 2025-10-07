// Copyright (C) 2010-2025 David S. Simmons

/*Calculates center of mass of a set of coordinates*/

#include "coordinate.h"

namespace std
{
  
Coordinate center_of_mass(Coordinate* coordinates, int n_coordinates);
  
Coordinate * center_of_mass_running(Coordinate* coordinates, int n_coordinates);
  
}