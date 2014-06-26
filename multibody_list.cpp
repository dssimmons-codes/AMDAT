/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_List class*/
/*Written by David S. Simmons*/

#include <fstream>
#include <stdlib.h>
#include "multibody_list.h"
#include "analysis.h"
#include "multibody_analysis.h"
#include "system.h"

using namespace std;

Multibody_List::Multibody_List();
{
  sys=0;
  
  n_times=0;
  n_bodies=-2;
  
  time_conversion.resize(n_times);
  multibodies.resize(n_times);
  
}





Multibody_List::Multibody_List(const Multibody_List & copy);
{
  sys=copy.sys;
  
  n_times=copy.n_times;
  n_bodies=copy.n_bodies;
  
  time_conversion = copy.time_conversion;
  multibodies=copy.multibodies;
  
}







Multibody_List::Multibody_List(int timecount);
{
  int timeii;
  
  sys=copy.sys;
  n_times=0;
  n_bodies=-2;
  
  time_conversion.resize(n_times);
  multibodies.resize(n_times);
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    multibodies[timeii].resize(0);
  }
}