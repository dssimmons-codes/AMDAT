// Copyright (C) 2010-2025 David S. Simmons

/*Relative_Displacement_Strings class - Identifies particles remaining in spacial promity to one anotherand converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#ifndef RELATIVE_DISPLACEMENT_STRINGS
#define RELATIVE_DISPLACEMENT_STRINGS

#include "multibody.h"
#include "dynamic_cluster_multibodies.h"
#include "neighbor_list.h"

namespace std{
  
class Relative_Displacement_Strings: public Dynamic_Cluster_Multibodies
{
    Neighbor_List * neighbor_list;
    
    float threshold;
    int steps_for_averaging;
        
    bool clustered_check(Trajectory*, Trajectory*, int, int);
    Coordinate get_imageoffset(Trajectory*, Trajectory*, int, int);
  public:
    
    Relative_Displacement_Strings();
    Relative_Displacement_Strings(const Relative_Displacement_Strings&);
    Relative_Displacement_Strings operator=(const Relative_Displacement_Strings&);
    
    Relative_Displacement_Strings(System * syst, int tgap, Neighbor_List* nlist, float thresh, int avgsteps);
  
};
  
  
  
}


#endif