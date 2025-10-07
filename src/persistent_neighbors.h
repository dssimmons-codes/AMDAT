// Copyright (C) 2010-2025 David S. Simmons

/*Persistent_Neighbors class - Identifies particles participating in stringlike cooperative rearrangements and converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#ifndef PERSISTENT_NEIGHBORS
#define PERSISTENT_NEIGHBORS

#include "multibody.h"
#include "dynamic_cluster_multibodies.h"
#include "neighbor_list.h"

namespace std{
  
class Persistent_Neighbors: public Dynamic_Cluster_Multibodies
{
    Neighbor_List * neighbor_list;

        
    bool clustered_check(Trajectory*, Trajectory*, int, int);
    Coordinate get_imageoffset(Trajectory*, Trajectory*, int, int);
  public:
    
    Persistent_Neighbors();
    Persistent_Neighbors(const Persistent_Neighbors&);
    ~Persistent_Neighbors(){};
    Persistent_Neighbors operator=(const Persistent_Neighbors&);
    
    Persistent_Neighbors(System * syst, int tgap, Neighbor_List* nlist);
  
};
  
  
  
}


#endif