// Copyright (C) 2010-2025 David S. Simmons

/*Comover_Multibodies class - Identifies particles remaining in spacial promity to one anotherand converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#ifndef COMOVER_MULTIBODIES
#define COMOVER_MULTIBODIES

#include "multibody.h"
#include "dynamic_cluster_multibodies.h"

namespace std{
  
class Comover_Multibodies: public Dynamic_Cluster_Multibodies
{
    float ** sigmatrix;		//stores particle sizes
    int n_atomtypes;
    
    float threshold;
        
    void allocate_sig_matrix(string);
    bool clustered_check(Trajectory*, Trajectory*, int, int);
    Coordinate get_imageoffset(Trajectory*, Trajectory*, int, int);
  public:
    
    Comover_Multibodies();
    Comover_Multibodies(const Comover_Multibodies&);
    ~Comover_Multibodies();
    Comover_Multibodies operator=(const Comover_Multibodies&);
    
    Comover_Multibodies(System * syst, int tgap, float thresh, string sigmatrixname);
  
};
  
  
  
}


#endif