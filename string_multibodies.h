


#ifndef STRING_MULTIBODIES
#define STRING_MULTIBODIES

#include "multibody.h"
#include "dynamic_cluster_multibodies.h"

namespace std{
  
class String_Multibodies: public Dynamic_Cluster_Multibodies
{
    float ** sigmatrix;		//stores particle sizes
    int n_atomtypes;
    
    float threshold;
        
    void allocate_sig_matrix(string);
    bool clustered_check(Trajectory*, Trajectory*, int, int);
    Coordinate get_imageoffset(Trajectory*, Trajectory*, int, int);
  public:
    
    String_Multibodies();
    String_Multibodies(const String_Multibodies&);
    ~String_Multibodies();
    String_Multibodies operator=(const String_Multibodies&);
    
    String_Multibodies(System * syst, int tgap, float thresh, string sigmatrixname);
  
};
  
  
  
}


#endif