/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Parent class for analysis of multibodies*/
/*Written by David S. Simmons*/


#include "analysis.h"

namespace std{

class Multibody_Analysis:public Analysis
{
  Multibody_List * multibody_list;
  
  int n_bodies;	//number of bodies in each multibody
  
  public:
    
    virtual void analyze(Multibody_List *){cout<<"Error: Multibody list targets not implemented for this analysis method.\n";};
    virtual void list_displacementkernel(int, int, int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Multibody*){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void postprocess_multibody_list(){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
}
}