/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Parent class for analysis of multibodies*/
/*Written by David S. Simmons*/


#include "analysis.h"

namespace std{

class Multibody_Analysis
{
  System* system
  
  Multibody_List * multibody_list;
  
  public:
    Multibody_Analysis();
    Multibody_Analysis(const Multibody_Analysis &);
    Multibody_Analysis operator=(const Multibody_Analysis &);
    Multibody_Analysis~(){};
    
    
    virtual void analyze(Multibody_List *){cout<<"Error: Multibody list targets not implemented for this analysis method.\n";};
    virtual void list_displacementkernel(int, int, int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Multibody*,int,int,int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void postprocess(){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
}
}