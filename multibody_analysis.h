/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Parent class for analysis of multibodies*/
/*Written by David S. Simmons*/

#ifndef MULTIBODY_ANALYSIS
#define MULTIBODY_ANALYSIS

#include <iostream>

#include "multibody_list.h"

namespace std{

class Multibody_Analysis
{
  protected:
    System* system;
  
    Multibody_List * multibody_list;
  
  public:
    Multibody_Analysis();
    Multibody_Analysis(const Multibody_Analysis &);
    Multibody_Analysis operator=(const Multibody_Analysis &);
    ~Multibody_Analysis(){};
    
    
    virtual void analyze(Multibody_List *){cout<<"Error: Multibody list targets not implemented for this analysis method.\n";};
    virtual void list_displacementkernel(int, int, int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Multibody*,int,int,int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void postprocess(){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    
    virtual void write(string)const{cout<<"Error: No standard write method implemented for this multibody_analysis method.\n";};	
};
}

#endif