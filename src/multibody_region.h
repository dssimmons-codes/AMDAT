// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_Region class - stores list of multibody objects that are within a given region*/
/*Written by David S. Simmons*/


#ifndef MULTIBODY_REGION
#define MULTIBODY_REGION

#include "coordinate.h"
#include "multibody_list.h"
#include "multibody_analysis.h"

namespace std{
  
class Multibody_Region: public Multibody_List, public Multibody_Analysis   
{
  private:
    
    Coordinate upperbound;
    Coordinate lowerbound;   

  public:
    

    
    Multibody_Region();
    Multibody_Region(System*syst);
    Multibody_Region(System*syst, Coordinate low, Coordinate high);
    Multibody_Region(const Multibody_Region &);
    Multibody_Region operator=(const Multibody_Region &);
    
    void analyze(Multibody_List * mblist);
    void list_displacementkerneyl(int,int,int){};
    void listkernel(Multibody*,int,int,int);
    void postprocess_list();
    
    void write(string) const;
  
};
  
}


#endif
