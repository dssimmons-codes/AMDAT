// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Method for parent class for analysis of multibodies*/
/*Written by David S. Simmons*/

#include "multibody_analysis.h"

using namespace std;

Multibody_Analysis::Multibody_Analysis()
{
  system = 0;
  multibody_list = 0;
}

Multibody_Analysis::Multibody_Analysis(const Multibody_Analysis & copy)
{
  system = copy.system;
  multibody_list = copy.multibody_list;
}


Multibody_Analysis Multibody_Analysis::operator=(const Multibody_Analysis & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    multibody_list = copy.multibody_list;
  }
  
  return *this;
}