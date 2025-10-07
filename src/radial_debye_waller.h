// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Radial Debye Waller class*/
/*Written by David S. Simmons*/

#ifndef RADIAL_DEBYE_WALLER
#define RADIAL_DEBYE_WALLER

#include "coordinate.h"
#include "analysis.h"
#include "system.h"

namespace std{

class Radial_Debye_Waller: public Analysis
{
    int time_index;		//index of timegap at which radial debye waller factor will be calculated
    int time_weighting;		//number of independent configuration spacing corresponding to this timegap
    float * msd;
    float * n_atoms;
    float * density;
    int n_bins;
    float bin_size;
    Coordinate center;
  
    int currenttime, nexttime;
    
  public: 
    Radial_Debye_Waller(System* sys, int timeii, int bincount, float maxrange, Coordinate c);
   
    Analysis_Type what_are_you(){Analysis_Type type = radial_debye_waller; return type;};		//virtual method to report the type of analysis
    
    void write(string);
    void write(ofstream&)const;
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void postprocess_list();
};
}

#endif