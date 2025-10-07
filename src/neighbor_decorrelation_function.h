// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate neighbor decorrelation function*/
/*Written by David S. Simmons*/

#ifndef NEIGHBOR_DECORRELATION_FUNCTION
#define NEIGHBOR_DECORRELATION_FUNCTION

#include "system.h"
#include <sstream>
#include "neighbor_list.h"

namespace std{

class Neighbor_Decorrelation_Function: public Analysis
{
  
    Neighbor_List * n_list;
  
    int n_times;
    float * ndf;
    float * weighting;
    float * timetable;
    void initialize(System*, Neighbor_List*);   
    
    int atomcount;
   
    
    
  public:
    Neighbor_Decorrelation_Function();			//default constructor
    Neighbor_Decorrelation_Function(const Neighbor_Decorrelation_Function &);		//copy constructor
    Neighbor_Decorrelation_Function(System*, Neighbor_List*);
    Neighbor_Decorrelation_Function operator = (const Neighbor_Decorrelation_Function &);	//assignment
    
    Analysis_Type what_are_you(){Analysis_Type type = neighbor_decorrelation_function; return type;};		//virtual method to report the type of analysis
    
    void write(string)const;
    void write(ofstream&)const;
    void set(System * sys, Neighbor_List* nlist){initialize(sys, nlist);};
    
    void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
    
    float show(int t)const{return ndf[t];};			//method to return one timestep of msd array
//	bool isThreadSafe(){return true;};
};
}

#endif
