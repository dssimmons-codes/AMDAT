/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#ifndef RADIAL_DISRIBUTION_FUNCTION
#define RADIAL_DISRIBUTION_FUNCTION

#include "system.h"
#include <sstream>

namespace std{

class Radial_Distribution_Function: public Analysis
{
    float max_distance;
    int n_bins;
    float bin_size;
    int n_times;
    float ** time_rdf;
    float * mean_rdf;
    int * n_atoms_i;
    int * n_atoms_j;
  
    void initialize(System*);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    
    /*internal calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    
  public:
    Radial_Distribution_Function();			//default constructor
    Radial_Distribution_Function(const Radial_Distribution_Function &);		//copy constructor
    Radial_Distribution_Function(System*sys, int n_bins, float max_distance=0);
    
    Radial_Distribution_Function operator = (const Radial_Distribution_Function &);	//assignment
    
    //Radial_Distribution_Function operator+ (const Radial_Distribution_Function &);
    
    //Analysis_Type what_are_you(){Analysis_Type type = mean_square_displacement; return type;};		//virtual method to report the type of analysis
    
    float * normalized();
    void write(string);
    void set(System * sys){initialize(sys);};
    
    void analyze(Trajectory_List *,Trajectory_List *); //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
    
    float show(int t)const{return msd[t];};			//method to return one timestep of msd array
//	bool isThreadSafe(){return true;};
};
}

#endif
