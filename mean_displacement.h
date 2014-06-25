/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate mean-displacement*/
/*Written by David S. Simmons*/

#ifndef MEAN_DISPLACEMENT
#define MEAN_DISPLACEMENT

#include "system.h"
#include "coordinate.h"
#include <sstream>

namespace std{

class Mean_Displacement: public Analysis
{
    int n_times;
    Coordinate * md;
    int * weighting;
    float * timetable;
    void initialize(System*);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    int atomcount;
    
    /*internal calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    
  public:
    Mean_Displacement();			//default constructor
    Mean_Displacement(const Mean_Displacement &);		//copy constructor
    Mean_Displacement(System*);
    Mean_Displacement operator = (const Mean_Displacement &);	//assignment
    
    Analysis_Type what_are_you(){Analysis_Type type = mean_displacement; return type;};		//virtual method to report the type of analysis
    
    float * normalized();
    void write(string);

    void set(System * sys){initialize(sys);};
    
    void preprocess();
     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
    
    Coordinate show(int t)const{return md[t];};			//method to return one timestep of msd array
};
}

#endif