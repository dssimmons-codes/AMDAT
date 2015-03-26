/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#ifndef MEAN_VELOCITY_UNSTEADY
#define MEAN_VELOCITY_UNSTEADY

#include "system.h"
#include <sstream>
#include "coordinate.h"

namespace std{

class Mean_Velocity_Unsteady: public Analysis
{
    int n_times;
    Coordinate * mean_velocity;
    int * weighting;
    void initialize(System*);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    
    /*internal calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    
  public:
    Mean_Velocity_Unsteady();			//default constructor
    Mean_Velocity_Unsteady(const Mean_Velocity_Unsteady &);		//copy constructor
    Mean_Velocity_Unsteady(System*);
    Mean_Velocity_Unsteady operator = (const Mean_Velocity_Unsteady &);	//assignment
    
    //Analysis_Type what_are_you(){Analysis_Type type = mean_square_displacement; return type;};		//virtual method to report the type of analysis
    
    float * normalized();
    void write(string);
    void set(System * sys){initialize(sys);};
    
    void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int){};
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();
//	bool isThreadSafe(){return true;};
};
}

#endif
