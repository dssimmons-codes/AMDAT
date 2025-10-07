/*Molecular Dynamics Analysis Toolkit*/
/*Class to calculate non-gaussian parameter*/
/*Written by David S. Simmons*/

#ifndef NON_GAUSSIAN_PARAMETER
#define NON_GAUSSIAN_PARAMETER

#include "system.h"
#include <string>
#include "mean_square_displacement.h"

namespace std{

class Non_Gaussian_Parameter: public Analysis
{
    int n_times;
    float * ngp;
    long int * weighting;
    Mean_Square_Displacement const * msd;
    float * timetable;
    int atomcount;
    float * n_atoms;

    //calculation variables
    int currenttime, nexttime, currenttimegap;
    
  public:
    Non_Gaussian_Parameter();
    Non_Gaussian_Parameter(const Non_Gaussian_Parameter &);
    Non_Gaussian_Parameter(System*, const Mean_Square_Displacement*);
    Non_Gaussian_Parameter operator=(const Non_Gaussian_Parameter &);
    
    Analysis_Type what_are_you(){Analysis_Type type = non_gaussian_parameter; return type;};		//virtual method to report the type of analysis
    
    void write(string)const;
    void write(ofstream&)const;
    int max()const;

     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void listkernel(Trajectory* , int ,int , int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii);
    void postprocess_bins();
    //bool isThreadSafe(){return true;};
};

}

#endif