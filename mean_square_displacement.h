/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#ifndef MEAN_SQUARE_DISPLACEMENT
#define MEAN_SQUARE_DISPLACEMENT

#include "system.h"
#include <sstream>

namespace std{

class Mean_Square_Displacement: public Analysis
{
    int n_times;
    float * msd;
    int * weighting;
    float * timetable;
    void initialize(System*);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    int atomcount;
    
    /*internal calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    
  public:
    Mean_Square_Displacement();			//default constructor
    Mean_Square_Displacement(const Mean_Square_Displacement &);		//copy constructor
    Mean_Square_Displacement(System*);
    Mean_Square_Displacement operator = (const Mean_Square_Displacement &);	//assignment
    
    void displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj);
    void atomkernel(Trajectory * traj);
    
    void displacementkernel(int timegap,int thisii, int nextii,int species_index, int molecule_index, int atom_type, int atom_index);
    void atomkernel(int species_index, int molecule_index, int atom_type, int atom_index);
    
    float * normalized();
    void write(string);
    void preprocess();
    void postprocess();
    void set(System * sys){initialize(sys);};
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
    
    float show(int t)const{return msd[t];};			//method to return one timestep of msd array
	bool isThreadSafe(){return true;};
};
}

#endif
