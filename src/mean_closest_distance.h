// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate radial_distribution_function from n^2 method*/
/*Written by David S. Simmons*/

#ifndef MEAN_CLOSEST_DISTANCE
#define MEAN_CLOSEST_DISTANCE

#include "system.h"
#include <sstream>
#include "analysis_onetime.h"

namespace std{

class Mean_Closest_Distance: public Analysis_Onetime
{

    int n_times;
    float * time_mcd;
    float mean_mcd;
    int * n_atoms_i;
    int * n_atoms_j;
    
    int total_center_atoms;
    
    float tempdistance;  //Tracks current shortest distance found
    
  public:
    Mean_Closest_Distance();			//default constructor    
    Mean_Closest_Distance(const Mean_Closest_Distance &);		//copy constructor
    Mean_Closest_Distance(System*sys, int timescheme);
    
    Mean_Closest_Distance operator = (const Mean_Closest_Distance &);	//assignment
    
    //Mean_Closest_Distance operator+ (const Mean_Closest_Distance &);
    
    void set(System*sys, int timescheme);
    
    Analysis_Type what_are_you(){Analysis_Type type = radial_distribution_function; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list();
    void write(string);
    void write(ofstream& output)const;
    

    
    
//	bool isThreadSafe(){return true;};
};
}

#endif
