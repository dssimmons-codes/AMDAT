// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate velocity autocorrelation function*/
/*Written by David S. Simmons*/

#ifndef VELOCITY_AUTOCORRELATION
#define VELOCITY_AUTOCORRELATION

#include <sstream>
#include "mean_square_displacement.h"

namespace std{

class Velocity_Autocorrelation
{
    int n_times;			//number of times in array
    float * vac;			//velocity autocorrelation function data
    double * fourier_space_vac_real;
    double * fourier_space_vac_imag;
    float time_unit;
    int fouriersize;			//size of array from fourier transform
    
    Mean_Square_Displacement * msd;	//mean square displacement object
    System * system;			//system for which vac is being calculated
    
    float * timetable;

    void calculate();
    
    bool linearity_error;		// track whether msd data has been found to be nonlinear
    

  public:
    Velocity_Autocorrelation(){};
    Velocity_Autocorrelation(Mean_Square_Displacement*);		//constuctor
    Velocity_Autocorrelation(Mean_Square_Displacement*,int);		//constuctor to include a timestep to end calculation at. This can be important if msd data consists of a number of linearly spaced time steps followed by a different time stepping (as in the case of linear blocks)
    void initialize(Mean_Square_Displacement*);
    void initialize(Mean_Square_Displacement*,int);
    void fourier_transform();
    void write(string)const;
    void write(ofstream&)const;
    void write_fourier(string)const;

};

}
#endif