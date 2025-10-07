// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*General class to hold a time correlation function*/
/*Written by David S. Simmons*/

#ifndef TIME_CORRELATION_FUNCTION
#define TIME_CORRELATION_FUNCTION
#include "system.h"
#include "coordinate.h"

namespace std{

class Space_Time_Correlation_Function: public Analysis
{

  protected:
    int n_bins;			//number of spacial bins in correlation function
    int n_times;		//number of times in correlation function
    float bin_size;		//size of bins
    float max_value;		//maximum value included in bins
    
    int * weighting;	//number of measurements contributing to correlation at each time
    
    float * timetable;		//table of times corresponding to correlation data
    
    float** correlation;		//discrete correlation function
    
    float** spatial_inverse;		//stores the spatial inverse transform of the correlation function. Default is that this is simply the spherical Fourier transform
    bool spatial_inverse_calculated;	//tracks whether the spatial inverse correlation function has been computed
    int n_wavenumbers;			//number of wavenumbers if the spatial inverse 
    
    void bin(int timestep, float distance);  //function to place datapoint in bin
    
    void postprocess_list();
    void spherical_fourier();
    
  public:
    ~Space_Time_Correlation_Function();

    Analysis_Type what_are_you(){Analysis_Type type = space_time_correlation_function; return type;};		//virtual method to report the type of analysis
    
    void clear_memory();
    
    Space_Time_Correlation_Function operator + (const Space_Time_Correlation_Function &)const;	//add two correlation functions
    //Space_Time_Correlation_Function operator = (Space_Time_Correlation_Function);
    
    void write(string filename)const;
    void write(ofstream& output)const;
    
    virtual void calculate_spatial_inverse(int n_wavenums);		//virtual method to calculate spatial inverse correlation
    void write_spatial_inverse(string filename)const;
    
    float show(int t, int b)const{return correlation[t][b];};		//return correlation at a given time and bin
    int show_n_bins()const{return n_bins;};				//return number of bins
    float show_bin_size()const{return bin_size;};			//return bin size
};
}

#endif
