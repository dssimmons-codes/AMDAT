/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*General class to hold a time correlation function*/
/*Written by David S. Simmons*/

#ifndef TIME_CORRELATION_FUNCTION
#define TIME_CORRELATION_FUNCTION
#include "system.h"
#include "coordinate.h"

namespace std{

class Space-Time_Correlation_Function: public Analysis
{

  protected:
    int n_bins;			//number of spacial bins in correlation function
    int n_times;		//number of times in correlation function
    float bin_size;		//size of bins
    float min_value;		//minimum value included in bins
    float max_value;		//maximum value included in bins
    
    int * weighting;	//number of measurements contributing to correlation at each time
    
    float * timetable;		//table of times corresponding to correlation data
    
    float** correlation;		//unnormalized discrete correlation function
    
    void bin(int timestep, float distance);  //function to place datapoint in bin
    
    void postprocess_list();
    
  public:
    ~Space-Time_Correlation_Function();
    
    
    
    void clear_memory();
    
    
    Space-Time_Correlation_Function operator + (const Space-Time_Correlation_Function &)const;	//add two correlation functions
    //Space-Time_Correlation_Function operator = (Space-Time_Correlation_Function);
    
    void write(string filename)const;
    void write(ofstream& output)const;
    
    float** spherical_fourier(string, int);
    
    
    float show(int t, int b)const{return correlation[t][b];};		//return correlation at a given time and bin
    int show_n_bins()const{return n_bins;};				//return number of bins
    float show_bin_size()const{return bin_size;};			//return bin size
};
}

#endif