/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*General class to hold a time correlation function*/
/*Written by David S. Simmons*/

#ifndef TIME_CORRELATION_FUNCTION
#define TIME_CORRELATION_FUNCTION
#include "system.h"
#include "coordinate.h"

namespace std{

class Time_Correlation_Function: public Analysis
{
  private:
    void displacementkernel(int,int, int, int, int, int, int){};
    void atomkernel(int, int, int, int){};
    void cellkernel(int, Atom_Trajectory, Atom_Trajectory){};
  protected:
    int n_bins;			//number of spacial bins in correlation function
    int n_times;		//number of times in correlation function
    float bin_size;		//size of bins
    float min_value;		//minimum value included in bins
    float max_value;		//maximum value included in bins
    int n_atoms_looped;		//counts how many atoms have been looped over in generating the correlation function
    double * n_atoms;				//more generally holds number of atoms for normalization at each time, when this atom count is different at different times.
    int * weighting;		//number of displacement measurements at each time
    float * timetable;		//table of times corresponding to correlation data
    
    int** correlation;		//unnormalized discrete correlation function
    
    void bin(int timestep, float distance);  //function to place datapoint in bin
    
  public:
    ~Time_Correlation_Function();
    void clear_memory();
    Time_Correlation_Function operator + (const Time_Correlation_Function &)const;	//add two correlation functions
    //Time_Correlation_Function operator = (Time_Correlation_Function);
    float ** normalized()const;		//returns time correlation function normalized by datapoints per time and by bin volume
    void write(string filename)const;
    float** spherical_fourier(string, int);
    float show(int t, int b)const{return normalized()[t][b];};		//return normalized correlation at a given time and bin
    int show_n_bins()const{return n_bins;};				//return number of bins
    float show_bin_size()const{return bin_size;};			//return bin size
};
}

#endif