/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class Vector_Self_Correlation: calculates the time autocorrelation function of vectors connecting user-specified atoms within the same molecule */
/*Written by David S. Simmons*/

#ifndef VECTOR_AUTOCORRELATION
#define VECTOR_AUTOCORRELATION

#include <string.h>
#include "analysis.h"
#include "system.h"

namespace std{
	
class Vector_Autocorrelation: public Analysis
{ 
    int n_times;		//number of times in correlation function
    float * timetable;		//table of times corresponding to correlation data
    float * correlation;	//correlation function
    float * orientational_correlation;
    int * weighting;		//weighting for normalization of correlation function
    float mean_length;
    
    //data specifying vectors (bonds) between atoms to be calculated
    int n_vectors;
    int * vector_specieslist;
    int * vector_type1list;
    int * vector_type2list;
    int * vector_index1list;
    int * vector_index2list;
    
    void initialize(System*, string);
    void calculate_mean_vector_length();

    
  public:
    
    Vector_Autocorrelation();
    Vector_Autocorrelation(System*,string);
    Vector_Autocorrelation(const Vector_Autocorrelation &);
    
    ~Vector_Autocorrelation();
    
    Vector_Autocorrelation operator = (const Vector_Autocorrelation &);
    
    void set(System*,string);
    void write(string filename);
    void list_displacementkernel(int,int,int);
    
};

}

#endif