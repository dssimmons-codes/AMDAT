/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#include <string>
#include "analysis.h"
#include "wave_vectors.h"
#include <string>
#include "system.h"

#ifndef CORRELATION_2D
#define CORRELATION_2D

namespace std{

class Correlation_2D:public Analysis	
{
    virtual void atomkernel(int,int,int,int){};				//dummy method
    virtual void displacementkernel(int,int,int,int,int,int,int){};	//dummy method
    
  protected:
    float ** correlation;
    float * timetable;
    
    Wave_Vectors const * wavevectors;
    
    int * timegap_weighting;
    int n_atoms_represented;			//holds number of atoms for normalization when this number is the same for all displacement times
    float * n_atoms;				//more generally holds number of atoms for normalization at each time.  It is this quantity which is ultimately used for normalization.
    
    int n_times;
    int firsttime;
    int lasttime;
    
    int n_spacebins;
    int first_bin_index;
    int last_bin_index;
    float bin_size;
    
    virtual void postprocess();
    virtual void postprocess_list();
    
  public:
    Correlation_2D();
    Correlation_2D(const Correlation_2D &);
    //~Correlation_2D();
    Correlation_2D operator =(const Correlation_2D &);
    
    Analysis_Type what_are_you(){Analysis_Type type = correlation_2d; return type;};		//virtual method to report the type of analysis
    
    Correlation_2D operator + (const Correlation_2D &) const;
    Correlation_2D operator - (const Correlation_2D &) const;
    Correlation_2D operator & (const Correlation_2D &) const;
    Correlation_2D operator | (const Correlation_2D &) const;
    
    void write(string filename) const;
    void write(ofstream& output) const;
};

}

#endif