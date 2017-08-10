/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate multibody mean gyration radius*/
/*Written by David S. Simmons*/

#ifndef BOND_AUTOCORRELATION_FUNCTION
#define BOND_AUTOCORRELATION_FUNCTION

#include <sstream>
#include <string>

#include "multibody_analysis.h"
#include "coordinate.h"

namespace std{

class Bond_Autocorrelation_Function: public Multibody_Analysis
{
    int n_times;
    float * baf;
    int * weighting;
    float * timetable;
    void initialize(System*);
    void initialize(System*, Coordinate dim);  
    void initialize(System*, int l_type, Coordinate dim); 
    int atomcount;
    
    Coordinate dimensions;
    
    typedef Coordinate(Bond_Autocorrelation_Function::*vector_prep)(Coordinate)const;
    vector_prep vprep;
    
    Coordinate prep_inplane(Coordinate)const;
    Coordinate prep_outofplane(Coordinate)const;
    
    typedef float(Bond_Autocorrelation_Function::*polynomial_choice)(float)const;
    polynomial_choice legendre_p;
    
    float legendre_1(float val)const;
    float legendre_2(float val)const;


    
    
  public:
    Bond_Autocorrelation_Function();			//default constructor
    Bond_Autocorrelation_Function(const Bond_Autocorrelation_Function &);		//copy constructor
    Bond_Autocorrelation_Function operator = (const Bond_Autocorrelation_Function &);	//assignment
    ~Bond_Autocorrelation_Function();
    
    Bond_Autocorrelation_Function(System*);
    Bond_Autocorrelation_Function(System*, Coordinate dim);
    
    Bond_Autocorrelation_Function(System*, int l_type, Coordinate dim);
    
    
    //Analysis_Type what_are_you(){Analysis_Type type = gyration_radius; return type;};		//virtual method to report the type of analysis
    
    void set(System * sys){initialize(sys);};
    
    void analyze(Multibody_List * mblist);
    void list_displacementkernel(int,int,int);
    void listkernel(Multibody *, int, int, int);
    void postprocess_list();
    
    void write(string) const;
    void write(ofstream&) const;
    
    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();

//	bool isThreadSafe(){return true;};
};
}

#endif
