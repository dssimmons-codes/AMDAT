/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate multibody mean gyration radius*/
/*Written by David S. Simmons*/

#ifndef GYRATION_RADIUS
#define GYRATION_RADIUS

#include <sstream>
#include <string>

#include "multibody_analysis.h"

namespace std{

class Gyration_Radius: public Multibody_Analysis
{
    float gyration_radius;
    int weighting;
    void initialize(System*);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

    
    
  public:
    Gyration_Radius();			//default constructor
    Gyration_Radius(const Gyration_Radius &);		//copy constructor
    Gyration_Radius operator = (const Gyration_Radius &);	//assignment
    ~Gyration_Radius(){};
    
    Gyration_Radius(System*);
    
    //Analysis_Type what_are_you(){Analysis_Type type = gyration_radius; return type;};		//virtual method to report the type of analysis
    
    void set(System * sys){initialize(sys);};
    
    void analyze(Multibody_List * mblist);
    void listkernel(Multibody *, int, int, int);
    void postprocess();
    
    void write(string);
    void write(ofstream&)const;
    
    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();
    
    float show()const{return gyration_radius;};			//method to return gyration radius
//	bool isThreadSafe(){return true;};
};
}

#endif
