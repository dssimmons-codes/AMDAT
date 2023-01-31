/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate list of particle displacements*/
/*Written by David S. Simmons*/


#ifndef DISPLACEMENT_LIST
#define DISPLACEMENT_LIST

#include "analysis.h"
#include "value_list.h"
#include <string>
#include <sstream>

namespace std{

/*need to add ability to set timegap index*/
  
class Displacement_List : public Value_List<float>, public Analysis
{
  int timegap;
  
  public:
    Displacement_List();			//default constructor
    Displacement_List(const Displacement_List &);		//copy constructor
    Displacement_List(System*);
    Displacement_List(System*, int);					//initialize with int specifying timegap
    Displacement_List operator = (const Displacement_List &);	//assignment
    
    Analysis_Type what_are_you(){Analysis_Type type = displacement_list; return type;};		//virtual method to report the type of analysis
    
    float * normalized();
    void write(string)const;
    void write(ofstream&)const;
    
    void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
//	bool isThreadSafe(){return true;};
};
}

#endif