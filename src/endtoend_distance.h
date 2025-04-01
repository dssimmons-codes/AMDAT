/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate multibody mean end-to-end distance*/
/*Written by Pierre Kawak*/

#ifndef ENDTOEND_DISTANCE
#define ENDTOEND_DISTANCE

#include <sstream>
#include <string>

#include "multibody_analysis.h"

namespace std{

class Endtoend_Distance: public Multibody_Analysis
{
    float endtoend_distance;
    float endtoend_distance_x;
    float endtoend_distance_y;
    float endtoend_distance_z;
    int weighting;
    void initialize(System*);


    float * ree_by_n;
    float * ree_x_by_n;
    float * ree_y_by_n;
    float * ree_z_by_n;
    int max_n;
    int * weighting_by_n;

  public:
    Endtoend_Distance();			//default constructor
    Endtoend_Distance(const Endtoend_Distance &);		//copy constructor
    Endtoend_Distance operator = (const Endtoend_Distance &);	//assignment
    ~Endtoend_Distance(){};

    Endtoend_Distance(System*);

    //Analysis_Type what_are_you(){Analysis_Type type = endtoend_distance; return type;};		//virtual method to report the type of analysis

    void set(System * sys){initialize(sys);};

    void analyze(Multibody_List * mblist);
    void listkernel(Multibody *, int, int, int);
    void postprocess();

    void write(string);
    void write(ofstream&)const;

    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();

    float show()const{return endtoend_distance;};			//method to return end-to-end distance
//	bool isThreadSafe(){return true;};
};
}

#endif
