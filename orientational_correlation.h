/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate mean orientational correlation between specified multibodies and a specified fixed vector*/
/*Written by David S. Simmons*/

#ifndef ORIENTATIONAL_CORRELATION
#define ORIENTATIONAL_CORRELATION

#include <sstream>
#include <string>

#include "multibody_analysis.h"

namespace std{


class Orientational_Correlation: public Multibody_Analysis
{
    float * correlation;
    float * weighting;
    float overall_correlation;
    int n_times;
    Coordinate correlated_vector;

public:
    Orientational_Correlation();
    Orientational_Correlation(const Orientational_Correlation &);
    Orientational_Correlation operator = (const Orientational_Correlation &);

    Orientational_Correlation(System*);
    Orientational_Correlation(System*, Coordinate);

    void analyze(Multibody_List * mblist);
    void list_displacementkernel(int,int,int){};
    void listkernel(Multibody *, int, int, int);
    void postprocess_list();

    void write(string) const;
    void write(ofstream&)const;

};

}

#endif