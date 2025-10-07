// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*MClass containing information about system composition*/
/*Written by Daniel A. Hunsicker*/


#ifndef COMPOSITION
#define COMPOSITION

#include "analysis.h"
#include "system.h"
#include <string>
#include "analysis_onetime.h"

namespace std{


class Composition: public Analysis_Onetime
{
    protected:
//    System * system;
    int n_atomtypes;
    int n_molecules;
    int n_times;
    int current_time;
    int current_total_atoms;
    int total_atoms;
    float volume;

    float * current_density;
    float average_density;
    float * time_average_comp;
    float** current_comp;




    public:
    Composition();
    Composition(System*,int n_xbins, int n_ybins, int n_zbins, float lx, float ly, float lz, int timescheme=-1);
    Composition(const Composition &);
    ~Composition();

    Composition operator=(const Composition &);
    
    Analysis_Type what_are_you(){Analysis_Type type = composition; return type;};
    
    void listkernel(Trajectory *);
    void timekernel(int);
    void postprocess_list();
    void write(string);
    void write(ofstream&);
    bool isThreadSafe(){return false;};
    
    
    

};
}
#endif // COMPOSITION
