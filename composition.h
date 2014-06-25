/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*MClass containing information about system composition*/
/*Written by Daniel A. Hunsicker*/


#ifndef COMPOSITION
#define COMPOSITION

#include "analysis.h"
#include "system.h"
#include <string>

namespace std{


class Composition : public Analysis
{
    System * system;
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
    Composition(System*,int n_xbins, int n_ybins, int n_zbins, float lx, float ly, float lz);
    Composition(const Composition &);
    ~Composition();

    Composition operator=(const Composition &);
    
    Analysis_Type what_are_you(){Analysis_Type type = composition; return type;};
    
     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void listkernel(Trajectory *);

    void write(string);
    bool isThreadSafe(){return false;};
    
    
    

};
}
#endif // COMPOSITION
