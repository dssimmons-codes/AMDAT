/*Intermolecular_Contacts class: identifies and stores list of particles that are contacting (within specified distance threshold) a particle in another molecule*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#ifndef INTERMOLECULAR_CONTACTS
#define INTERMOLECULAR_CONTACTS

#include "analysis.h"
#include "everytime_trajectory_list.h"
#include "system.h"

namespace std{

class Intermolecular_Contacts: public Everytime_Trajectory_List
{
  private:
    float threshold;

  public:
    Fast_Particles();
    Fast_Particles(System * sys);
    Fast_Particles(System * sys, float distance_threshold);	//construct to find  contacting particles based on user-specificed distance cutoff

    Fast_Particles operator = (const Fast_Particles &);	//assignment
    Fast_Particles(const Fast_Particles &);		//copy constructor

    void set(System *);		//construct to find fast_particles based on distance cutoff
    void set(System * sys, float distance_threshold);	//set to find particles based on user-specificed distance cutoff

    Analysis_Type what_are_you(){Analysis_Type type = intermolecular_contacts; return type;};		//virtual method to report the type of analysis

    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list(){};
    void bin(int, float);

    void write(string);
    void write(ofstream& output)const;
};

}

#endif
