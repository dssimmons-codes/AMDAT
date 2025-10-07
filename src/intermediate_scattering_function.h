// Copyright (C) 2010-2025 David S. Simmons

/*Written by David S. Simmons*/

#include "wave_density.h"
#include "correlation_2d.h"

namespace std{

class Intermediate_Scattering_Function:public Correlation_2D
{
    Wave_Density * wavedensity1;
    Wave_Density * wavedensity2;
    void atomkernel(int,int,int,int){};


  public:
    Intermediate_Scattering_Function();
    ~Intermediate_Scattering_Function(); //destructor
    Intermediate_Scattering_Function(const Intermediate_Scattering_Function &); //copy constructor
    Intermediate_Scattering_Function operator = (const Intermediate_Scattering_Function &); //assignment operator

    Intermediate_Scattering_Function(System * sys, Wave_Density * wd1, Wave_Density * wd2, bool parallel=1);		//calculates for all wavenumbers
    Intermediate_Scattering_Function(System * sys, Wave_Density * wd1, Wave_Density * wd2, int wn1, bool parallel=1);		//calculates for a single wavenumber
    Intermediate_Scattering_Function(System * sys, Wave_Density * wd1, Wave_Density * wd2, int wn1, int wn2, bool parallel=1);		//calculates for a range of wavenumbers
    Intermediate_Scattering_Function(System * sys, int timegapii, Wave_Density * wd1, Wave_Density * wd2, bool parallel=1);		//calculates for all wavenumbers at a particular timegap.  Setting timegapii=0 gives the structure factor.
    //Intermediate_Scattering_Function(System* sys, Particle_List * particle_list, Wave_Density * wd2, int wn1, int wn2);		/*Method to calculate the intermediate scattering function at all timegaps, for a range of wavenumbers bounded by wn1 and wn2, based on a time-dependent particle list given by particle_list.  It will base the wavevectors and system off of the wavedensity that is provided.*/

    
    Analysis_Type what_are_you(){Analysis_Type type = intermediate_scattering_function; return type;};		//virtual method to report the type of analysis

    //void displacementkernel(int timegapii, int thisii, int nextii, int species_index, int molecule_index, int atom_type, int atom_index);

    void displacementkernel(int timegapii, int thisii, int nextii, Trajectory * traj);

};

}
