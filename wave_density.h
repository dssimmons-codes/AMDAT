/*Density_Fourier Class definition.  This class calculates and stores the spacial fourier transform of the number density for each wavevector passed from a Wave_Vectors object.*/
/*Written by David S. Simmons*/

#include <complex>
#include "wave_vectors.h"
#include "system.h"
#include "analysis.h"
#include <string>
#include "particle_list.h"

#ifndef WAVE_DENSITY
#define WAVE_DENSITY
namespace std{

class Wave_Density:public Analysis
{
    Wave_Vectors const * wavevectors;	//object storing sorted wavevectors
    complex<double> *** density;	//data array of wave-space densities
    int first_wavenumber_index;		//lower index of wavevector bins to request from Wave_Vectors object
    int last_wavenumber_index;		//upper index of wavevector bins to request from Wave_Vectors object
    
    int n_wavenumbers;			//number of bins in which wavevectors are placed; for a spherically symmetric system these correspond to wavenumbers
    int n_times;			//number of configurations in trajectory
    int n_atoms_looped;			//counts number of atoms looped over to generate density
    int * n_atoms;
    bool sf;				//changes time definition for density calculation of static S(q)
    //calculation variable
    int currenttime;
    string atomset;			//stores set of atoms for which this wave density is calculated, for later information
    
    
    void displacementkernel(int,int,int,int,int,int,int){};
    void clear_memory();

  public:
    Wave_Density();
    ~Wave_Density();									//destructor
    Wave_Density(System * sys, const Wave_Vectors * wv, bool stf=0);
    Wave_Density(System * sys, const Wave_Vectors * wv, int inner, int outer);		//constructor
    Wave_Density(const Wave_Density &);							//copy constructor
    Wave_Density operator = (const Wave_Density &);					//equality operator
    
    void set(System * sys, const Wave_Vectors * wv);
    void set(System * sys, const Wave_Vectors * wv, int inner, int outer);		
    
    int show_n_atoms(int timeii)const{return n_atoms[timeii];};
    const Wave_Vectors * show_wavevectors()const{return wavevectors;};		//return pointer to Wave_Vector object
    int show_first()const{return first_wavenumber_index;};
    int show_last()const{return last_wavenumber_index;};

    complex<double> show_density(int time, int wavenumberindex, int vectorindex)const{return density[time][wavenumberindex][vectorindex];};
    int show_n_atoms_looped()const{return n_atoms_looped;};
    
    void atomlist_kernel(int expii, int species_index, int molecule_index, int atom_type, int atom_index, Particle_List*);
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int){};
    void listkernel(Trajectory *);
    void postprocess_list(){};
    
    void postprocess();
    void postprocess_list(Particle_List*);
    
    void atomkernel(Trajectory * traj);
    
};

}

#endif