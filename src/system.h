/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*System class: holds information on system-wide properties.*/
/*Written by David S. Simmons*/

#ifndef SYSTEM
#define SYSTEM
//#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>

#include "molecule.h"
#include "boolean_list.h"
#include "trajectory_set.h"
#include "trajectory.h"
#include "analysis.h"
#include "multibody_analysis.h"
#include "tokenize.h"
#include "vector_map.h"


namespace std {

class Multibody_Set;

template <class valType>
class Value_List;

class System
{
  protected:
    bool np;			//true if system is np (not constant volume), false if it is nv (constant volume)

    int n_atoms;				//total number of atoms in system
    int n_species;				//number of molecular species
    int n_atomtypes;				//number of types of atoms in system
    string * species_name;
    vector <string> atomtype_name;
    int * atoms_per_species;
    int total_molecules;			//total number of molecules in system
    int * n_molecules;				//array of number of molecules of each species
    int * moleculecountholder;
    int total_trajectories;
    string timetype;				//stores type of time spacing
    int n_timesteps;				//number of timesteps
    int n_exponential_steps;			//number of exponential timesteps taken before starting next exponential block
    int n_exponentials;				//number of exponential time blocks
    int n_timegaps;				//number of distinct timegaps between configuration times
    float * timelist;				//array of trajectory times
    float time_unit;				//numerical unit of time
    float exp_base;				//base of exponential scaling
    bool frt;					//controls whether first time step is zero or one; this has very large implications on the way that displacements are calculated throughout the code
    int first_exponent;				//set value of first exponent in exponential timescheme.  This is usually either zero or one, and it not count the zeroth time, if one exists.
   

    Vector_Map<string,Value_List<float>*> value_lists;      //stores value lists

    Molecule ** molecules;			//array of molecules in rows by type

    Atom_Trajectory ** atomlist;		//direct list of trajectories by atomID
    Molecule ** moleculelist;			//direct list of molecules by molecule ID
    vector <Trajectory *> trajectorylist;		//direct list of trajectories by trajectory ID

    int n_multibodies;
    //vector <Multibody_Set*> multibody_sets;	//store sets of multibodies defined by user
    Vector_Map <string, Multibody_Set*> multibody_sets;
    Vector_Map <string, Trajectory_Set*> trajectory_sets;   //stores sets of trajectories defined by user


    bool unwrapped;				//are unwrapped coordinates defined?
    bool wrapped;				//are wrapped coordinate defined?
    bool boxified;				//has system been boxified?

    Coordinate * box_size;			//size of box
    Coordinate ** box_boundary;			//boundaries of box

    float * rho;					//mean system density

    void create_molecules(int **);		//method to create molecule objects
    void create_timelist();			//calculate timelist based on input parameters

    int displacement_limit;			//sets limit on number of time-datapoints used for linear time-scaled displacement times


    
    /*methods to read in various trajectory formats*/
    void read_trajectory(string trajectory_type, vector<string> file_in, string fileline);
    void count_atoms(int **);					//method to calculate number of atoms
    void xyz_prep(vector<string> file_in, string fileline);
    void xyz_prep_withlog(vector<string> file_in, string fileline);
    void read_xyz(string);					//method to read in trajectory data in standard xyz format
    void read_xyz(string, string);					//method to read in trajectory data in standard xyz format with an extra xyz structure file listing order of molecular blocks in trajectory
    void read_position(string stem);				//method to read in sequence of binary trajectory files, with each file corresponding to one exponential block
    void xtc_prep(vector<string> file_in, string fileline);
    int ** read_gro(string);				//method to read .gro file to determine trajectory file structure
    void read_xtc_format(string filename, int** atomidentifier);
    void custom_prep(vector<string> file_in, string fileline);
    void custom_byid_prep(vector<string> file_in, string fileline);
    void read_custom(string);
    void read_custom(string, string);
    void read_custom_byid(string);
    
    void custom_manual_prep(vector<string> file_in, string fileline);
    void read_custom_manual(string, string);
    //void read_custom_manual(string, string, string); TODO add this later
    //TODO: should implement manual for the by_id versions as well



    bool floatCompare(float, float);
    
    /*String-handling methods*/
    Tokenize tokenize;
    bool in_string_array(string * tokens, int array_size, string target);
    bool in_string_array(vector <string> tokens, string target);
    int find_in_string_array(string * tokens, int array_size, string target);
    int find_in_string_array(vector <string> tokens, string target);

  public:
    System();							//default constructor that queries use for necessary values
    System(vector<string> file_in ,bool ensemble=0);
    //~System();
    void clear_memory();

    /*Method to add velocity data from separate custom file*/
    void read_velocity_byid(string);
    
    /*------Methods to perform some sort of operation on entire system-------*/
    void unwrap();				//send command to all molecules to unwrap atom trajectories
    void wrap();				//send command to all molecules to wrap atom trajectories
    void set_limit(int limit){displacement_limit=limit;};	//set limit on number of time points looped over for displacement times within linear time scaling
    void boxify();						//method to check every coordinate to see if it is within the system boundaries; if not, it treats it as an image coordinate and replaces it with the 'real' coordinate


    /*Methods to handle value lists*/
    Value_List<float>* find_value_list(string,bool allow_nofind=0)const;
    void add_value_list(Value_List<float>*, string);
    void delete_value_list(string);

    
    /*------Methods to handle multibodies and multibody_sets--------*/
    Multibody_Set* create_multibody_set(string setname, int n_args, string * args);
    Multibody_Set* create_multibody_set();			//creates a multibody_set containing a multibody for each molecule in the system, with each multibody containing all the trajectories in the corresponding molecule
    Multibody_Set* create_multibody_set(int speciesii);	//creates a multibody_set containing a multibody for each molecule of a given species, with each multibody containing all the trajectories in the corresponding molecule
    Multibody_Set* create_multibody_set(int speciesii, int type);//creates a multibody_set containing a multibody for each molecule of a given species, with each multibody containing all the trajectories of a specified type in the corresponding molecule
    Multibody_Set* create_multibody_set(int speciesii, int n_trajectories, int * type, int * index);//creates a multibody_set containing a multibody for each molecule of a given species, with each multibody containing n_trajectories specified by arrays providing the type and index of each trajectory.

    void add_multibody_set(string,Multibody_Set*);
    Multibody_Set* find_multibody_set(string, bool allow_nofind=0)const;
    void delete_multibody_set(string);


    Trajectory_Set* create_trajectory_set(string setname, string multibodysetname, string traj_typename,  bool centertype);	//centertype = 1 for centroid or 0 for COM
    Trajectory_Set* find_trajectory_set(string setname, bool allow_nofind)const;
    void add_trajectory_set(string trajectory_set_name,Trajectory_Set* trajectory_set);
    void add_trajectories (Trajectory_Set * new_trajectories);

    /*--------Methods to return information about box/system---------*/
    Coordinate size()const{return box_size[0];};		//return size of system;
    Coordinate size(int time){return box_size[time];};
    const Coordinate* time_dependent_size()const{return box_size;};
    Coordinate min_box_dimensions()const;
    const Coordinate* boundaries()const{return box_boundary[0];};		//return boundaries of system
    const Coordinate * boundaries (int time) const {return box_boundary[time];};
     Coordinate ** time_dependent_boundaries() const {return box_boundary;};
    float show_rho()const{return rho[0];};
    float show_rho(int timeii)const{return rho[timeii];};

    /*--------Methods to return trajectories and information about trajectories---------*/
    int show_n_atoms()const{return n_atoms;};		//return total number of atoms in system
    int show_n_molecules()const{return total_molecules;};
    int show_n_molecules(int speciesii) const {return n_molecules[speciesii];};
    int show_n_trajectories()const{return trajectorylist.size();};
    Atom_Trajectory* show_atom(int atomID)const{return atomlist[atomID];};		//return pointer to atom via atomID
    Molecule* show_molecule(int moleculeID)const{return moleculelist[moleculeID];};	//return pointer to molecule via moleculeID
    Trajectory* show_trajectory(int trajectoryID)const{return trajectorylist[trajectoryID];};	//return pointer to trajectory via trajectoryID
    Molecule* show_molecule(int species_index, int molecule_index)const{return &molecules[species_index][molecule_index];};		//return molecule from master list of molecules
    
    int show_species_index(string) const;
    int show_atomtype_index(string) const;
    bool atomtype_exists (string) const;
    int add_atomtype(string);
    int show_n_atomtypes()const{return n_atomtypes;};	//return number of atomtypes including only atomtypes existing at time of trajectory file read in
    int show_n_expanded_atomtypes()const {return atomtype_name.size();};	//return number of trajectory types including newly defined trajectories
    
    Coordinate show_unwrapped(int species_index, int molecule_index, int atom_type, int atom_index, int timestep) const;	//returns a single unwrapped coordinate of a given atom
    


    /*-------------Methods to return information about time scheme------------*/
    int show_n_timesteps()const{return n_timesteps;};		//return number of timesteps in system
    int* timegap_weighting(bool fullblock=1)const;			//return number of time spacings corresponding to each timegap
    int timegap_weighting(int,bool fullblock=1)const;			//return number of time spacings corresponding to given timegap
    int show_n_timegaps()const{return n_timegaps;};			//return number of timegaps
    int show_n_exponentials()const{return n_exponentials;};		//return number of exponential blocks in time scheme
    float show_exp_base()const{return exp_base;};			//return exponential base of time scheme
    int show_n_exponential_steps()const{return n_exponential_steps;};	//return number of exponential steps per block in time scheme
    float * displacement_times() const;	//calculate array of times corresponding to displacement timesteps
    float displacement_times(int) const;	//calculate array of times corresponding to displacement timesteps
    float show_time(int tii)const{return timelist[tii];};	//return time spacing corresponding to  time spacing index
    void show_times(int ntimes, int * timeindices, float * times)const;	//return array of times given array of time indices
    float show_time_unit()const{return time_unit;};		//return basic unit of time in time scheme
    bool show_frt()const{return frt;};			//show whether first timestep is 0 or 1
    int blockstart(int blockii)const{return blockii*n_exponential_steps;};		//return index corresponding to start of block
    int block(int blockii,int* timelist)const;	//calculates indices of times in block and stores in timelist and number of indices is returned value
    void big_block(int blockii, int* timelist)const;	//calculates  indices of all times that lead to a valid timegap from start of block, including those in latter blocks. Stores this list in timelist.
    int big_blocksize(int blockii)const{return n_exponential_steps+n_exponentials-blockii;};	//return number of timegaps accessible from start of givem block
    string show_timetype()const{return timetype;};		//returns type of time scheme
    int show_displacement_limit(){return displacement_limit;};

    /*---------Loop methods that take analysis object as argument---------*/

    /*loops over time spacings for use with system loops*/
    void displacement_loop(Analysis*, Trajectory *, bool fullblock=1)const;
    void displacement_loop(Analysis*, Trajectory *, int, bool fullblock=1)const;

    /*loops over time spacings for use with new trajectory lists*/
    void displacement_list(Analysis*, bool fullblock=0)const;
    void displacement_list(Analysis*, int timegap, bool fullblock=0)const;
    void displacement_list(Analysis*, int timegap, int firstblock, int lastblock, bool fullblock=0)const;
    
    void displacement_list(Multibody_Analysis*, bool fullblock=0)const;
    void displacement_list(Multibody_Analysis*, int timegap, bool fullblock=0)const;
    void displacement_list(Multibody_Analysis*, int timegap, int firstblock, int lastblock, bool fullblock=0)const;

    /*Methods to perform analyses on subsets of atoms in system*/
    void loop_species_moleculecom(Analysis* analysis, int species_index)const;
    void loop_atom_species(Analysis* analysis,int species_index, int atomtype, int atomindex)const;
    void loop_type_molecule(Analysis* analysis, int species_index, int molecule, int atomtype)const;
    void loop_molecule(Analysis* analysis, int species_index, int molecule)const;
    void loop_species(Analysis* analysis, int species_index)const;
    void loop_type_species(Analysis* analysis,int species_index, int atomtype)const;
    void loop_type_system(Analysis* analysis, int atomtype)const;
    void loop_system(Analysis* analysis)const;

    /*Methods to write new trajectory files from system data*/
    void write_starr()const;
    void write_single_particle(int trajii, string filename)const;

};

}

#endif
