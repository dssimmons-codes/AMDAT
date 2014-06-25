/*String class - Identifies strings in fluid*/
/*Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#ifndef STRINGS
#define STRINGS


#include <string>
#include "analysis.h"

namespace std {

class Strings: public Analysis
{
    
  int n_timegaps;			//number of timegaps to be considered
  int timegap;			//timegap index
  float timetable;			//stores timegap information				
  int trajectories_considered;		
  float ** sigmatrix;			//stores particle sizes
  int n_atomtypes;
  
  int maxstrings;		//maximum number of strings per timegap per time
  int maxstringatoms;		//maximum number of atoms per string
  
  int n_times;			//number of times at which string lists are assembled for each timegap
  int n_trajectories;
  float threshold;			//maximum distance threshold for string definition
  int ** stringID;			//list of string ID for each atom in system at each start time under consideration. The first index is the timegap, the second index is the starting time, the third index is the unique atomID of the atom.
  int * stringcount;			//next open string ID at given timegap and starttime
  int *** stringlist;			//list of strings, with each string consisting of the atomIDs of the contituent atoms. First index is timegap, second is initial time, third is stringID, and fourth is index of atom within string.
  int ** stringsize;			//size of each string. First index is timegap, second is initial time, third is stringID.
  
  float mean_length;			//average length of strings as a function of timegap
  float mean_strings;			//mean number of strings as a function of timegap
  float * length_distribution;		//distibution of string lengths at each timegap
  float order_parameter;		//order parameter comprised of fraction of trajectories in strings
  
  
  /*these are internal processing variables and have no fixed value*/
  int thistime;				//stores current time
  int nexttime;				//stores later time based on timegap
  int current_timegap;
  
  void create_string(int,int,int);	//method to create new string
  void merge_strings(
  int,int,int);	//method to combine two strings
  
  void defragment_strings();		//eliminate empty strings and defragment string arrays
  void calculate_mean_strings();		//find mean number of strings
  void calculate_mean_length();		//find mean string lengths
  void stringlength_distribution();		//calculate distribution of string lengths
  void fraction_in_strings();			//calculate fraction of trajectories in strings
  
  void allocate_sig_matrix(string);
  
  public:

    Strings(System * sys,  int timegap1, float thresh, string sigmatrixname, int maximum_strings= 1000, int maximum_stringatoms = 200);
    void preprocess(){};
    void initialize(System* sys, int timegap1, float thresh, string sigmatrixname, int maximum_strings, int maximum_stringatoms);
    void write(string)const;		//write string statistics
    void write_trajectories(string);	//write strings to trajectory files, with each time block getting its own file
    void write_jump(string, int)const;		//write all two timeframes for each string to file, with the first frame being the first frame of the block and the second frame being specified by the value of last_timegap.
    
    void analyze(Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*);
    void postprocess_list();
    
    void construct_trajectory_list(Trajectory_List*);
};


}

#endif