


#ifndef STRING_MULTIBODIES
#define STRING_MULTIBODIES

#include "multibody.h"
#include "provisional_multibodies.h"

namespace std{
  
class String_Multibodies: public Analysis: public Provisional_Multibodies
{
    int timegap;
    
    float ** sigmatrix;		//stores particle sizes
    int n_atomtypes;
    
    float threshold;
    
    int trajectories_considered;
    
    //data members that change during calculation and have no permanent value
    int * stringID;
    int trajindex;
    vector<bool> string_validity;	//during string computation, stores 1 if a string is good and 0 if it should be discarded because it has been merged into another string
    int mass_switch_ID(int oldID,int newID);	//change all stringIDs with value oldID to newID
    int n_trajectories;
    
    
    void allocate_sig_matrix(string);
  public:
    
    String_Multibodies();
    String_Multibodies(const String_Multibodies&);
    ~String_Multibodies();
    String_Multibodies operator=(const String_Multibodies&);
    
    String_Multibodies(System * syst, int tgap, float thresh, string sigmatrixname);
    
    
    
    void analyze(Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*,int,int,int);
    void postprocess_list();
    
  
}
  
  
  
}


#endif