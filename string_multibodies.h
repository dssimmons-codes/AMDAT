


#ifndef STRING_MULTIBODIES
#define STRING_MULTIBODIES

#include "multibody.h"

namespace std{
  
class String_Multibodies: public Analysis
{
    int timegap;
    float ** sigmatrix;		//stores particle sizes
    float threshold;
    
    int trajectories_considered;
  
    vector<vector<Multibody>> string_multibodies;
    
    //data members that change during calculation and have no permanent value
    int * stringID;
    int trajindex;
    vector<bool> string_validity;	//during string computation, stores 1 if a string is good and 0 if it should be discarded because it has been merged into another string
    int mass_switch_ID(int oldID,int newID);	//change all stringIDs with value oldID to newID
    int n_trajectories;
    
    //need to add the following
    float mean_length;
    float mean_strings;
    float * length_distribution;
    float order_parameter;		
    
  public:
    void analyze(Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*,int,int,int);
    void postprocess_list();
    
  
}
  
  
  
}


#endif