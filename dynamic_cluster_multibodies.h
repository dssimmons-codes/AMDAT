/*AMDAT: Amorphous Molecular Dynamics Analysis Toolkit*/
/*Dynamic_Cluster_Multibodies class: builds multibodies based on any pairwise dynamic criterion*/
/*David S Simmons*/

#ifndef DYNAMIC_CLUSTER_MULTIBODIES
#define DYNAMIC_CLUSTER_MULTIBODIES

#include "analysis.h"
#include "provisional_multibodies.h"

namespace std{

class Dynamic_Cluster_Multibodies: public Analysis, public Provisional_Multibodies
{
  protected:
    int timegap;
    int trajectories_considered;
    
    //data members that change during calculation and have no permanent value
    int * multibodyID;
    int trajindex;
    int n_trajectories;
    int max_trajectories;
    vector<bool>string_validity;
    
  
    int mass_switch_ID(int oldID,int newID);	//change all stringIDs with value oldID to newID
    virtual bool clustered_check(Trajectory*, Trajectory*, int, int){return false;};	//virtual method used to determine whether two trajectories are in the same multibody at a given time 
    
  public:
    Dynamic_Cluster_Multibodies();
    Dynamic_Cluster_Multibodies(const Dynamic_Cluster_Multibodies&);
    ~Dynamic_Cluster_Multibodies();
    Dynamic_Cluster_Multibodies operator=(const Dynamic_Cluster_Multibodies&);
    
    Dynamic_Cluster_Multibodies(System*, int);
    
    void analyze(Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*,int,int,int);
    void postprocess_list(){};

};

}

#endif