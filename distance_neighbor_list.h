/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Distance_Neighbor_List class - builds neighbor list from distance thresholding*/
/*Written by David S. Simmons*/



#include "neighbor_list.h"
#include "analysis_onetime.h"

#ifndef DISTANCE_NEIGHBOR_LIST
#define DISTANCE_NEIGHBOR_LIST

namespace std
{

class Distance_Neighbor_List: public Neighbor_List, public Analysis_Onetime
{
    float threshold;
    
    float ** sigmatrix;		//stores particle sizes
    int n_atomtypes;
    
    void allocate_sig_matrix(string sig_file);
  
  public:
    
    Distance_Neighbor_List();
    Distance_Neighbor_List(const Distance_Neighbor_List&);
    ~Distance_Neighbor_List();
    Distance_Neighbor_List operator=(const Distance_Neighbor_List&);
    
    Distance_Neighbor_List(System* sys, float thresh, string sigmatrixname, int firsttime=-1, int lasttime=-1);
    
    Analysis_Type what_are_you(){Analysis_Type type = distance_neighbor_list; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list(){};
    
    void write(string){};
    void write(ofstream& output)const{};
  
};

  
}

#endif