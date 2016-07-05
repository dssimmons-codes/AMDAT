/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Voronoi_Neighbor_List class - builds neighbor list from distance thresholding*/
/*Written by David S. Simmons*/



#include "neighbor_list.h"
#include "analysis_onetime.h"
#include "voro++.hh"

#ifndef VORONOI_NEIGHBOR_LIST
#define VORONOI_NEIGHBOR_LIST

namespace std
{

class Voronoi_Neighbor_List: public Neighbor_List, public Analysis_Onetime
{
    container * voronoi;
    int n_x,n_y,n_z;

  public:
    
    Voronoi_Neighbor_List();
    Voronoi_Neighbor_List(const Voronoi_Neighbor_List&);
    ~Voronoi_Neighbor_List();
    Voronoi_Neighbor_List operator=(const Voronoi_Neighbor_List&);
    
    Voronoi_Neighbor_List(System* sys, int firsttime=-1, int lasttime=-1);
    
    Analysis_Type what_are_you(){Analysis_Type type = distance_neighbor_list; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int){};
    void postprocess_list();
    
    void write(string){};
    void write(ofstream& output)const{};
  
};

  
}

#endif