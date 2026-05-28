// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_TrajList_Intersection class - stores list of multibody objects that contain beads within a given trajectory_list*/
/*Written by David S. Simmons*/


#ifndef MULTIBODY_TRAJLIST_INTERSECTION
#define MULTIBODY_TRAJLIST_INTERSECTION

#include "coordinate.h"
#include "multibody_list.h"
#include "multibody_analysis.h"
#include "trajectory_list.h"

namespace std{

class Multibody_TrajList_Intersection: public Multibody_List, public Multibody_Analysis
{
  private:

    int threshold;
    Trajectory_List * trajectory_list;

  public:


    Multibody_TrajList_Intersection();
    Multibody_TrajList_Intersection(System*syst);
    Multibody_TrajList_Intersection(System*syst, Trajectory_List* trajlist, int thresh);
    Multibody_TrajList_Intersection(const Multibody_TrajList_Intersection &);
    Multibody_TrajList_Intersection operator=(const Multibody_TrajList_Intersection &);

    void analyze(Multibody_List * mblist);
    void list_displacementkerneyl(int,int,int){};
    void listkernel(Multibody*,int,int,int);
    void postprocess_list();

    void write(string) const;

};

}


#endif
