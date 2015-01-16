/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Analysis class identifying system edge as function of time*/
/*Written by David S. SImmons*/


#ifndef EDGE_DETECTOR_TIMEDEPENDENT
#define EDGE_DETECTOR_TIMEDEPENDENT

#include "analysis.h"
#include "system.h"
#include <string>
#include "coordinate.h"

namespace std{


class Edge_Detector_Timedependent : public Analysis
{
    int n_times;
    Coordinate * time_dependent_edge;
    Coordinate unit_vector;
    
    float max_magnitude;


    public:
    Edge_Detector_Timedependent();
    Edge_Detector_Timedependent(System*,Coordinate);
    Edge_Detector_Timedependent(const Edge_Detector_Timedependent &);
    ~Edge_Detector_Timedependent();

    Edge_Detector_Timedependent operator=(const Edge_Detector_Timedependent &);
    
    Analysis_Type what_are_you(){Analysis_Type type = edge_detector_timedependent; return type;};
    
     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void listkernel(Trajectory *,int timegapii, int thisii, int nextii);

    void write(string);
    bool isThreadSafe(){return false;};
    

};
}
#endif // COMPOSITION
