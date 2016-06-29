/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/



/*Need to add mimplementation of methods and child classes for actually creating these objects - two basic approaches are distance thresholding and voronoi*/

#include "boolean_list.h"
#include "system.h"
#include "trajectory_list.h"

namespace std
{
  
class Neighbor_List
{
  protected:
    System * syst;
    
    vector<vector<vector<Trajectory*>>> neighbors;	//indices are time, base trajectory, neighbor trajectory
    mutable Boolean_List * included;	//array of boolean lists specifying which trajectories are in value list at each time: [internal_time]
    int n_times;
    
    bool threshold_neighborcount(int, int, bool, int)const;
    bool threshold_neighborcount(int, int, int, int)const;
    void update_size();
  public:
    Neighbor_List();
    Neighbor_List(const Neighbor_List&);
    ~Neighbor_List();
    Neighbor_List operator=(const Neighbor_List&);
    
    bool is_neighbor(int timeii, int trajii, Trajectory* trajcheck)const;	//returns true of trajcheck is a neighbor of trajectory indexed by trajii at time timeii; false otherwise


    //Not yet implemented
    vector<Trajectory*> persistent_neighbors(int trajii, int time1){};	//returns vector of trajectories in a particle's neighborlist at a given time
    
    vector<Trajectory*> persistent_neighbors(int trajii, int time1, int time2){};	//returns vector of trajectories in a particle's neighborlist at both of two times
    
    neighborloop(Analysis* analysis, int timii, int trajii){};
    
};
}