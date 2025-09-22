/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/



/*Need to add mimplementation of methods and child classes for actually creating these objects - two basic approaches are distance thresholding and voronoi*/

#include <vector>
#include "boolean_list.h"
#include "system.h"
#include "trajectory_list.h"
#include "value_list.h"

#ifndef NEIGHBOR_LIST
#define NEIGHBOR_LIST

namespace std
{
  

  
class Neighbor_List: public Value_List<float>
{
  protected:
    vector<vector<vector<Trajectory*>>> neighbors;	//indices are time, base trajectory, neighbor trajectory
    //mutable Boolean_List * included;	//array of boolean lists specifying which trajectories are in value list at each time: [internal_time]
    
    vector<bool> computed_times;
    

  public:
    Neighbor_List();
    Neighbor_List(const Neighbor_List&);
    ~Neighbor_List();
    Neighbor_List operator=(const Neighbor_List&);
    
    Neighbor_List(System*sys);
    
    bool is_neighbor(int timeii, int trajii, Trajectory* trajcheck)const;	//returns true if trajcheck is a neighbor of trajectory indexed by trajii at time timeii; false otherwise

    vector<Trajectory*> show_neighbors(int trajii, int time1)const;	//returns vector of trajectories in a particle's neighborlist at a given time
    vector<Trajectory*> persistent_neighbors(int trajii, int time1, int time2)const;	//returns vector of trajectories in a particle's neighborlist at both of two times
    int n_persistent_neighbors(int trajii, int time1, int time2)const; //returns number of trajectories in a particle's neighborlist at both of two times
    
    void read_bonds(string filename); //read a list of bonds in a custom LAMMPS format and make each bonded atom a neighbor of its bonded atoms (i.e. neighbors denote bonds)
    
    void neighborloop(Analysis* analysis, int timii, int trajii){};
    int show_n_neighbors(int timeii, int trajii)const{return neighbors[timeii][trajii].size();};	//returns the number of neighbors for a given trajectory
    
    void write_statistics(string filename, int n_moments)const;
    
};
}

#endif
