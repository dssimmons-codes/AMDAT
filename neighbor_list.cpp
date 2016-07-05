/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/


#include "neighbor_list.h"


using namespace std;




Neighbor_List::Neighbor_List()
{
  syst=0;
  included = new Boolean_List[1];
  time_conversion = new int [1];
  defined_times = new bool [1];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
    time_conversion[timeii]=0;
  }
}

Neighbor_List::Neighbor_List(const Neighbor_List& copy):Value_List<int>(copy)
{
  
  neighbors=copy.neighbors;
  computed_times=copy.computed_times;
}


Neighbor_List Neighbor_List::operator=(const Neighbor_List& copy)
{
  if(this!=&copy)
  {

    Value_List::operator=(copy);
    neighbors=copy.neighbors;
    computed_times=copy.computed_times;
    
  }
  return *this;
}


Neighbor_List::~Neighbor_List()
{
}


Neighbor_List::Neighbor_List(System * sys):Value_List<int>(sys)
{
  neighbors.resize(n_times);
  for(int timeii=0;timeii<n_times;timeii++)
  {
    neighbors[timeii].resize(syst->show_n_trajectories());
  }
  computed_times.resize(n_times,false);
}


bool Neighbor_List::is_neighbor(int timeii, int trajii, Trajectory* trajcheck)const
{
  bool check = false;
  
  int neighborii;
  int n_neighbors=neighbors[timeii][trajii].size();
  if(!computed_times[timeii])
  {
    cout<<"\Error: neighbor list has not been computed for time " << timeii << ".\n";
    exit(0);
  }
  else if(!(included[timeii])(trajii))
  {
    cout<<"\nWarning: trajectory " << trajii << " in trajectory list is not included in the neighbor list.";
    check=false;
  }
  else
  {
    for(neighborii=0;neighborii<n_neighbors;neighborii++)
    {
      check=check||(trajcheck==neighbors[timeii][trajii][neighborii]);
    }
  }
  
  return check;
}


void Neighbor_List::update_size()
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
    if(neighbors[timeii].size()<syst->show_n_trajectories())
    {
      neighbors[timeii].resize(syst->show_n_trajectories());
    }
  }
}




vector<Trajectory*> Neighbor_List::persistent_neighbors(int trajii, int time1)
{
  
  vector<Trajectory*> temp;
  
  if(computed_times[time1]&&included[time1](trajii))
  {
    temp=neighbors[time1][trajii];
  }
  
  return temp;
}



vector<Trajectory*> persistent_neighbors(int trajii, int time1, int time2)
{
  
}