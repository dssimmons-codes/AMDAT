/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/


#include "neighbor_list.h"


using namespace std;




Neighbor_List::Neighbor_List()
{
  syst=0;
  included = new Boolean_List[1];
}

Neighbor_List::Neighbor_List(const Neighbor_List& copy)
{
  syst=copy.syst;
  n_times=copy.n_times;
  included = new Boolean_List [n_times];
  for(int timeii=0; timeii<n_times; timeii++)
  {
    included[timeii]=copy.included[timeii];
  }
  neighbors=copy.neighbors;
}


Neighbor_List Neighbor_List::operator=(const Neighbor_List&)
{
  if(this!=&copy)
  {
    delete [] included;
    syst=copy.syst;
    n_times=copy.n_times;
    included = new Boolean_List [n_times];
    for(int timeii=0; timeii<n_times; timeii++)
    {
      included[timeii]=copy.included[timeii];
    }
    neighbors=copy.neighbors;
  }
  return *this;
}


Neighbor_List::~Neighbor_List()
{
   delete [] included;
}


bool Neighbor_List::is_neighbor(int timeii, int trajii, Trajectory* trajcheck)const
{
  bool check = false;
  
  int neighborii;
  int n_neighbors=neighbors[timeii][trajii].size();
  
  for(neighborii=0;neighborii<n_neighbors;neighborii++)
  {
    check=check||(trajcheck==neighbors[timeii][trajii][neighborii]);
  }
  
  return check;
}