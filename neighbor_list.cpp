


using namespace std;




Neighbor_List()
{
  syst=0;
  included = new Boolean_List[1];
}

Neighbor_List(const Neighbor_List& copy)
{
  syst=copy.syst;
  if(syst!=0)
  {
    included=system->
  }
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