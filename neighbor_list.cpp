/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/


#include "neighbor_list.h"
#include "version.h"


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

Neighbor_List::Neighbor_List(const Neighbor_List& copy):Value_List<float>(copy)
{
  
  neighbors=copy.neighbors;
  computed_times=copy.computed_times;
}


Neighbor_List Neighbor_List::operator=(const Neighbor_List& copy)
{
  if(this!=&copy)
  {
    Value_List<float>::operator=(copy);
    neighbors=copy.neighbors;
    computed_times=copy.computed_times;
    
  }
  return *this;
}


Neighbor_List::~Neighbor_List()
{
}


Neighbor_List::Neighbor_List(System * sys):Value_List<float>(sys)
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




vector<Trajectory*> Neighbor_List::persistent_neighbors(int trajii, int time1)
{
  
  vector<Trajectory*> temp;
  
  if(computed_times[time1]&&included[time1](trajii))
  {
    temp=neighbors[time1][trajii];
  }
  
  return temp;
}



vector<Trajectory*> Neighbor_List::persistent_neighbors(int trajii, int time1, int time2)
{
  
}


void Neighbor_List::write_statistics(string filename, int n_moments)const
{
  int timeii, trajii, binii;
  int maximum = max();
  int n_values=0;
  
  int sizeii, momentii;
  float* moments;
  moments=new float [n_moments];
  
  cout << "\nWriting value dist and statistics to file.";
  ofstream output(filename.c_str());
  output << "Value list statistics created by AMDAT v." << VERSION << "\n";
   
  vector<int> dist;
  
  dist.resize(maximum+1,0);
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    if(computed_times[timeii])
    {
      n_values+=included->show_n_included();
      for(trajii=0;trajii<values[timeii].size();trajii++)
      {
	if(included[timeii](trajii))
	{
	  dist[values[timeii][trajii]]++;
	}
      }
    }
  }
  
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    moments[momentii]=0;
  }
  
  
  for(binii=0;binii<dist.size();binii++)
  {
    moments[0]+=dist[binii];
    for(momentii=1;momentii<n_moments;momentii++)
    {
      moments[momentii]+=pow(float(binii),float(momentii))*float(dist[binii])/float(n_values);
    }
  }
  
  output << "Value\t";
  for(binii=0;binii<dist.size();binii++)
  {
    output << binii << "\t";
  }
  
    output<<"\nFrequency\t";
  
  for(binii=0;binii<dist.size();binii++)
  {
    
    output<<float(dist[binii])/float(n_values)<<"\t";
  }
  
  output<<"\n\nMoment\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<momentii<<"\t";
  }
  
  output << "\nValue\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<moments[momentii]<<"\t";
  }
  
  output << "\n\nTotal_Values\t" << n_values;
  
}
