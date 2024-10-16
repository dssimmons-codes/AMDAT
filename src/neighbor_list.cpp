/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Neighbor_List class - stores a time-dependent neighbor list for a set of trajectories*/
/*Written by David S. Simmons*/


#include "neighbor_list.h"
#include "version.h"
#include "tokenize.h"

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




vector<Trajectory*> Neighbor_List::show_neighbors(int trajii, int time1)const
{
  
  vector<Trajectory*> temp;
  
  if(computed_times[time1]&&included[time1](trajii))
  {
    temp=neighbors[time1][trajii];
  }
  
  return temp;
}


vector<Trajectory*> Neighbor_List::persistent_neighbors(int trajii, int time1, int time2)const
{
  
  vector<Trajectory*> temp;
  int n1ii, n2ii;
  bool check;
  
  if(computed_times[time1]&&included[time1](trajii))
  {
    for(n1ii=0;n1ii<neighbors[time1][trajii].size();n1ii++)
    {
      check=false;
      for(n2ii=0;n2ii<neighbors[time2][trajii].size();n2ii++)
      {
	if(neighbors[time2][trajii][n2ii]==neighbors[time1][trajii][n1ii])
	{
	  check=true;
	}
      }
      if(check)
      {
	temp.push_back(neighbors[time1][trajii][n1ii]);
      }
    }
  }
  
  return temp;
}


int Neighbor_List::n_persistent_neighbors(int trajii, int time1, int time2)const
{
  
  int n_count=0;
  int n1ii, n2ii;
  bool check;
  
  if(computed_times[time1]&&included[time1](trajii))
  {
    for(n1ii=0;n1ii<neighbors[time1][trajii].size();n1ii++)
    {
      check=false;
      for(n2ii=0;n2ii<neighbors[time2][trajii].size();n2ii++)
      {
        if(neighbors[time2][trajii][n2ii]==neighbors[time1][trajii][n1ii])
        {
          check=true;
        }
      }
      if(check)
      {
        n_count++;
      }
    }
  }
  
  return n_count;
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

void Neighbor_List::read_bonds(string filename)
{
  string line;
  vector <string> args;
  int n_bonds_listed;
  int trajid1, trajid2;
  Tokenize tokenize;

  ifstream filexyz(filename.c_str());
  ifstream * fileobject = &filexyz;

  for(int timestepii=0; timestepii<n_times; timestepii++)
  {
    /*read in header information from  file*/
    line = "";
    getline(*fileobject,line);		//read in "ITEM: TIMESTEP" line
    getline(*fileobject,line);		//read in timestep line
    getline(*fileobject,line);		//read in "ITEM: NUMBER OF ENTRIES" line
    line = "";
    getline(*fileobject,line);		//read in number of entries
    args = tokenize(line);
    n_bonds_listed = atoi(args[0].c_str());

    getline(*fileobject,line);		//read in box bounds header
    getline(*fileobject,line);		//read in box bounds
    getline(*fileobject,line);		//read in box bounds
    getline(*fileobject,line);		//read in box bounds
    getline(*fileobject,line);		//read in header

     for(int bondii=0; bondii<n_bonds_listed; bondii++)
     {
       line = "";
       getline(*fileobject,line);   //read in next bond
       args = tokenize(line);
       trajid1 = atoi(args[0].c_str())-1;   //get trajectory index from bond line; assume trajectory index is lammps atom number - 1. This will work provided that the trajectory file itself was sorted in the order of atom IDs. Otherwise, this method will nominally work for bond statistics but will not correctly associated bonds with their correct atoms.
       trajid2 = atoi(args[1].c_str())-1;  //same for second atom in bond

       /*make sure that both bonded atoms are in the boolean list of atoms with neighbors*/
       (included[timestepii])(trajid1,1);
       (included[timestepii])(trajid2,1);

       /*add the two atoms to each others' neighbor lists*/
       neighbors[timestepii][trajid1].push_back(syst->show_trajectory(trajid2));
       neighbors[timestepii][trajid2].push_back(syst->show_trajectory(trajid1));
     }
  }

  values.resize(neighbors.size());
  int timeii,trajii;

  for(timeii=0;timeii<neighbors.size();timeii++)
  {
    if(computed_times[timeii])
    {
      values[timeii].resize(neighbors[timeii].size());
      for(trajii=0;trajii<neighbors[timeii].size();trajii++)
      {
        if(included[timeii](trajii))
        {
          values[timeii][trajii]=neighbors[timeii][trajii].size();
        }
      }
    }
  }
}
