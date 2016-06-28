/*AMDAT: Amorphous Molecular Dynamics Analysis Toolkit*/
/*Methods for Dynamic_Cluster_Multibodies class: builds multibodies based on any pairwise dynamic criterion*/
/*David S Simmons*/


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibody.h"
#include "version.h"
#include "dynamic_cluster_multibodies.h"
#include "system.h"


using namespace std;

Dynamic_Cluster_Multibodies::Dynamic_Cluster_Multibodies():Provisional_Multibodies(), Analysis()
{
  timegap = -1;
  trajectories_considered=0;
  multibodyID=new int [1];
}


Dynamic_Cluster_Multibodies::Dynamic_Cluster_Multibodies(const Dynamic_Cluster_Multibodies&copy):Provisional_Multibodies(copy), Analysis(copy)
{
    timegap=copy.timegap;
    trajectories_considered=copy.trajectories_considered;
    multibodyID=new int [1];
    
    if(system!=0)
    {
      time_conversion=new int [system->show_n_timesteps()];
      for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
      {
	time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
      }
    }
    else
    {
      time_conversion = new int [1];
    }
}


Dynamic_Cluster_Multibodies::~Dynamic_Cluster_Multibodies()
{
  delete [] multibodyID;
}


Dynamic_Cluster_Multibodies Dynamic_Cluster_Multibodies::operator=(const Dynamic_Cluster_Multibodies& copy)
{
  if(this!=&copy)
  {
    Provisional_Multibodies::operator=(copy);
    Analysis::operator=(copy);
    timegap=copy.timegap;
    trajectories_considered=copy.trajectories_considered;
    
    if(system!=0)
    {
      time_conversion=new int [system->show_n_timesteps()];
      for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
      {
	time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
      }
    }
    else
    {
      time_conversion = new int [1];
    }
  }
  return *this;
}

Dynamic_Cluster_Multibodies::Dynamic_Cluster_Multibodies(System*syst, int tgap)
{
  system=syst;
  timegap=tgap;
  time_conversion=new int [system->show_n_timesteps()];
  for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
  }
  n_times=system->show_n_exponentials();
  multibodyID = new int [1];
  multibodyID[0]=-1;
  

}


void Dynamic_Cluster_Multibodies::analyze(Trajectory_List*t_list)
{
	trajectory_list=t_list;
	system->displacement_list(this, timegap);
	postprocess_list();
}


 void Dynamic_Cluster_Multibodies::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	trajectories_considered+=(trajectory_list[0]).show_n_trajectories(thisii);
	trajindex=0;
	n_trajectories = trajectory_list->show_n_trajectories(thisii);
	max_trajectories=system->show_n_trajectories();
	int currenttime = int(float(thisii)/float(system->show_n_exponential_steps()));
	multibodies.resize(multibodies.size()+1);
	delete [] multibodyID;
	
	//assign memory to track array of string associated with each particle
	multibodyID = new int [max_trajectories];
	for(int trajii=0;trajii<max_trajectories;trajii++)
	{
	  multibodyID[trajii]=-1;
	}
	trajindex=-1;
	
	//loop over all trajectories, adding, growing, or concatenating strings associated with each
	(trajectory_list[0]).listloop(this,timegapii,thisii,nextii);
	
	//eliminated discarded strings
	for(int multibodyii=multibodies[currenttime].size()-1;multibodyii>=0;multibodyii--)
	{
	  if(string_validity[multibodyii]=0)
	  {
	    (multibodies[currenttime]).erase(multibodies[currenttime].begin()+multibodyii);
	  }
	}
}



void Dynamic_Cluster_Multibodies::listkernel(Trajectory* trajectory1, int timegapii, int thisii, int nextii)
{
    int trajii;
    bool clustered;
    Trajectory * trajectory2;
    int trajectory1ID, trajectory2ID;
    int currenttime = int(float(thisii)/float(system->show_n_exponential_steps()));
    trajectory1ID=trajectory1->show_trajectory_ID();
    
    
    trajindex++;
    
    for(trajii=0;trajii<n_trajectories;trajii++)
    {
      trajectory2 = (*trajectory_list)(thisii, trajii);
      trajectory2ID=trajectory2->show_trajectory_ID();
      if(trajectory2!=trajectory1)
      {
	
	clustered=clustered_check(trajectory1, trajectory2, thisii, nextii);
	if(clustered)		//if the two trajectories are identified to be in the same cluster
	{
	  if(multibodyID[trajectory1ID]==-1 && multibodyID[trajectory2ID]==-1)	//if neither atom is in a multibody
	  {
	    /*create new string*/
	    multibodies[currenttime].emplace_back(system);
	    string_validity.push_back(true);
	    multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory1);
	    multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory2);
	    multibodyID[trajindex] =  multibodies[currenttime].size()-1;
	    multibodyID[trajii] = multibodies[currenttime].size()-1;
	  }
	  else if(multibodyID[trajindex]==-1 && multibodyID[trajii]>=0)	//if atom 2 is in a string but atom 1 is not
	  {
	    multibodies[currenttime][multibodyID[trajii]].add_body(trajectory1);
	    multibodyID[trajindex] =  multibodyID[trajii];
	  }
	  else if(multibodyID[trajindex]>=0 && multibodyID[trajii]==-1)	//if atom 1 is in a string but atom 2 is not
	  {
	    multibodies[currenttime][multibodyID[trajindex]].add_body(trajectory2);
	    multibodyID[trajii] =  multibodyID[trajindex];
	  }
	  else if(multibodyID[trajindex]>=0 && multibodyID[trajii]>0)		//if both atoms are in strings but not the same string
	  {
	    multibodies[currenttime][multibodyID[trajindex]].absorb_multibody(multibodies[currenttime][multibodyID[trajii]]);	//merge second multibody into first
	    string_validity[multibodyID[trajii]]=false;	//note that second multibody should be erased
	    mass_switch_ID(multibodyID[trajii],multibodyID[trajindex]);	//switch multibody IDs for trajectories in second multibody to accord with merge
	  }
	}

      }
    }
    if(multibodyID[trajectory1ID]==-1)	//create single-body multibody for multibody 1 if this multibody does not belong to a cluster (allows for multibodies of size 1, such that all trajectories in the list will be in one of the multibodies - this can later be removed with thresholding operations
    {
      multibodies[currenttime].emplace_back(system);
      string_validity.push_back(true);
      multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory1);
      multibodyID[trajindex]= multibodies[currenttime].size()-1;
    }
}

int Dynamic_Cluster_Multibodies::mass_switch_ID(int oldID, int newID)
{
  int changed=0;
  for(int trajii=0;trajii<max_trajectories;trajii++)
  {
    if(multibodyID[trajii]== oldID)
    {
      multibodyID[trajii]==newID;
      changed++;
    }
  }
  return changed;
}

