// Copyright (C) 2010-2025 David S. Simmons

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
  max_trajectories=0;
  
  multibodyID=new int[1];
  multibodyID[1]=-1;
  imageindex = new Coordinate [1];
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
    max_trajectories=copy.max_trajectories;
    
    if(system!=0)
    {
      time_conversion=new int [system->show_n_timesteps()];
      for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
      {
	time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
      }
      max_trajectories=system->show_n_trajectories();
  
      multibodyID = new int [max_trajectories];
      imageindex = new Coordinate [max_trajectories];
  
      for(int trajii=0;trajii<max_trajectories;trajii++)
      {
	multibodyID[trajii]=-1;
      }
    }
    else
    {
      time_conversion = new int [1];
      multibodyID = new int [1];
      multibodyID[0]=-1;
      imageindex = new Coordinate [1];
    }
  }
  return *this;
}

Dynamic_Cluster_Multibodies::Dynamic_Cluster_Multibodies(System*syst, int tgap)
{
  system=syst;
  timegap=tgap;
  time_conversion=new int [system->show_n_timesteps()];
  for(int timeii=0;timeii<system->show_n_timesteps()-1;timeii++)
  {
    time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
  }
  time_conversion[system->show_n_timesteps()-1]=time_conversion[system->show_n_timesteps()-2];
  
  n_times=system->show_n_exponentials();
  multibodyID = new int [1];
  multibodyID[0]=-1;
  imageindex = new Coordinate [1];
  
  max_trajectories=system->show_n_trajectories();
  
  multibodyID = new int [max_trajectories];
  imageindex = new Coordinate [max_trajectories];
  
  for(int trajii=0;trajii<max_trajectories;trajii++)
  {
    multibodyID[trajii]=-1;
  }
  

}


void Dynamic_Cluster_Multibodies::analyze(Trajectory_List*t_list)
{
	trajectory_list=t_list;
	multibodies.resize(0);
	system->displacement_list(this, timegap);
	postprocess_list();
}


 void Dynamic_Cluster_Multibodies::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	trajectories_considered+=(trajectory_list[0]).show_n_trajectories(thisii);
	n_trajectories = trajectory_list->show_n_trajectories(thisii);
	int currenttime = int(float(thisii)/float(system->show_n_exponential_steps()));
	multibodies.resize(multibodies.size()+1);
	delete [] multibodyID;
	delete [] imageindex;
	
	//assign memory to track array of string associated with each particle
	multibodyID = new int [max_trajectories];
	imageindex = new Coordinate [max_trajectories];
	
	for(int trajii=0;trajii<max_trajectories;trajii++)
	{
	  multibodyID[trajii]=-1;
	}
	
	nn_bodies=0;
	
	//loop over all trajectories, adding, growing, or concatenating strings associated with each
	multibody_validity.clear();
	(trajectory_list[0]).listloop(this,timegapii,thisii,nextii);
	//eliminated discarded strings
	for(int multibodyii=multibodies[currenttime].size()-1;multibodyii>=0;multibodyii--)
	{
	  if(!multibody_validity[multibodyii])
	  {
	    (multibodies[currenttime]).erase(multibodies[currenttime].begin()+multibodyii);
	  }
	}
}



void Dynamic_Cluster_Multibodies::listkernel(Trajectory* trajectory1, int timegapii, int thisii, int nextii)
{
  //thisii is the system time index and should be used for calls to objects that use a system time index. current time is the internal time index and should be used for calls to internal data members
  
    int trajii;
    bool clustered;
    Trajectory * trajectory2;
    int trajectory1ID, trajectory2ID;
    int currenttime = int(float(thisii)/float(system->show_n_exponential_steps()));
    trajectory1ID=trajectory1->show_trajectory_ID();
    Coordinate imageoffset, imagemodifier;
    
    
    for(trajii=0;trajii<n_trajectories;trajii++)
    {
      trajectory2 = (*trajectory_list)(thisii, trajii);
      trajectory2ID=trajectory2->show_trajectory_ID();
      if(trajectory2!=trajectory1)
      {
	
	clustered=clustered_check(trajectory1, trajectory2, thisii, nextii);
	if(clustered)		//if the two trajectories are identified to be in the same cluster
	{
	  imageoffset=get_imageoffset(trajectory1, trajectory2, thisii, nextii);
	  if(multibodyID[trajectory1ID]==-1 && multibodyID[trajectory2ID]==-1)	//if neither atom is in a multibody
	  {
	    //cout << "\n\nnew\t"<<trajectory1ID<<"\t"<<trajectory2ID<<"\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	    /*create new multibody*/
	    imageindex[trajectory1ID].set(0,0,0);
	    imageindex[trajectory2ID]=imageoffset;
	    multibodies[currenttime].emplace_back(system,thisii);
	    multibody_validity.push_back(true);
	    multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory1,imageindex[trajectory1ID]);
	    multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory2,imageindex[trajectory2ID]);
	    multibodyID[trajectory1ID] = multibodies[currenttime].size()-1;
	    multibodyID[trajectory2ID] = multibodies[currenttime].size()-1;
	    nn_bodies+=2;
	    //cout << "\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	  }
	  else if(multibodyID[trajectory1ID]==-1 && multibodyID[trajectory2ID]>=0)	//if atom 2 is in a multibody but atom 1 is not
	  {
	    //cout << "\n\n1to2\t"<<trajectory1ID<<"\t"<<trajectory2ID<<"\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	    imageindex[trajectory1ID]=imageindex[trajectory2ID]-imageoffset;
	    multibodies[currenttime][multibodyID[trajectory2ID]].add_body(trajectory1,imageindex[trajectory1ID]);
	    multibodyID[trajectory1ID] =  multibodyID[trajectory2ID];
	    nn_bodies+=1;
	    //cout << "\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	  }
	  else if(multibodyID[trajectory1ID]>=0 && multibodyID[trajectory2ID]==-1)	//if atom 1 is in a multibody but atom 2 is not
	  {
	    //cout << "\n\n2to1\t"<<trajectory1ID<<"\t"<<trajectory2ID<<"\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	    imageindex[trajectory2ID]=imageindex[trajectory1ID]+imageoffset;
	    multibodies[currenttime][multibodyID[trajectory1ID]].add_body(trajectory2,imageindex[trajectory2ID]);
	    multibodyID[trajectory2ID] =  multibodyID[trajectory1ID];
	    nn_bodies+=1;
	    //cout << "\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	  }
	  else if(multibodyID[trajectory1ID]>=0 && multibodyID[trajectory2ID]>=0 && multibodyID[trajectory1ID]!=multibodyID[trajectory2ID])		//if both atoms are in multibodies but not the same multibody
	  {
	    //cout << "\n\nmerged\t"<<trajectory1ID<<"\t"<<trajectory2ID<<"\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	    imagemodifier=imageindex[trajectory1ID]+imageoffset-imageindex[trajectory2ID];
	    multibodies[currenttime][multibodyID[trajectory1ID]].absorb_multibody(multibodies[currenttime][multibodyID[trajectory2ID]],imagemodifier);	//merge second multibody into first, with appropriate correction applied for relative image indices for bodies in second multibody
	    for(int bodyii=0;bodyii<(multibodies[currenttime][multibodyID[trajectory2ID]]).show_n_bodies();bodyii++)	//loop over bodies in second multibody
	    {
	      imageindex[(multibodies[currenttime][multibodyID[trajectory2ID]]).show_body_ID(bodyii)]+=imagemodifier;	//update relative image flags of bodies in second multibody to accord with merge
	    }
	    //multibodies[currenttime][multibodyID[trajectory2ID]].clear();
	    multibody_validity[multibodyID[trajectory2ID]]=false;	//note that second multibody should be erased
	    mass_switch_ID(multibodyID[trajectory2ID],multibodyID[trajectory1ID]);	//switch multibody IDs for trajectories in second multibody to that of first multibody to accord with merge
	    //cout << "\t"<<multibodyID[trajectory1ID]<<"\t"<<multibodyID[trajectory2ID];
	    
	  }
	}
      }
    }
    


    if(multibodyID[trajectory1ID]==-1)	//create single-body multibody for multibody 1 if this multibody does not belong to a cluster (allows for multibodies of size 1, such that all trajectories in the list will be in one of the multibodies - this can later be removed with thresholding operations
    {
      //cout << "\n\nsinglet\t"<<trajectory1ID<<"\t"<<multibodyID[trajectory1ID];
      imageindex[trajectory1ID].set(0,0,0);
      multibodies[currenttime].emplace_back(system,thisii);
      multibody_validity.push_back(true);
      multibodies[currenttime][multibodies[currenttime].size()-1].add_body(trajectory1,imageindex[trajectory1ID]);
      multibodyID[trajectory1ID]= multibodies[currenttime].size()-1;
      nn_bodies+=1;
    }


}

int Dynamic_Cluster_Multibodies::mass_switch_ID(int oldID, int newID)
{
  int changed=0;

  for(int trajii=0;trajii<max_trajectories;trajii++)
  {
    if(multibodyID[trajii]== oldID)
    {
      multibodyID[trajii]=newID;
      changed++;
    }
  }
  
  //cout<<"\tnchanged "<<changed;
  return changed;
}

