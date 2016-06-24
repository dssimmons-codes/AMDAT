

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibodies.h"
#include "version.h"












 void String_Multibodies::analyze(Trajectory_List*t_list)
{
	trajectory_list=t_list;
	system->displacement_list(this, timegap);
	postprocess_list();
}



 void String_Multibodies::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	trajectories_considered+=(trajectory_list[0]).show_n_trajectories(thisii);
	trajindex=0;
	n_trajectories = trajectory_list->show_n_trajectories(thisii);
	int currenttime = int(float(thisii)/float(system->show_n_exponential_steps()));
	multibodies.push(back);
	delete [] stringID;
	
	//assign memory to track array of string associated with each particle
	stringID = new int [n_trajectories];
	for(int trajii=0;trajii<n_trajectories;trajii++)
	{
	  stringID[trajii]=-1;
	}
	trajindex=-1;
	
	//loop over all trajectories, adding, growing, or concatenating strings associated with each
	(trajectory_list[0]).listloop(this,thistime);
	
	//eliminated discarded strings
	for(int stringii=multibodies[currenttime].size()-1;stringii>=0;stringii--)
	{
	  if(string_validity[stringii]=0)
	  {
	    multibodies[currenttime][stringii].erase();
	  }
	}
}



void String_Multibodies::listkernel(Trajectory* trajectory1, int timegapii, int thisii, int nextii)
{
    int trajii;
    Trajectory * trajectory2;
    
    int currenttime = int(float(thistime)/float(system->show_n_exponential_steps()));
    
    trajindex++;
    
    for(trajii=0;trajii<n_trajectories;trajii++)
    {
      trajectory2 = (*trajectory_list)(thistime, trajii);
      if(trajectory2!=trajectory1)
      {
	
	
	distance = (trajectory1->show_coordinate(thisii)-trajectory2->show_coordinate(nextii)).length_unwrapped(system->size(thistime));
	if(distance<threshold*sigmatrix[trajtype1][trajtype2])
	{
	  if(stringID[currenttime][trajectory1ID]==-1 && stringID[currenttime][trajectory2ID]==-1)	//if neither atom is in a string
	  {
	    /*create new string*/
	    multibodies[currenttime].push_back(system);
	    string_validity.push_back(true);
	    multibodies[currenttime][multibodies[currenttime].size()].add_body(trajectory1);
	    multibodies[currenttime][multibodies[currenttime].size()].add_body(trajectory2);
	    stringID[trajindex] =  multibodies[currenttime].size()-1;
	    stringID[trajii] = multibodies[currenttime].size()-1;
	  }
	  else if(stringID[trajindex]==-1 && stringID[trajii]>=0)	//if atom 2 is in a string but atom 1 is not
	  {
	    multibodies[currenttime][stringID[trajii]].addbody(trajectory1);
	    stringID[trajindex] =  stringID[trajii];
	  }
	  else if(stringID[trajindex]>=0 && stringID[trajii]==-1)	//if atom 1 is in a string but atom 2 is not
	  {
	    multibodies[currenttime][stringID[trajindex]].addbody(trajectory2);
	    stringID[trajii] =  stringID[trajindex];
	  }
	  else if(stringID[trajindex]>=0 && stringID[trajii]>0)		//if both atoms are in strings but not the same string
	  {
	    multibodies[currenttime][stringID[trajindex]].absorb_multibody(multibodies[currenttime][stringID[trajii]]);
	    string_validity[stringID[trajii]]=false;
	    mass_switch_ID(stringID[trajii],stringID[trajindex]);
	  }
	}
	
      }
    }
}


int String_Multibodies mass_switch_ID(int oldID, int newID)
{
  int changed=0;
  for(int trajii=0;trajii<n_trajectories;trajii++)
  {
    if(stringID[trajii]== oldID)
    {
      stringID[trajii]==newID;
      changed++;
    }
  }
  return changed;
}

