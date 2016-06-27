

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibodies.h"
#include "version.h"









String_Multibodies::String_Multibodies()
{
  system=0;
  timegap=-1;
  threshold=-1;
  n_times=0;
  time_conversion=new int [1];
  n_atomtypes=0;
  
}


String_Multibodies(const String_Multibodies& copy):Provisional_Multibodies(copy)
{
    system=copy.system;
    timegap=copy.timegap;
    threshold=thresh;
    n_times=copy.n_times;
    n_atomtypes=copy.n_atomtypes;
    
    if(system!=0)
    {
      time_conversion=new int [system->show_n_timesteps()];
      for(timeii=0;timeii<system->show_n_timesteps();timeii++)
      {
	time_conversion[timeii]=int(float(timeii-sys->show_frt())/float(system->show_n_exponential_steps()));
      }
    }
    else
    {
      time_conversion = new int [1];
    }
    
    basename=copy.basetime;
    for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      for(int type2ii=0;type2ii<n_atomtypes;type2ii++)
      {
	sigmatrix[typeii][type2ii]=copy.sigmatrix[typeii][type2ii];
      }
    }
}

~String_Multibodies():~Provisional_Multibodies();
{
  delete [] stringID;
  
}


String_Multibodies operator=(const String_Multibodies& copy)
{
  if(this!=&copy)
  {
    system=copy.system;
    timegap=copy.timegap;
    threshold=thresh;
    n_times=copy.n_times;
    n_atomtypes=copy.n_atomtypes;
    
    if(system!=0)
    {
      time_conversion=new int [system->show_n_timesteps()];
      for(timeii=0;timeii<system->show_n_timesteps();timeii++)
      {
	time_conversion[timeii]=int(float(timeii-sys->show_frt())/float(system->show_n_exponential_steps()));
      }
    }
    else
    {
      time_conversion = new int [1];
    }
    
    basename=copy.basetime;
    for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      for(int type2ii=0;type2ii<n_atomtypes;type2ii++)
      {
	sigmatrix[typeii][type2ii]=copy.sigmatrix[typeii][type2ii];
      }
    }
  }
  
  return *this;
}


String_Multibodies(System * syst, int tgap, float thresh, string sigmatrixname)
{
  system=syst;
  timegap=tgap;
  threshold=thresh;
  n_times=system->show_n_exponentials();
  time_conversion=new int [system->show_n_timesteps()];
  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=int(float(timeii-sys->show_frt())/float(system->show_n_exponential_steps()));
  }
  basename="";
  allocate_sig_matrix(sigmatrixname);
}




/*allocate matrix of particle sizes and assign values*/
void String_Multibodies::allocate_sig_matrix(string sig_file)
{
  Tokenize tokenize;
  
    string line;
    line = "";
    int sig_tokens=0;
    string * sig_ARGS;
    int matsize, lineii, argii, type1index, type2index;
    int * type_index;
    
    n_atomtypes = system->show_n_expanded_atomtypes();
    sig_ARGS =new string [n_atomtypes+1];

    string * species_name;
    species_name = new string [0];

    ifstream file(sig_file.c_str());

    
    sigmatrix=new float* [n_atomtypes];
	for(lineii=0;lineii<n_atomtypes;lineii++)
	{
	  sigmatrix[lineii]=new float[n_atomtypes];
	  for(argii=0;argii<n_atomtypes;argii++)
	  {
	    sigmatrix[lineii][argii]=0;
	  }
	}
    
    
    if (file.is_open())
    {
       
        //get first line of matrix
        getline (file,line);
        sig_tokens = tokenize(line, sig_ARGS);
        matsize = sig_tokens-1;
	
	if(matsize>n_atomtypes)
	{
	  cout<<"\nError: matrix size is greater than number of trajectory types.\n";
	  exit(0);
	}
	
		
	type_index = new int [matsize];
	
	if(!system->atomtype_exists(sig_ARGS[0]))
	{
	  cout << "\nError: Trajectory type " << sig_ARGS[0] << " not found.\n";
	  exit(0);
	}
        type_index[0] = system->show_atomtype_index(sig_ARGS[0]);

	for(lineii=1;lineii<matsize;lineii++)
	{
	  if(file.eof())
	  {
	    cout<<"\nError: Sigma matrix is not square.\n";
	    exit(0);
	  }
	  getline (file,line);
	  sig_tokens = tokenize(line, sig_ARGS);
	  if(sig_tokens!=matsize+1)
	  {
	    cout<<"\nError: Sigma matrix is not square.\n";
	    exit(0);
	  }
	  

	  if(!system->atomtype_exists(sig_ARGS[0]))
	  {
	    cout << "\nError: Trajectory type " << sig_ARGS[0] << " not found.\n";
	    exit(0);
	  }
          type_index[0] = system->show_atomtype_index(sig_ARGS[0]);
	}
	
	file.clear();
	file.seekg(0,ios::beg);
	
	for(lineii=0;lineii<matsize;lineii++)
	{
	  getline (file,line);
	  sig_tokens = tokenize(line, sig_ARGS);
	  for(argii=1; argii<= matsize; argii++)
	  {
	    sigmatrix[type_index[lineii]][type_index[argii-1]]=atof(sig_ARGS[argii].c_str());
	  }
	}
    }
    else
    {
        cout << "\nError: sigma data file not opened succesfully.\n";
        exit(1);
    }
    file.close();


}



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
	multibodies.push_back();
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


