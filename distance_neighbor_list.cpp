/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Distance_Neighbor_List class - builds neighbor list from distance thresholding*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "distance_neighbor_list.h"

using namespace std;



Distance_Neighbor_List::Distance_Neighbor_List():Neighbor_List()
{
  threshold=0;
  system=0;
  n_times=1;
  n_atomtypes=1;
  sigmatrix=new float* [1];
  sigmatrix[0] = new float [1];
}



Distance_Neighbor_List::Distance_Neighbor_List(const Distance_Neighbor_List& copy):Neighbor_List(copy)
{
  int typeii;
  
  threshold=copy.threshold;
  system=copy.system;
  n_atomtypes=copy.n_atomtypes;
  sigmatrix = new float* [n_atomtypes];
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    sigmatrix[typeii] = new float [n_atomtypes];
  }
    
 
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    for(int type2ii=0;type2ii<n_atomtypes;type2ii++)
    {
	sigmatrix[typeii][type2ii]=copy.sigmatrix[typeii][type2ii];
    }
  }
}




Distance_Neighbor_List Distance_Neighbor_List::operator=(const Distance_Neighbor_List& copy)
{
  int typeii;
  
  if(this!=&copy)
  {
    delete [] included;
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    
    delete [] sigmatrix;
    
    n_atomtypes=copy.n_atomtypes;
    threshold=copy.threshold;
    system=syst=copy.system;
    n_times=copy.n_times;
    included = new Boolean_List [n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
      included[timeii]=copy.included[timeii];
    }
  
    sigmatrix = new float* [n_atomtypes];
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      sigmatrix[typeii] = new float [n_atomtypes];
    }
    
    
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      for(int type2ii=0;type2ii<n_atomtypes;type2ii++)
      {
	  sigmatrix[typeii][type2ii]=copy.sigmatrix[typeii][type2ii];
      }
    }
    
  }
  neighbors=copy.neighbors;
  computed_times=copy.computed_times;
  

  return *this;
}


Distance_Neighbor_List::~Distance_Neighbor_List()
{
    for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    delete [] sigmatrix;
}


Distance_Neighbor_List::Distance_Neighbor_List(System* sys, float thresh, string sigmatrixname, int firsttime, int lasttime):Neighbor_List(sys)
{
  int startoffset=firsttime;
  int endoffset;
  
  threshold=thresh;
  system=sys;
  n_times=system->show_n_timesteps();
  included = new Boolean_List [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(sys);
  }
  allocate_sig_matrix(sigmatrixname);
  computed_times.clear();
  
  if(firsttime==-1)
  {
    computed_times.clear();
    computed_times.resize(n_times,true);
  }
  else
  {
    if(lasttime==-1)
    {
      endoffset=firsttime;
    }
    else
    {
      endoffset=lasttime;
    }
  
    int n_blocks = system->show_n_exponentials();
    int blocksize = system->show_n_exponential_steps();
    int blockstart;
  
    if(firsttime>blocksize||lasttime>blocksize)
    {
      cout<<"\nError: time offsets selected for neighbor list construction must fit within the block.\n";
      exit(0);
    }
  
    for(int blockii=0;blockii<n_blocks;blockii++)
    {
      blockstart = blockii*blocksize;
      computed_times[blockstart]=true;
      for(int offsetii=startoffset; offsetii<=endoffset;offsetii++)
      {
	computed_times[blockstart+offsetii]=true;
      }
    }
  }
  
}


void Distance_Neighbor_List::timekernel2(int timeii)
{
  if(computed_times[timeii])
  {
   trajectory_list->listloop(this,0, timeii, 0);
  }
}

void Distance_Neighbor_List::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  int trajectory1ID;
  trajectory1ID=current_trajectory->show_trajectory_ID();
  (included[thisii])(trajectory1ID,1);
  
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}

void Distance_Neighbor_List::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  int trajectory1ID, trajectory2ID, trajtype1, trajtype2;
  trajectory1ID=traj1->show_trajectory_ID();
  trajectory2ID=traj2->show_trajectory_ID();
  
  trajtype1 = traj1->show_type()-1;
  trajtype2 = traj2->show_type()-1;
  
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    if(distance<(threshold*sigmatrix[trajtype1][trajtype2]))
    {
      neighbors[thisii][trajectory1ID].push_back(traj2);
    }
  }
}



/*allocate matrix of particle sizes and assign values*/
void Distance_Neighbor_List::allocate_sig_matrix(string sig_file)
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