/*Methods for Relative_Displacement_Strings class - Identifies particles participating in stringlike cooperative rearrangements and converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibody.h"
#include "version.h"
#include "relative_displacement_strings.h"
#include "system.h"


using namespace std;



Relative_Displacement_Strings::Relative_Displacement_Strings():Dynamic_Cluster_Multibodies()
{
  replacement_threshold=0;
  neighbor_threshold=0;
  n_atomtypes=1;
  sigmatrix=new float*[1];
  sigmatrix[0]=new float[1];
  
}


Relative_Displacement_Strings::Relative_Displacement_Strings(const Relative_Displacement_Strings& copy):Dynamic_Cluster_Multibodies(copy)
{
  int typeii, type2ii;
  
    neighbor_threshold=copy.neighbor_threshold;
    replacement_threshold=copy.replacement_threshold;
    n_times=copy.n_times;
    
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    delete [] sigmatrix;
    
    n_atomtypes=copy.n_atomtypes;
    sigmatrix = new float* [n_atomtypes];
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      sigmatrix[typeii] = new float [n_atomtypes];
    }
    
    basename=copy.basename;
    for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      for(int type2ii=0;type2ii<n_atomtypes;type2ii++)
      {
	sigmatrix[typeii][type2ii]=copy.sigmatrix[typeii][type2ii];
      }
    }
}

Relative_Displacement_Strings::~Relative_Displacement_Strings()
{
  for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    delete [] sigmatrix;
  
}


Relative_Displacement_Strings Relative_Displacement_Strings::operator=(const Relative_Displacement_Strings& copy)
{
  if(this!=&copy)
  {
    int typeii, type2ii;
    
    system=copy.system;
    timegap=copy.timegap;
    neighbor_threshold=copy.neighbor_threshold;
    replacement_threshold=copy.replacement_threshold;
    n_times=copy.n_times;
    
    for(typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    delete [] sigmatrix;
    
    n_atomtypes=copy.n_atomtypes;
    
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
    
    basename=copy.basename;
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


Relative_Displacement_Strings::Relative_Displacement_Strings(System * syst, int tgap, float n_thresh,float r_thresh, string sigmatrixname):Dynamic_Cluster_Multibodies(syst,tgap)
{
  system=syst;
  neighbor_threshold=n_thresh;
  replacement_threshold=n_thresh;
  allocate_sig_matrix(sigmatrixname);
}




/*allocate matrix of particle sizes and assign values*/
void Relative_Displacement_Strings::allocate_sig_matrix(string sig_file)
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


bool Relative_Displacement_Strings::clustered_check(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  bool check;
  int trajtype1, trajtype2;
  float initial_separation,distance;
  
  trajtype1 = trajectory1->show_type()-1;
  trajtype2 = trajectory2->show_type()-1;
  
  initial_separation=(trajectory1->show_coordinate(thisii)-trajectory2->show_coordinate(thisii)).length_unwrapped(system->size(thisii));
  
  distance = (trajectory1->show_coordinate(thisii)-trajectory2->show_coordinate(nextii)).length_unwrapped(system->size(thisii));
  
  check= (initial_separation<(neighbor_threshold*sigmatrix[trajtype1][trajtype2]))&&(distance<(replacement_threshold*initial_separation));
  return check;
}



Coordinate Relative_Displacement_Strings::get_imageoffset(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  return (trajectory1->show_coordinate(thisii)).closest_image(trajectory2->show_coordinate(thisii),system->size(thisii));
}

