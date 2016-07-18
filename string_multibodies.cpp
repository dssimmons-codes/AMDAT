/*Methods for String_Multibodies class - Identifies particles participating in stringlike cooperative rearrangements and converts them to multibodies*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "multibody.h"
#include "version.h"
#include "string_multibodies.h"
#include "system.h"


using namespace std;



String_Multibodies::String_Multibodies():Dynamic_Cluster_Multibodies()
{
  threshold=0;
  n_atomtypes=1;
  sigmatrix=new float*[1];
  sigmatrix[0]=new float[1];
  
}


String_Multibodies::String_Multibodies(const String_Multibodies& copy):Dynamic_Cluster_Multibodies(copy)
{
  int typeii, type2ii;
  
    threshold=copy.threshold;
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

String_Multibodies::~String_Multibodies()
{
  for(int typeii=0;typeii<n_atomtypes;typeii++)
    {
      delete [] sigmatrix[typeii];
    }
    delete [] sigmatrix;
  
}


String_Multibodies String_Multibodies::operator=(const String_Multibodies& copy)
{
  if(this!=&copy)
  {
    int typeii, type2ii;
    
    system=copy.system;
    timegap=copy.timegap;
    threshold=copy.threshold;
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


String_Multibodies::String_Multibodies(System * syst, int tgap, float thresh, string sigmatrixname):Dynamic_Cluster_Multibodies(syst,tgap)
{
  system=syst;
  threshold=thresh;
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
    int lineii, argii, type1index, type2index;
        int * type_index;
    int matsize;

    n_atomtypes = system->show_n_expanded_atomtypes();
    sig_ARGS =new string [n_atomtypes+1];

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
       /*Learn atomtypes and check for squareness*/
       
       
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
          type_index[lineii] = system->show_atomtype_index(sig_ARGS[0]);
	}
	
	/*Now that atom types are known and squareness is checked, read in sigma values*/
	
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
    file.close();
    }
    else
    {
        cout << "\nError: sigma data file not opened succesfully.\n";
        exit(1);
    }
    
    cout<<"\nSigma matrix used: ";
    for (int indii=0;indii<n_atomtypes;indii++)
    {
      cout <<"\n" << indii << "\t";
      for (lineii=0;lineii<n_atomtypes;lineii++)
      {
	cout<<sigmatrix[indii][lineii]<<"\t";
      }
    }
    cout <<"\n";


}


bool String_Multibodies::clustered_check(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  bool check;
  int trajtype1, trajtype2;
  float distance;
  
  trajtype1 = trajectory1->show_type()-1;
  trajtype2 = trajectory2->show_type()-1;
  
  distance = (trajectory1->show_coordinate(thisii)-trajectory2->show_coordinate(nextii)).length_unwrapped(system->size(thisii));
  check= (distance<(threshold*sigmatrix[trajtype1][trajtype2]));
  return check;
}



Coordinate String_Multibodies::get_imageoffset(Trajectory* trajectory1, Trajectory* trajectory2, int thisii, int nextii)
{
  return (trajectory1->show_coordinate(thisii)).closest_image(trajectory2->show_coordinate(thisii),system->size(thisii));
}

