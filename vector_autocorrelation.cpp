/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for Vector_Autocorrelation*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "version.h"
#include "vector_autocorrelation.h"
#include "tokenize.h"

using namespace std;


/*Default constructor*/

Vector_Autocorrelation::Vector_Autocorrelation()
{
  system=0;
  n_times=0;
  n_vectors = 0;
  correlation=new float [n_times];
  weighting=new int[n_times];

  vector_specieslist = new int [n_vectors];
  vector_type1list = new int [n_vectors];
  vector_type2list = new int [n_vectors];
  vector_index1list = new int [n_vectors];
  vector_index2list = new int [n_vectors];
}


/*----------------------------------------------------------------------------------------*/


/*Copy constructor*/


Vector_Autocorrelation::Vector_Autocorrelation(const Vector_Autocorrelation & copy)
{
  system=copy.system;
  n_times=copy.n_times;;

  timetable = system->displacement_times();

  correlation=new float [n_times];
  orientational_correlation=new float [n_times];
  weighting=new int[n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    orientational_correlation[timeii]=copy.orientational_correlation[timeii];
    correlation[timeii]=copy.correlation[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  vector_specieslist = new int [n_vectors];
  vector_type1list = new int [n_vectors];
  vector_type2list = new int [n_vectors];
  vector_index1list = new int [n_vectors];
  vector_index2list = new int [n_vectors];

  for (int vectorii=0;vectorii<n_vectors;vectorii++)
  {
    vector_specieslist[vectorii]=copy.vector_specieslist[vectorii];
    vector_type1list[vectorii]=copy.vector_type1list[vectorii];
    vector_index1list[vectorii]=copy.vector_index1list[vectorii];
    vector_type2list[vectorii]=copy.vector_type2list[vectorii];
    vector_index2list[vectorii]=copy.vector_index2list[vectorii];
  }
}


/*---------------------------------------------------------------------------------------*/


/*Destructor*/


Vector_Autocorrelation::~Vector_Autocorrelation()
{
  delete [] correlation;
  delete [] weighting;
  delete [] vector_specieslist;
  delete [] vector_type1list;
  delete [] vector_type2list;
  delete [] vector_index1list;
  delete [] vector_index2list;
}


/*---------------------------------------------------------------------------------------*/


/*Equality operator*/

Vector_Autocorrelation Vector_Autocorrelation::operator = (const Vector_Autocorrelation & copy)
{
  if(this!=&copy)
  {

    system=copy.system;
  n_times=copy.n_times;;

  timetable = system->displacement_times();

  delete [] correlation;
  delete [] weighting;
  delete [] vector_specieslist;
  delete [] vector_type1list;
  delete [] vector_type2list;
  delete [] vector_index1list;
  delete [] vector_index2list;

  correlation=new float [n_times];
  weighting=new int[n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii]=copy.correlation[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  vector_specieslist = new int [n_vectors];
  vector_type1list = new int [n_vectors];
  vector_type2list = new int [n_vectors];
  vector_index1list = new int [n_vectors];
  vector_index2list = new int [n_vectors];

  for (int vectorii=0;vectorii<n_vectors;vectorii++)
  {
    vector_specieslist[vectorii]=copy.vector_specieslist[vectorii];
    vector_type1list[vectorii]=copy.vector_type1list[vectorii];
    vector_index1list[vectorii]=copy.vector_index1list[vectorii];
    vector_type2list[vectorii]=copy.vector_type2list[vectorii];
    vector_index2list[vectorii]=copy.vector_index2list[vectorii];
  }

  }
  return *this;
}


/*---------------------------------------------------------------------------------------*/


/*Full constructor that actually calculates autocorrelation data*/


Vector_Autocorrelation::Vector_Autocorrelation(System* sys, string bondfilename)
{

  initialize(sys,bondfilename);

  calculate_mean_vector_length();

  system->displacement_list(this,bool(0));

  for(int timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii]/=(float(weighting[timeii])*mean_length*mean_length);
    orientational_correlation[timeii]/=(float(weighting[timeii]));
  }

}


/*---------------------------------------------------------------------------------------*/


/*Calculates autocorrelation data*/


void Vector_Autocorrelation::set(System* sys, string bondfilename)
{

  delete [] correlation;
  delete [] weighting;
  delete [] vector_specieslist;
  delete [] vector_type1list;
  delete [] vector_type2list;
  delete [] vector_index1list;
  delete [] vector_index2list;


  initialize(sys,bondfilename);

  calculate_mean_vector_length();

  system->displacement_list(this,bool(0));

  for(int timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii]/=(float(weighting[timeii])*mean_length*mean_length);
    orientational_correlation[timeii]/=(float(weighting[timeii]));
  }

}



void Vector_Autocorrelation::initialize(System* sys, string bondfilename)
{
  string line;
  int n_args;
  string args [128];
  ifstream input;

  system=sys;

  n_times=system->show_n_timegaps();
  timetable = system->displacement_times();

  correlation=new float [n_times];
  orientational_correlation=new float [n_times];
  weighting=new int[n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    orientational_correlation[timeii]=0;
    correlation[timeii]=0;
    weighting[timeii]=0;
  }


  input.open(bondfilename.c_str());
  if (!input.is_open()){cout << "Error opening file " << bondfilename << ".\n"; exit(1);}

  line = "";
  getline(input,line);
  n_args = tokenize(line, args);

  n_vectors = atoi(args[0].c_str());

  vector_specieslist = new int [n_vectors];
  vector_type1list = new int [n_vectors];
  vector_type2list = new int [n_vectors];
  vector_index1list = new int [n_vectors];
  vector_index2list = new int [n_vectors];

  for (int vectorii=0;vectorii<n_vectors;vectorii++)
  {
    line = "";
    getline(input,line);
    n_args = tokenize(line, args);

    if(n_args!=5)
    {
      cout << "Number of arguments in bond file line is incorrect. Format is <species> <atom 1 species> <atom 1 index> <atom 2 species> <atom 2 index>.\n";
      exit(0);
    }
    else
    {
      vector_specieslist[vectorii]=system->show_species_index(args[0]);
      vector_type1list[vectorii]=system->show_atomtype_index(args[1]);
      vector_index1list[vectorii]=atoi(args[2].c_str());
      vector_type2list[vectorii]=system->show_atomtype_index(args[3]);
      vector_index2list[vectorii]=atoi(args[4].c_str());

      if(vector_specieslist[vectorii]==-1)
      {
	cout << "Species " << args[0] << " does not exist. Please check vector list for errors.\n";
	exit(0);
      }

      if(vector_type1list[vectorii]==-1)
      {
	cout << "Atom type " << args[1] << " does not exist. Please check vector list for errors.\n";
	exit(0);
      }

      if(vector_type2list[vectorii]==-1)
      {
	cout << "Atom type " << args[3] << " does not exist. Please check vector list for errors.\n";
	exit(0);
      }

      if(vector_index1list[vectorii]>=system->show_molecule(vector_specieslist[vectorii], 0)->typecount(vector_type1list[vectorii]))
      {
	cout << "Atom index "<< vector_index1list[vectorii] << " greater than number of atoms of type " << args[1] << " in species " << args[0] <<". Please check vector list for errors.\n";
	exit(0);
      }

      if(vector_index2list[vectorii]>=system->show_molecule(vector_specieslist[vectorii], 0)->typecount(vector_type2list[vectorii]))
      {
	cout << "Atom index "<< vector_index2list[vectorii] << " greater than number of atoms of type " << args[3] << " in species " << args[0] <<". Please check vector list for errors.\n";
	exit(0);
      }
    }
  }
}


void Vector_Autocorrelation::calculate_mean_vector_length()
{
  mean_length = 0;
  int vector_count=0;
  for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    for(int vectorii=0;vectorii<n_vectors;vectorii++)
    {
      for(int moleculeii=0;moleculeii<system->show_n_molecules(vector_specieslist[vectorii]);moleculeii++)
      {
	mean_length += (system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type2list[vectorii], vector_index2list[vectorii])->show_coordinate(timeii)-system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type1list[vectorii], vector_index1list[vectorii])->show_coordinate(timeii)).length();
	vector_count++;
      }
    }
  }
  mean_length/=float(vector_count);
}



void Vector_Autocorrelation::list_displacementkernel(int timegapii,int thisii,int nextii)
{
  int vectorii,moleculeii;
  Coordinate thiscoord1, thiscoord2, nextcoord1, nextcoord2;
  for(vectorii=0;vectorii<n_vectors;vectorii++)
  {
    for(moleculeii=0;moleculeii<system->show_n_molecules(vector_specieslist[vectorii]);moleculeii++)
    {
      thiscoord1 = system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type1list[vectorii], vector_index1list[vectorii])->show_coordinate(thisii);
      thiscoord2 = system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type2list[vectorii], vector_index2list[vectorii])->show_coordinate(thisii);
      nextcoord1 = system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type1list[vectorii], vector_index1list[vectorii])->show_coordinate(nextii);
      nextcoord2 = system->show_molecule(vector_specieslist[vectorii], moleculeii)->show_atom_trajectory(vector_type2list[vectorii], vector_index2list[vectorii])->show_coordinate(nextii);

      correlation[timegapii]+=(nextcoord2-nextcoord1)&(thiscoord2-thiscoord1);
      orientational_correlation[timegapii]+=(nextcoord2-nextcoord1).unit_vector()&(thiscoord2-thiscoord1).unit_vector();
      weighting[timegapii]++;
    }
  }
}


/*Method to write correlation data to file*/

void Vector_Autocorrelation::write(string filename)
{
  int timeii, vectorii;

  cout << "\nWriting vector time-autocorrelation function to file "<<filename<<".\n";

  ofstream output(filename.c_str());

  output << "Vector time-autocorrelation data created by AMDAT v." << VERSION << "\n";

  output << "\ntime\ttotal_correlation\torientational_correlation\n";

  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<correlation[timeii]<<"\t"<<orientational_correlation[timeii]<<"\n";
  }

  output << "Calculated for the following vectors:\n";

  output<<"species\tatomtype1\tatomindex1\tatomtype2\tatomindex2\n";

  for(vectorii=0;vectorii<n_vectors;vectorii++)
  {
    output << vector_specieslist[vectorii]<<"\t"<<vector_type1list[vectorii]<<"\t"<<vector_index1list[vectorii]<<"\t"<<vector_type2list[vectorii]<<"\t"<<vector_index2list[vectorii]<<"\n";
  }

}


