/*Methods for Particle_List class - stores identifying data for a list of atoms*/
/*Amorphous Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#include "particle_list.h"
#include <fstream>

using namespace std;

Particle_List::Particle_List(int n_times, int capacity)
{
  int timeii;

  n_times = 0;
  capacity = 0;
  species_index = new int * [n_times];
  molecule_index = new int * [n_times];
  type_index = new int * [n_times];
  atom_index = new int * [n_times];
  n_atoms = new int [n_times];
  trajectoryIDs = new int * [n_times];

  for(timeii=0;timeii<n_times;timeii++)
  {
    species_index[timeii] = new int [capacity];
    molecule_index[timeii] = new int [capacity];
    type_index[timeii] = new int [capacity];
    atom_index[timeii] = new int [capacity];
    trajectoryIDs[timeii] = new int [capacity];
    n_atoms[timeii] = 0;
  }

}


void Particle_List::clear()
{
  int timeii;

  for(timeii=0;timeii<n_times;timeii++)
  {
    delete species_index[timeii];
    delete molecule_index[timeii];
    delete type_index[timeii];
    delete atom_index[timeii];
  }

  delete [] species_index;
  delete [] molecule_index;
  delete [] type_index;
  delete [] atom_index;
  delete [] n_atoms;
}



void Particle_List::atomloop(Analysis * analysis)
{
  int timeii;
  int atomii;

  for(timeii=0;timeii<n_times;timeii++)
  {
    for(atomii=0;atomii<n_atoms[timeii];atomii++)
    {
      analysis->atomlist_kernel(timeii,species_index[timeii][atomii],molecule_index[timeii][atomii],type_index[timeii][atomii],atom_index[timeii][atomii],this);
    }
  }
}



void Particle_List::atomloop(Analysis * analysis, int timeii)
{
  int atomii;

  for(atomii=0;atomii<n_atoms[timeii];atomii++)
  {
    analysis->atomlist_kernel(timeii,species_index[timeii][atomii],molecule_index[timeii][atomii],type_index[timeii][atomii],atom_index[timeii][atomii],this);
  }

}


void Particle_List::write_count(string filename)
{
  double avg_atoms=0;
  int timeii;

  /*calculate average number of atoms*/
  for(timeii=0;timeii<n_times;timeii++)
  {
    avg_atoms += double(n_atoms[timeii])/double(n_times);
  }

  ofstream output(filename.c_str());
  output << "Average_atoms:\t"<< avg_atoms<<"\n";
  output << "Count_List\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << "\t" << n_atoms[timeii] << "\n";
  }
}


int Particle_List::show_atom(int timeii, int atomii)
{
	int tempID;
	if(timeii>n_times)
	{
		tempID=-2;
	}
	else
	{
		if(atomii<n_atoms[timeii])
		{
			tempID=trajectoryIDs[timeii][atomii];
		}
		else
		{
			tempID=-1;
		}
	}
	return tempID;
}
