/*Methods for Particle_List_Exptime class - stores identifying data for a list of atoms at the start time of eaxh exponential time block*/
/*Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#include "particle_list_exptime.h"
#include <fstream>

using namespace std;

Particle_List_Exptime::Particle_List_Exptime(int n_times, int capacity)
{
  int timeii;
  
  n_times = 0;
  capacity = 0;
  species_index = new int * [n_times];
  molecule_index = new int * [n_times];
  type_index = new int * [n_times];
  atom_index = new int * [n_times];
  n_atoms = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    species_index[timeii] = new int [capacity];
    molecule_index[timeii] = new int [capacity];
    type_index[timeii] = new int [capacity];
    atom_index[timeii] = new int [capacity];
    n_atoms[timeii] = 0;
  }
}



void Particle_List_Exptime::atomloop(Analysis * analysis)
{
  int expii, timeii;
  int atomii;
  for(expii=0;expii<n_times;expii++)
  {
    timeii = (expii*system->show_n_exponential_steps()) + system->show_frt();
    for(atomii=0;atomii<n_atoms[expii];atomii++)
    {
      cout << species_index[expii][atomii] << "\t" <<molecule_index[expii][atomii]<<"\t"<<type_index[expii][atomii]<<"\t"<<atom_index[expii][atomii]<<"\n";cout.flush();
      analysis->atomlist_kernel(timeii,species_index[expii][atomii],molecule_index[expii][atomii],type_index[expii][atomii],atom_index[expii][atomii],this);
    }
  }
}




void Particle_List_Exptime::atomloop(Analysis * analysis, int timeii)
{
  int expii;
  int atomii;

  expii = int(float(timeii - system->show_frt())/float(system->show_n_exponential_steps()));
  for(atomii=0;atomii<n_atoms[expii];atomii++)
  {
    analysis->atomlist_kernel(timeii,species_index[expii][atomii],molecule_index[expii][atomii],type_index[expii][atomii],atom_index[expii][atomii],this);
  }
}




int Particle_List_Exptime::show_n_atoms(int timeindex)
{
  int expii = int(float(timeindex - system->show_frt())/float(system->show_n_exponential_steps()));
  return n_atoms[expii];
}




void Particle_List_Exptime::write_count(string filename)
{
  double avg_atoms=0;
  int expii, timeii;
  /*calculate average number of atoms*/
  for(expii=0;expii<n_times;expii++)
  {
    avg_atoms += double(n_atoms[expii])/double(n_times);
  }

  ofstream output(filename.c_str());
  output << "Average_atoms:\t"<< avg_atoms<<"\n";
  output << "Count_List\n";
  for(expii=0;expii<n_times;expii++)
  {
    timeii = (expii*system->show_n_exponential_steps()) + system->show_frt();
    output <<system->show_time(timeii)<<"\t"<< n_atoms[expii] << "\n";
  }
}



int Particle_List_Exptime::show_n_atoms(int timeindex, int atomii)
{
	int timeii = int(float(timeindex - system->show_frt())/float(system->show_n_exponential_steps()));
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