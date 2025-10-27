// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_Bead_Region class methods - stores list of multibody objects that are within a given region*/
/*Written by David S. Simmons*/


#include "multibody_bead_region.h"
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "system.h"

using namespace std;

Multibody_Bead_Region::Multibody_Bead_Region()
{
  int timeii;

  n_bodies=-2;

  n_times = 0;

  time_conversion = new int [1];

  for(timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=timeii;
  }
}


 Multibody_Bead_Region::Multibody_Bead_Region(System* syst)
{
  int timeii;

  system = sys = syst;

  n_bodies=-2;

  n_times = system->show_n_timesteps();

  time_conversion = new int [system->show_n_timesteps()];

  for(timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=timeii;
  }
  multibodies.resize(system->show_n_timesteps());
}


Multibody_Bead_Region::Multibody_Bead_Region(System* syst, Coordinate low, Coordinate high, int thresh)
{
  int timeii;

  system = sys = syst;
  lowerbound = low;
  upperbound = high;
  threshold = thresh;
  n_bodies=-2;

  n_times = system->show_n_timesteps();

  time_conversion = new int [system->show_n_timesteps()];

  for(timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=timeii;
  }
  multibodies.resize(system->show_n_timesteps());
}

Multibody_Bead_Region::Multibody_Bead_Region(const Multibody_Bead_Region & copy)
{

    system=sys=copy.sys;

    n_times=copy.n_times;
    n_bodies=copy.n_bodies;
    threshold=copy.threshold;

     time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
    multibodies=copy.multibodies;
   upperbound=copy.upperbound;
  lowerbound=copy.lowerbound;


}


Multibody_Bead_Region Multibody_Bead_Region::operator=(const Multibody_Bead_Region & copy)
{
  if(this!=&copy)
  {

    delete [] time_conversion;
    system=sys=copy.sys;

    n_times=copy.n_times;
    n_bodies=copy.n_bodies;
    threshold=copy.threshold;

     time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
    multibodies=copy.multibodies;
   upperbound=copy.upperbound;
  lowerbound=copy.lowerbound;
  }

  return *this;
}



void Multibody_Bead_Region::analyze(Multibody_List * m_list)
{
  int timeii;

  multibody_list=m_list;

  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    multibody_list->listloop(this,0, timeii, 0);
  }

  postprocess_list();

}


void Multibody_Bead_Region::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  int n_bodies;
  int bodies_in=0;

  n_bodies = current_multibody->show_n_bodies();

  for(int bodyii=0;bodyii<n_bodies;bodyii++)
  {
      bodies_in+=((*current_multibody)(bodyii)->show_coordinate(thisii)).within(lowerbound,upperbound);
  }

  if(bodies_in>=threshold)
  {
    multibodies[convert_time(thisii)].push_back(current_multibody);
  }
}

void Multibody_Bead_Region::postprocess_list()
{
  check_n_bodies();
}


void Multibody_Bead_Region::write(string filename)const
{
  	float avg_multibodies=0;
	int timeii;
	int realtimeii;

	/*calculate average number of atoms*/
	for(timeii=0;timeii<n_times;timeii++)
	{
		avg_multibodies += (multibodies[timeii].size())/n_times;
	}

	ofstream output(filename.c_str());
	output << "Multibody list statistics created by AMDAT v." << amdat::build::SEMVER << "\n";
	output << "Average_multibodies in list:\t"<< avg_multibodies<<"\n";
	output << "Count_List\n";
	for(timeii=0;timeii<n_times;timeii++)
	{
		output << system->show_time(timeii) << "\t" << multibodies[timeii].size() << "\n";
	}
}
