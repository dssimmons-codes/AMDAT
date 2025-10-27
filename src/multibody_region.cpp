// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_Region class methods - stores list of multibody objects that are within a given region*/
/*Written by David S. Simmons*/


#include "multibody_region.h"
#include <iostream>
#include <stdlib.h>
#include "generated/version.h"
#include "system.h"

using namespace std;

Multibody_Region::Multibody_Region()
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


 Multibody_Region::Multibody_Region(System* syst)
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


Multibody_Region::Multibody_Region(System* syst, Coordinate low, Coordinate high)
{
  int timeii;
  
  system = sys = syst;
  lowerbound = low;
  upperbound = high;
  
  n_bodies=-2;
  
  n_times = system->show_n_timesteps();
  
  time_conversion = new int [system->show_n_timesteps()];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=timeii;
  }
  multibodies.resize(system->show_n_timesteps());
}

Multibody_Region::Multibody_Region(const Multibody_Region & copy)
{
  
    system=sys=copy.sys;

    n_times=copy.n_times;
    n_bodies=copy.n_bodies;

     time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
    multibodies=copy.multibodies;
   upperbound=copy.upperbound;
  lowerbound=copy.lowerbound;
 
 
}


Multibody_Region Multibody_Region::operator=(const Multibody_Region & copy)
{
  if(this!=&copy)
  {

    delete [] time_conversion;
    system=sys=copy.sys;

    n_times=copy.n_times;
    n_bodies=copy.n_bodies;

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



void Multibody_Region::analyze(Multibody_List * m_list)
{
  int timeii;
  
  multibody_list=m_list;
  
  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    multibody_list->listloop(this,0, timeii, 0);
  }
  
  postprocess_list();
  
}


void Multibody_Region::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  if((current_multibody->show_coordinate(thisii)).within(lowerbound,upperbound))
  {
    multibodies[convert_time(thisii)].push_back(current_multibody);
  }
}

void Multibody_Region::postprocess_list()
{
  check_n_bodies();
}


void Multibody_Region::write(string filename)const
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