// Copyright (C) 2010-2025 David S. Simmons

/*Methods for Exptime_Trajectory_List class - stores a list of trajectories that is the same for all times*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include "exptime_trajectory_list.h"

using namespace std;

Exptime_Trajectory_List::Exptime_Trajectory_List()
{

	sys=0;
	system=const_cast<System*>(sys);

	capacity = 0;

	n_times=0;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	included = new Boolean_List [n_times];
	time_conversion = new int [n_times];



}

Exptime_Trajectory_List::Exptime_Trajectory_List(System* syst, int capacity)
{
  int timeii;

  sys=syst;
  system=const_cast<System*>(sys);

  if(capacity==0){capacity=system->show_n_atoms()+system->show_n_molecules();}

  n_times=system->show_n_exponentials();
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(system);
  }

  time_conversion=new int [system->show_n_timesteps()];
  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=int(float(timeii - sys->show_frt())/float(system->show_n_exponential_steps()));
  }
}



void Exptime_Trajectory_List::reset(System* syst, int cap)
{
  int timeii;

  sys=syst;
  system=const_cast<System*>(sys);

  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] n_trajectories;
  delete [] time_conversion;

  if(cap==0){capacity=system->show_n_atoms()+system->show_n_molecules();}
  else{capacity=cap;}

  n_times=system->show_n_exponentials();
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(system);
  }

  time_conversion=new int [system->show_n_timesteps()];
  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=int(float(timeii - sys->show_frt())/float(system->show_n_exponential_steps()));
  }
}

void Exptime_Trajectory_List::write_count(string filename)const
{
	double avg_trajectories=0;
	int timeii;
	int realtimeii;

	/*calculate average number of atoms*/
	for(timeii=0;timeii<n_times;timeii++)
	{
		avg_trajectories += double(n_trajectories[timeii])/double(n_times);
	}

	ofstream output(filename.c_str());
	output << "Average_trajectories:\t"<< avg_trajectories<<"\n";
	output << "Count_List\n";
	for(timeii=0;timeii<n_times;timeii++)
	{
		realtimeii = (timeii*system->show_n_exponential_steps()) + system->show_frt();
		output << system->show_time(realtimeii) << "\t" << n_trajectories[timeii] << "\n";
	}
}
