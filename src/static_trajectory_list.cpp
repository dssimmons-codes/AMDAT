/*Methods for Static_Trajectory_List class - stores a list of trajectories that is the same for all times*/
/*Amorphous Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Written by David S. Simmons*/

#include "static_trajectory_list.h"
#include "version.h"
using namespace std;


Static_Trajectory_List::Static_Trajectory_List()
{

	capacity = 0;
    n_atomtypes = 0;
	n_times=1;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	included = new Boolean_List [n_times];
	for(int timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
	}

	//set all system times to convert to 0 internal time
	time_conversion=new int [0];
}

Static_Trajectory_List::Static_Trajectory_List(System* syst, int capacity)
{

	sys=syst;
	system=const_cast<System*>(sys);

	if(capacity==0){capacity=system->show_n_atoms()+system->show_n_molecules();}
    n_atomtypes = system->show_n_atomtypes();
	n_times=1;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	included = new Boolean_List [n_times];
	for(int timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
		included[timeii].set(system);
	}

	time_conversion=new int [system->show_n_timesteps()];
	for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
	{
	  time_conversion[timeii]=0;
	}
}



void Static_Trajectory_List::reset(System* syst, int capacity)
{
	sys=syst;
	system=const_cast<System*>(sys);

	for(int timeii=0;timeii<n_times;timeii++)
	{
		delete [] trajectories[timeii];
	}
	delete [] trajectories;
	delete [] n_trajectories;
	delete [] time_conversion;
	delete [] included;
	if(capacity==0){capacity=system->show_n_trajectories();}

	n_atomtypes = system->show_n_atomtypes();
	n_times=1;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	included = new Boolean_List [n_times];
	for(int timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
		included[timeii].set(system);
	}

	time_conversion=new int [system->show_n_timesteps()];
	for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
	{
	  time_conversion[timeii]=0;
	}
}


void Static_Trajectory_List::atomkernel(Trajectory * traj)
{
  int trajid;
  trajectories[0][n_trajectories[0]] = traj;
  trajid=trajectories[0][n_trajectories[0]]->show_trajectory_ID();
  (included[0])(trajid,bool(1));
  n_trajectories[0]++;
}




void Static_Trajectory_List::set(System * syst, Trajectory_Set * trajectory_set)
{
  int trajii, timeii;
  
  sys=syst;
  system=const_cast<System*>(sys);
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
	delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] n_trajectories;
  delete [] time_conversion;
  delete [] included;
  
  n_atomtypes = system->show_n_atomtypes();
  n_times=1;
  
  if(capacity==0){capacity=system->show_n_trajectories();}
  
  
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  n_trajectories [0] = trajectory_set->show_n_trajectories();
  included = new Boolean_List [n_times];
  included[0].set(system);
  
  time_conversion=new int [system->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=0;
  }
  
  trajectories[0] = new Trajectory * [n_trajectories[0]];
  
  for(trajii=0;trajii<n_trajectories[0];trajii++)
  {
    trajectories[0][trajii]=trajectory_set->show_trajectory(trajii);
    (included[0])(trajectory_set->show_trajectory(trajii)->show_trajectory_ID(),1);
  }
  

}


