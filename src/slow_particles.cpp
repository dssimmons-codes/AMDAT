/*Methods for Slow_Particles class: identifies and stores list of slow particles at start of each exponential block*/
/*Amorphous Molecular dynamics analysis toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#include "slow_particles.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;

Slow_Particles::Slow_Particles()
{
	int timeii;
  
	displacement_time_index=-1;
	maxdistance = -1;
	n_times = 0;
  
	capacity = 0;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
	}
  
	atoms_considered = 0;
}



Slow_Particles::Slow_Particles(System * sys, Gaussian_Comparison * gc)
{
	int timeii;

	system = sys;
	gaussian_comparison = gc;
	displacement_time_index = gaussian_comparison->show_time_index();
	maxdistance = gaussian_comparison->show_slowboundary();
	capacity=system->show_n_atoms()+system->show_n_molecules();
	n_times = system->show_n_exponentials();
  
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
	}
	atoms_considered = 0;
}



void Slow_Particles::set(System * sys, Gaussian_Comparison * gc)
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
		delete [] trajectories[timeii];
	}
	delete [] trajectories;
	delete [] n_trajectories;
  
	system = sys;
	gaussian_comparison = gc;
	displacement_time_index = gaussian_comparison->show_time_index();
	maxdistance = gaussian_comparison->show_slowboundary();
	capacity=system->show_n_atoms()+system->show_n_molecules();
	n_times = system->show_n_exponentials();
  
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
	}
  
	atoms_considered = 0;
}


void Slow_Particles::atomkernel(Trajectory * traj)
{
	system->displacement_loop(this, traj ,displacement_time_index,bool(0));
	atoms_considered++;
}



void Slow_Particles::displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj)
{
	float square_displacement;		//particle displacement
	int expindex;			//index of exponential block
  
	square_displacement = pow(traj->distance(thisii,nextii),2);  //calculate particle displacement
  
	expindex = int((float(thisii)-float(system->show_frt()))/float(system->show_n_exponential_steps()));	//calculate which exponential block this corresponds to
	if(square_displacement < maxdistance)
	{
		if(n_trajectories[expindex]==capacity){cout<<"Error: particle list memory allocation full.\n";exit(1);}
		trajectories[expindex][n_trajectories[expindex]]=traj;
		n_trajectories[expindex]++;
	}
}

