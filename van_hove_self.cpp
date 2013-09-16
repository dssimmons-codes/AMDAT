/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for Van_Hove_Self: a class for self part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#include "van_hove_self.h"
#include <stdlib.h>
#include <iostream>
#include "version.h"

using namespace std;


Van_Hove_Self::Van_Hove_Self(System* sys, int bin_count, float value_min, float value_max)
{
	initialize(sys, bin_count, value_min, value_max);
}



Van_Hove_Self::Van_Hove_Self()
{
	system = 0;
	n_bins = 0;
	min_value=0;
	max_value=0;
	bin_size=0;
	n_times=0;
	timetable=0;

	correlation = new int * [1];
	weighting = new int [0];
	correlation[0]=new int [0];
	n_atoms = new double [0];
}



void Van_Hove_Self::set(System* sys, int bin_count, float value_max, float value_min)
{
	clear_memory();
	initialize(sys, bin_count, value_max, value_min);
}


void Van_Hove_Self::initialize(System* sys, int bin_count, float value_max, float value_min)
{
	int timeii;
	int binii;

	n_bins = bin_count;
	min_value = value_min;
	system = sys;

	if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
	else max_value=value_max;

	bin_size = (max_value-min_value)/float(n_bins);

	n_times = system->show_n_timegaps();

	n_atoms_looped=0;

	timetable = system->displacement_times();

	 //allocate memory for van hove self-correlation function and weighting and initialize to zero
	correlation = new int * [n_times];
	n_atoms = new double [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		correlation[timeii]=new int [n_bins];
		n_atoms[timeii] = 0;
		for(binii=0;binii<n_bins;binii++)
		{
			correlation[timeii][binii]=0;
		}
	}
}



void Van_Hove_Self::preprocess()
{
	weighting = system->timegap_weighting();			//get timegap weighting from system
}



void Van_Hove_Self::preprocess_list(Particle_List* particle_list)
{
	weighting = system->timegap_weighting(bool(0));			//get timegap weighting from system
}


void Van_Hove_Self::displacementkernel(int timegap,int thisii, int nextii, Trajectory * traj)
{
	float distance = (traj->show_unwrapped(nextii)-traj->show_unwrapped(thisii)).length();
	bin(timegap,distance);
}



void Van_Hove_Self::atomkernel(Trajectory * traj)
{
	system->displacement_loop(this, traj);
	n_atoms_looped++;
}




void Van_Hove_Self::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;

	weighting = system->timegap_weighting(bool(0));			//get timegap weighting from system
	system->displacement_list(this);
}

void Van_Hove_Self::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	currenttime=thisii;
	nexttime = nextii;
	timegap = timegapii;
	n_atoms[timegap]+=double(trajectory_list[0].show_n_trajectories(currenttime))/double(weighting[timegap]);
	(trajectory_list[0]).listloop(this,currenttime);
}

void Van_Hove_Self::listkernel(Trajectory* current_trajectory)
{
	float distance = (current_trajectory->show_unwrapped(nexttime)-current_trajectory->show_unwrapped(currenttime)).length();
	bin(timegap,distance);


}



/*
void Van_Hove_Self::list_displacementkernel(int timegapii, int thisii, int nextii, Particle_List* particle_list)
{
	nexttime = nextii;
	timegap = timegapii;

	particle_list->atomloop(this,thisii);
}
*/


void Van_Hove_Self::postprocess()
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
		n_atoms[timeii] = double(n_atoms_looped);
	}

}

