// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Van_Hove_Self: a class for self part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#include "van_hove_self.h"
#include <stdlib.h>
#include <iostream>
#include "generated/version.h"

using namespace std;


Van_Hove_Self::Van_Hove_Self(System* sys, int bin_count, float value_max)
{
	initialize(sys, bin_count, value_max);
}



Van_Hove_Self::Van_Hove_Self()
{
	system = 0;
	n_bins = 0;
	max_value=0;
	bin_size=0;
	n_times=0;
	timetable=0;

	correlation = new float * [1];
	weighting = new int [0];
	correlation[0]=new float [0];
}



void Van_Hove_Self::set(System* sys, int bin_count, float value_max)
{
	clear_memory();
	initialize(sys, bin_count, value_max);
}


void Van_Hove_Self::initialize(System* sys, int bin_count, float value_max)
{
	int timeii;
	int binii;

	n_bins = bin_count;
	system = sys;

	if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
	else max_value=value_max;

	bin_size = (max_value)/float(n_bins);

	n_times = system->show_n_timegaps();


	timetable = system->displacement_times();

	 //allocate memory for van hove self-correlation function and weighting and initialize to zero
	correlation = new float * [n_times];
	weighting = new int [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		correlation[timeii]=new float [n_bins];
		weighting[timeii] = 0;
		for(binii=0;binii<n_bins;binii++)
		{
			correlation[timeii][binii]=0;
		}
	}
}






void Van_Hove_Self::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;

	system->displacement_list(this);
	postprocess_list();
}

void Van_Hove_Self::list_displacementkernel(int timegapii, int thisii, int nextii)
{

	//weighting[timegapii]+=trajectory_list[0].show_n_trajectories(thisii);
	(trajectory_list[0]).listloop(this, timegapii, thisii, nextii);
}

void Van_Hove_Self::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
	float distance = (current_trajectory->show_unwrapped(nextii)-current_trajectory->show_unwrapped(thisii)).length();
	bin(timegapii,distance);
	weighting[timegapii]++;


}


