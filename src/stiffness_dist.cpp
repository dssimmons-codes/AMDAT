
#include "stiffness_dist.h"
#include <fstream>
#include "system.h"
#include <math.h>
#include <iostream>
#include "version.h"

using namespace std;


Stiffness_Dist::Stiffness_Dist(System* sys, int bins, float maxvalue, float t)
{
	int binii;
	
	system = sys;
  
	n_bins = bins+1;
	binsize = maxvalue / bins;
	time_index = t;
  
	time = system->displacement_times()[time_index];
	timeweighting = system->timegap_weighting()[time_index];
  
	distribution = new float [n_bins];
  
	atomcount = 0;
	
	for(binii=0;binii<n_bins;binii++)
	{
		distribution[binii]=0;
	}
  
	mean = 0;
	square_term = 0;
	log_mean=0;
	log_square_term=0;
}


/*---------Methods to perform analysis with trajectory list-------------*/


void Stiffness_Dist::analyze(Trajectory_List* t_list)
{
	trajectory_list = t_list;
	system->displacement_list(this, time_index);
	postprocess_list();
}


void Stiffness_Dist::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	atomcount += trajectory_list[0].show_n_trajectories(thisii)/float(timeweighting);
	(trajectory_list[0]).listloop(this,thisii);
}


void Stiffness_Dist::listkernel(Trajectory * current_trajectory)
{
	float stiffness;		//variable to temporarily hold mean-square displacement of atom
	int index;
  
	stiffness = pow(current_trajectory->distance(currenttime,nexttime),-2.0);  //determine msd of this atom

	mean+=stiffness/float(timeweighting);
	square_term += stiffness*stiffness/float(timeweighting);
	
	log_mean+=log(stiffness)/float(timeweighting);
	log_square_term += log(stiffness)*log(stiffness)/float(timeweighting);
	
	index = int(stiffness/binsize);
	if(index>=n_bins){index=n_bins-1;}		//use top bin as slough bin
	distribution[index]+=1.0/float(timeweighting);
	
}

void Stiffness_Dist::postprocess_list()
{
	int binii;
  
	mean /= atomcount;
	variance = square_term/atomcount-mean*mean;
	
	log_mean/=atomcount;
	log_variance = log_square_term/atomcount-log_mean*log_mean;
  
	for(binii=0;binii<n_bins;binii++)
	{
		distribution[binii]/=atomcount;
	}
}



/*-------Method to write stiffness distribution and statistics to data file-------*/

void Stiffness_Dist::write(string filename)const
{
	int binii;
  
	ofstream output (filename.c_str());		//open correlation file
	cout << "\nWriting distribution of molecular stiffnesses to file " <<filename<<"." ;
  
	output << "Stiffness distribution data created by AMDAT v." << VERSION << "\n"; 
	output << "MSD time = " << time << "\n";
	output << "Mean_stiffness\t" << mean <<"\n";
	output << "Variance\t" << variance<<"\n\n";
	output << "Log-mean_stiffness\t"<<log_mean<<"\n";
	output << "Log-variance\t"<<log_variance<<"\n";
  
	for(binii=0;binii<n_bins;binii++)
	{
		output << (float(binii)+0.5)*binsize << "\t" << distribution[binii] << "\n"; 
	}
  
	output << "Last bin includes overflow.";
  
	output.close();
  
}

void Stiffness_Dist::write(ofstream& output)const
{
	int binii;
  
	cout << "\nWriting distribution of molecular stiffnesses to file." ;
  
	output << "Stiffness distribution data created by AMDAT v." << VERSION << "\n"; 
	output << "MSD time = " << time << "\n";
	output << "Mean_stiffness\t" << mean <<"\n";
	output << "Variance\t" << variance<<"\n\n";
	output << "Log-mean_stiffness\t"<<log_mean<<"\n";
	output << "Log-variance\t"<<log_variance<<"\n";
  
	for(binii=0;binii<n_bins;binii++)
	{
		output << (float(binii)+0.5)*binsize << "\t" << distribution[binii] << "\n"; 
	}
  
	output << "Last bin includes overflow.";
  
	output.close();
  
}
