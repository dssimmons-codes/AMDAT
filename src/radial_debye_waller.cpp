/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Radial Debye Waller class*/
/*Written by David S. Simmons*/

#include "radial_debye_waller.h"
#include "atom_trajectory.h"
#include "molecule.h"
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "version.h"

#define PI 3.14159265

using namespace std;

Radial_Debye_Waller::Radial_Debye_Waller(System* sys, int timeii, int bincount, float maxrange, Coordinate c)
{
  int binii;

  system = sys;
  time_index = timeii;
  n_bins = bincount;
  bin_size = maxrange / float(n_bins);

  n_atoms = new float [n_bins];
  msd = new float [n_bins];
  density = new float [n_bins];
  for(binii=0;binii<n_bins;binii++)
  {
    msd[binii]=0;
    n_atoms[binii]=0;
  }

  time_weighting = system->timegap_weighting(bool(0))[timeii];

  center = c;
}




/*-------Methods to do analysis with trajectory list--------*/

void Radial_Debye_Waller::analyze(Trajectory_List * t_list)
{
	trajectory_list=t_list;
	system->displacement_list(this);
	postprocess_list();
}

void Radial_Debye_Waller::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	(trajectory_list[0]).listloop(this,currenttime);
}

void Radial_Debye_Waller::listkernel(Trajectory* current_trajectory)
{
	float radius;
	float distance;
	int bin;
	Coordinate coordinate1, coordinate2;

	coordinate1 = current_trajectory->show_coordinate(currenttime);
	coordinate2 = current_trajectory->show_coordinate(nexttime);

	radius = (coordinate1-center).length();
	distance = (coordinate2-coordinate1).length();
	bin = int(radius/bin_size);
	if(bin>=n_bins){bin=n_bins-1;}
	n_atoms[bin]++;
	msd[bin] += distance*distance;
}

void Radial_Debye_Waller::postprocess_list()
{
	int binii;
	float shellvolume;
	for(binii=0;binii<n_bins;binii++)
	{
		msd[binii]/=n_atoms[binii];
		n_atoms[binii]/=time_weighting;
		shellvolume = 4/3*PI*pow((float(binii)+1.0)*bin_size,3)-4/3*PI*pow((float(binii))*bin_size,3);
		density[binii] = n_atoms[binii]/shellvolume;
	}
}


/*-----------Method to write radial Debye-Waller factor to file------*/


void Radial_Debye_Waller::write(string filename)
{
int binii;

  cout << "\nWriting radial debye waller data to file " << filename << ".";

  ofstream output(filename.c_str());

  output << "Radial Debye-Waller data created by AMDAT v." << VERSION << "\n";
  output <<"time="<< system->show_time(time_index) <<"\n";
  output <<"radius\tn_atoms\tdensity\tmsd\n";

  for(binii=0;binii<n_bins;binii++)
  {
    output << (float(binii)+.5)*bin_size<<"\t"<<n_atoms[binii]<<"\t"<<density[binii]<<"\t"<<msd[binii]<<"\n";
  }
}


void Radial_Debye_Waller::write(ofstream& output)const
{
int binii;

  cout << "\nWriting radial debye waller data to file.";

  output << "Radial Debye-Waller data created by AMDAT v." << VERSION << "\n";
  output <<"time="<< system->show_time(time_index) <<"\n";
  output <<"radius\tn_atoms\tdensity\tmsd\n";

  for(binii=0;binii<n_bins;binii++)
  {
    output << (float(binii)+.5)*bin_size<<"\t"<<n_atoms[binii]<<"\t"<<density[binii]<<"\t"<<msd[binii]<<"\n";
  }
}
