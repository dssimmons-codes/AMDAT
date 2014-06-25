/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate distribution of particle displacements to user specified power*/
/*Written by David S. Simmons*/

#include "displacement_distribution.h"
#include <fstream>
#include "system.h"
#include <math.h>
#include <iostream>
#include "version.h"

using namespace std;


/*constructor*/
Displacement_Distribution::Displacement_Distribution(System* sys, float pow, int bins, float maxvalue, float t)
{
  int binii;

  system = sys;
  power=pow;

  n_bins = bins+1;
  binsize = maxvalue / bins;
  time_index = t;

  time = system->displacement_times()[time_index];
  timeweighting = system->timegap_weighting()[time_index];

  distribution = new float [n_bins];

  for(binii=0;binii<n_bins;binii++)
  {
    distribution[binii]=0;
  }

  mean = 0;
  square_term = 0;

  atomcount = 0;

}


/*----------------------------------------------------------------*/
/*Methods to do calculation on list of trajectories*/

void Displacement_Distribution::analyze(Trajectory_List* t_list)
{
  trajectory_list = t_list;
  system->displacement_list(this, time_index);
  postprocess_list();
}


void Displacement_Distribution::list_displacementkernel(int timegapii, int thisii, int nextii)
{
  currenttime=thisii;
  nexttime=nextii;
  atomcount += trajectory_list[0].show_n_trajectories(thisii)/float(timeweighting)                                                                          ;
  (trajectory_list[0]).listloop(this,thisii);
}


void Displacement_Distribution::listkernel(Trajectory * current_trajectory)
{
  float msd;		//variable to temporarily hold mean-square displacement of atom
  int index;

  msd = pow(current_trajectory->distance(currenttime,nexttime),2.0);  //determine msd of this atom

  mean += msd/float(timeweighting);
  square_term += msd*msd/float(timeweighting);

  index = int(msd/binsize);
  if(index>=n_bins){index=n_bins-1;}		//use top bin as slough bin

  distribution[index]+=1.0/float(timeweighting);
}

void Displacement_Distribution::postprocess_list()
{
  int binii;

  mean /= atomcount;
  variance = square_term/atomcount-mean*mean;

  for(binii=0;binii<n_bins;binii++)
  {
    distribution[binii]/=atomcount;
  }
}


/*----------------------------------------------------------------*/
/*Method to write distribution and statistics to file*/

void Displacement_Distribution::write(string filename)const
{
  int binii;

  ofstream output (filename.c_str());		//open correlation file
  cout << "\nWriting distribution of particle displacement magnitudes to the "<<power<<" power to file " <<filename<<"." ;

  output << "Particle displacement distribution to the "<< power <<  " power created by AMDAT v." << VERSION << "\n";
  output << "Last bin includes overflow.";
  output << "time = " << time << "\n";
  output << "Mean\t" << mean <<"\n";
  output << "Variance\t" << variance<<"\n\n";

  for(binii=0;binii<n_bins;binii++)
  {
    output << (float(binii)+0.5)*binsize << "\t" << distribution[binii] << "\n";
  }

  output.close();

}
