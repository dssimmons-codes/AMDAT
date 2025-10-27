// Copyright (C) 2010-2025 David S. Simmons

/*Written by David S. Simmons*/

#include "intermediate_scattering_function.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "generated/version.h"

using namespace std;


/*variable "parallel" sets whether multiple equally-spaced pairs between the same pair of exponential blocks are looped over, or if timegaps between blocks are taken only between the initial element of each block*/

Intermediate_Scattering_Function::Intermediate_Scattering_Function()
{
  //primitive data members
  n_times=1;
  n_spacebins=0;
  first_bin_index=0;
  last_bin_index=0;
  bin_size=0.0;
  firsttime=0;
  lasttime=0;
  n_atoms_represented=0;

  //pointers
  system=0;
  wavedensity1=0;
  wavedensity2=0;
  wavevectors=0;
  trajectory_list=0;

  //arrays
  timegap_weighting = new int [n_times];
  n_atoms = new float [n_times];
  correlation = new float * [n_times];
  timetable = new float [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    timegap_weighting[timeii]=0;
    n_atoms[timeii]=0;
    correlation[timeii]=0;
    timetable[timeii]=0;
    correlation[timeii] = new float[n_spacebins];
  }

}

/*Method to calculated the intermediate scattering function for all wavenumbers and all timegaps*/
Intermediate_Scattering_Function::Intermediate_Scattering_Function(System* sys, Wave_Density * wd1, Wave_Density * wd2, bool parallel)
{
  int timeii, binii;
  int wavenumberii;
  system = sys;
  wavedensity1 = wd1;
  wavedensity2 = wd2;
  if(wavedensity1->show_wavevectors()==wavedensity2->show_wavevectors())
  {
    wavevectors = wavedensity1->show_wavevectors();
    bin_size = wavedensity1->show_wavevectors()->show_delta_wavenumber();
  }
  else
  {

    cout  << "\nError. Wave vectors for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_system()==wavedensity2->show_system())
  {
    system = wavedensity1->show_system();
  }
  else
  {
    cout << "\nError. Systems for density distributions do not match!\n";
    exit(1);
  }

  if(wavedensity1->show_first()==wavedensity2->show_first())
  {
    first_bin_index = wavedensity1->show_first();
  }
  else
  {
    cout << "\nError. Wave vector bins for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_last()==wavedensity2->show_last())
  {
    last_bin_index = wavedensity1->show_last();
  }
  else
  {
    cout << "\nError. Wave vector bins for density distributions do not match!";
    exit(1);
  }

    n_atoms_represented=wavedensity1->show_n_atoms_looped();



  n_spacebins = last_bin_index-first_bin_index+1;
  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();


  timegap_weighting = system->timegap_weighting(parallel);



  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = new float[n_spacebins];
    n_atoms[timeii]=0;
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }

  /*Calculate intermediate scattering function*/
  Trajectory * traj=0;
  system->displacement_loop(this,traj,parallel);


  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    for(binii=0;binii<n_spacebins;binii++)
    {
            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));//normalize by number of atoms
    }
  }
}



/*----------------------------------------------------------------------------------*/


/*Method to calculate the intermediate scattering function for a specific wavenumber, wn, for all time gaps*/
Intermediate_Scattering_Function::Intermediate_Scattering_Function(System* sys, Wave_Density * wd1, Wave_Density * wd2, int wn, bool parallel)
{
  int timeii, wavenumberii, binii;
  system = sys;
  wavedensity1 = wd1;
  wavedensity2 = wd2;

  if(wavedensity1->show_wavevectors()==wavedensity2->show_wavevectors())
  {
    wavevectors = wavedensity1->show_wavevectors();
    bin_size = wavedensity1->show_wavevectors()->show_delta_wavenumber();
  }
  else
  {
    cout  << "Error. Wave vectors for density distributions do not match!\n";
    exit(1);
  }

  if(wavedensity1->show_system()==wavedensity2->show_system())
  {
    system = wavedensity1->show_system();
  }
  else
  {
    cout << "Error. Systems for density distributions do not match!\n";
    exit(1);
  }

  first_bin_index = wn;
  last_bin_index = wn;

  if(first_bin_index<wavedensity1->show_first()||first_bin_index<wavedensity2->show_first())
  {
    cout << "Error: requested wavenumber below range of that contained in fourier density object.";
    exit(1);
  }

  if(first_bin_index>wavedensity1->show_last()||first_bin_index>wavedensity2->show_last())
  {
    cout << "\nError: requested wavenumber above range of that contained in fourier density object.";
    exit(1);
  }

  n_atoms_represented=wavedensity1->show_n_atoms_looped();


  n_spacebins = last_bin_index-first_bin_index+1;
  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();

  timegap_weighting = system->timegap_weighting(parallel);

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = new float[n_spacebins];
    n_atoms[timeii]=0;
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }

  /*Calculate intermediate scattering function*/
  Trajectory * traj=0;
  system->displacement_loop(this,traj,parallel);


  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    for(binii=0;binii<n_spacebins;binii++)
    {
            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));		//normalize by number of atoms
    }
  }
}



/*----------------------------------------------------------------------------------*/



/*Method to calculate the intermediate scattering function at all timegaps, for a range of wavenumbers bounded by wn1 and wn2*/
Intermediate_Scattering_Function::Intermediate_Scattering_Function(System* sys, Wave_Density * wd1, Wave_Density * wd2, int wn1, int wn2, bool parallel)
{
  int timeii, wavenumberii, binii;
  system = sys;
  wavedensity1 = wd1;
  wavedensity2 = wd2;

  if(wavedensity1->show_wavevectors()==wavedensity2->show_wavevectors())
  {
    wavevectors = wavedensity1->show_wavevectors();
    bin_size = wavedensity1->show_wavevectors()->show_delta_wavenumber();
  }
  else
  {
    cout  << "\nError. Wave vectors for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_system()==wavedensity2->show_system())
  {
    system = wavedensity1->show_system();
  }
  else
  {
    cout << "Error. Systems for density distributions do not match!\n";
    exit(1);
  }

  first_bin_index = wn1;
  last_bin_index = wn2;

  if(first_bin_index<wavedensity1->show_first()||first_bin_index<wavedensity2->show_first())
  {
    cout << "Error: requested wavenumber below range of that contained in fourier density object.";
    exit(1);
  }

  n_atoms_represented=wavedensity1->show_n_atoms_looped();

  if(first_bin_index>wavedensity1->show_last()||first_bin_index>wavedensity2->show_last())
  {
    cout << "Error: requested wavenumber above range of that contained in fourier density object.";
    exit(1);
  }

  n_spacebins = last_bin_index-first_bin_index+1;
  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;

  timetable = system->displacement_times();
  timegap_weighting = system->timegap_weighting(parallel);

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = new float[n_spacebins];
    n_atoms[timeii]=0;
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
  /*Calculate intermediate scattering function*/
  Trajectory * traj=0;
  system->displacement_loop(this,traj,parallel);


  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    for(binii=0;binii<n_spacebins;binii++)
    {

            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));	//normalize by number of atoms

    }
  }
}



/*----------------------------------------------------------------------------------*/


/*Method to calculated the intermediate scattering function for all wavenumbers and at a given timegap.  timegapii=0 gives the structure factor.*/
Intermediate_Scattering_Function::Intermediate_Scattering_Function(System* sys, int timegapii, Wave_Density * wd1, Wave_Density * wd2, bool parallel)
{
  int timeii, wavenumberii, binii;
  system = sys;
  wavedensity1 = wd1;
  wavedensity2 = wd2;

  if(wavedensity1->show_wavevectors()==wavedensity2->show_wavevectors())
  {
    wavevectors = wavedensity1->show_wavevectors();
    bin_size = wavedensity1->show_wavevectors()->show_delta_wavenumber();
  }
  else
  {
    cout  << "\nError. Wave vectors for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_system()==wavedensity2->show_system())
  {
    system = wavedensity1->show_system();
  }
  else
  {
    cout << "\nError. Systems for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_first()==wavedensity2->show_first())
  {
    first_bin_index = wavedensity1->show_first();
  }
  else
  {
    cout << "\nError. Wave vector bins for density distributions do not match!";
    exit(1);
  }

  if(wavedensity1->show_last()==wavedensity2->show_last())
  {
    last_bin_index = wavedensity1->show_last();
  }
  else
  {
    cout << "\nError. Wave vector bins for density distributions do not match!";
    exit(1);
  }

  n_atoms_represented=wavedensity1->show_n_atoms_looped();

  n_spacebins = last_bin_index-first_bin_index+1;
  n_times = system->show_n_timegaps();
  firsttime = timegapii;
  lasttime = timegapii;
  timetable = system->displacement_times();


  timegap_weighting = system->timegap_weighting(parallel);

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = new float[n_spacebins];
    n_atoms[timeii]=0;
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }

  /*Calculate intermediate scattering function*/
  Trajectory * traj=0;
  system->displacement_loop(this,traj,timegapii,parallel);

  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    for(binii=0;binii<n_spacebins;binii++)
    {
            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));//normalize by number of atoms
    }
  }
}

/*----------------------------------------------------------------------------------*/

#ifdef NEVER

/*Method to calculate the intermediate scattering function at all timegaps, for a range of wavenumbers bounded by wn1 and wn2, based on a time-dependent particle list given by particle_list.  It will base the wavevectors and system off of the wavedensity that is provided.*/
Intermediate_Scattering_Function::Intermediate_Scattering_Function(System* sys, Particle_List * particle_list, Wave_Density * wd2, int wn1, int wn2)
{

  int timeii, wavenumberii, binii;
  system = sys;
  wavedensity2 = wd2;

  wavevectors = wavedensity2->show_wavevectors();
  wavedensity1 = new Wave_Density;
  wavedensity1->set(system,wavevectors,wn1,wn2);

  wavedensity1->atomlist(particle_list);
  first_bin_index = wn1;
  last_bin_index = wn2;

  if(first_bin_index<wavedensity1->show_first()||first_bin_index<wavedensity2->show_first())
  {
    cout << "Error: requested wavenumber below range of that contained in fourier density object.";
    exit(1);
  }

  n_atoms_represented = 0;

  if(first_bin_index>wavedensity1->show_last()||first_bin_index>wavedensity2->show_last())
  {
    cout << "Error: requested wavenumber above range of that contained in fourier density object.";
    exit(1);
  }

  n_spacebins = last_bin_index-first_bin_index+1;
  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;

  timetable = system->displacement_times();
  timegap_weighting = system->timegap_weighting(false);

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = new float[n_spacebins];
    n_atoms[timeii]=0;
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }

//  cout <<"\n"<< n_spacebins;
  /*Calculate intermediate scattering function*/

  Trajectory * traj=0;
  system->displacement_loop(this,traj,bool(0));



  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
//    cout << "\n"<<n_atoms[timeii]<<"\t"<<timegap_weighting[timeii];
    for(binii=0;binii<n_spacebins;binii++)
    {

            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));//normalize by number of atoms

    }
  }
}

#endif

Intermediate_Scattering_Function::~Intermediate_Scattering_Function()
{
  delete [] n_atoms;
  for(int timeii=firsttime;timeii<=lasttime;timeii++)
  {
      delete [] correlation[timeii];
  }
  delete [] correlation;
}

Intermediate_Scattering_Function::Intermediate_Scattering_Function(const Intermediate_Scattering_Function & copy)
{
 system = copy.system;
 wavedensity1 = copy.wavedensity1;
 wavedensity2 = copy.wavedensity2;
 timegap_weighting = copy.timegap_weighting;
 timetable=copy.timetable;
 wavevectors=copy.wavevectors;
 for(int timeii=firsttime;timeii<=lasttime;timeii++)
 {
       n_atoms[timeii] = copy.n_atoms[timeii];
    for(int binii=0;binii<n_spacebins;binii++)
    {

            correlation[timeii][binii]/=float(n_atoms[timeii]*float(timegap_weighting[timeii]));//normalize by number of atoms

    }
 }


}

Intermediate_Scattering_Function Intermediate_Scattering_Function::operator = (const Intermediate_Scattering_Function & copy)
{
  if(this!=&copy)
  {
    delete [] n_atoms;
    for(int timeii=firsttime;timeii<=lasttime;timeii++)
    {
	delete [] correlation[timeii];
    }
    delete [] correlation;

    system = copy.system;
    wavedensity1 = copy.wavedensity1;
    wavedensity2 = copy.wavedensity2;
    timegap_weighting = copy.timegap_weighting;
    timetable=copy.timetable;
    wavevectors=copy.wavevectors;

    for(int timeii=firsttime;timeii<=lasttime;timeii++)
    {
	  n_atoms[timeii] = copy.n_atoms[timeii];
	for(int binii=0;binii<n_spacebins;binii++)
	{
	  correlation[timeii][binii] = copy.correlation[timeii][binii];		//normalize by number of atoms
	}
    }
  }
  return *this;
}

/*----------------------------------------------------------------------------------*/

/*displacementkernel method to be called by displacementloop methods in System class*/
void Intermediate_Scattering_Function::displacementkernel(int timegapii, int thisii, int nextii, Trajectory * traj)
{
  int wavenumberii;
  int wavevectorii;
  int vectorcount;

  for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
  {
    vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index);
    for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
    {
      correlation[timegapii][wavenumberii] += (wavedensity1->show_density(thisii,wavenumberii,wavevectorii).real()*wavedensity2->show_density(nextii,wavenumberii,wavevectorii).real() + wavedensity1->show_density(thisii,wavenumberii,wavevectorii).imag()*wavedensity2->show_density(nextii,wavenumberii,wavevectorii).imag())/(float(vectorcount));
    }
  }
  n_atoms[timegapii]+=float(wavedensity1->show_n_atoms(thisii))/float(timegap_weighting[timegapii]);
}
