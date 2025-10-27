// Copyright (C) 2010-2025 David S. Simmons

/*Density_Fourier Class methods.  This class calculates and stores the spacial fourier transform of the number density for each wavevector passed from a Wave_Vectors object.*/
/*Written by David S. Simmons*/

#include <math.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "wave_density.h"
#include "generated/version.h"

using namespace std;

Wave_Density::Wave_Density()
{
  system=0;
  trajectory_list=0;
  n_wavenumbers = 0;
  n_times = 0;
  first_wavenumber_index = 0;
  last_wavenumber_index = 0;
  wavevectors=0;

  density = new complex<double> ** [n_times];
  n_atoms = new int[n_times];
  for(int tii=0;tii<n_times;tii++)
  {
    n_atoms[tii]=0;
    density[tii] = new complex<double> * [n_wavenumbers];
    for(int wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(int vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        ((density[tii][wavenumberii-first_wavenumber_index][vectorii]).real(0));
        ((density[tii][wavenumberii-first_wavenumber_index][vectorii]).imag(0));
      }
    }
  }
}

/*Constructor to perpare to calculate wave density for all wavenumbers*/
Wave_Density::Wave_Density(System * sys, const Wave_Vectors * wv, bool stf)
{
  int tii, wavenumberii, vectorii;
  sf=stf;
  system = sys; //set system for which wave density will be calculated
  trajectory_list=0;
  wavevectors = wv;

  first_wavenumber_index = 0;
  n_wavenumbers = wavevectors->show_n_wavenumbers();
  last_wavenumber_index = n_wavenumbers-1;

  n_times = system->show_n_timesteps();		//get number of times in trajectory from system object
  /*allocate memory and zero out densities*/
  density = new complex<double> ** [n_times];
  n_atoms = new int [n_times];
  for(tii=0;tii<n_times;tii++)
  {

    n_atoms[tii]=0;
    density[tii] = new complex<double> * [n_wavenumbers];
    for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].real(0));
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].imag(0));
      }
    }
  }

  /*zero out number of atoms looped over*/
  n_atoms_looped = 0;
}

/*Constructor to prepare to calculate wave density over a set range of wavenumber indices*/
Wave_Density::Wave_Density(System * sys, const Wave_Vectors * wv, int inner, int outer)
{
  int tii, wavenumberii, vectorii;

  system = sys; //set system for which wave density will be calculated
  trajectory_list=0;
  wavevectors = wv;

  /*define range of wavevector 'bin' indices to pull from Wave_Vector object; for a spherical geometry this corresponds to the indices of the desired wavenumbers*/
  first_wavenumber_index = inner;
  last_wavenumber_index = outer;
  n_wavenumbers = last_wavenumber_index - first_wavenumber_index + 1;

  n_times = system->show_n_timesteps();		//get number of times in trajectory from system object

  /*allocate memory and zero out densities*/
  density = new complex<double> ** [n_times];
  n_atoms = new int[n_times];
  for(tii=0;tii<n_times;tii++)
  {
    n_atoms[tii]=0;
    density[tii] = new complex<double> * [n_wavenumbers];
    for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
	    density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].real(0));
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].imag(0));
      }
    }
  }

  /*zero out number of atoms looped over*/
  n_atoms_looped = 0;
}

Wave_Density::Wave_Density(const Wave_Density & copy)
:Analysis(copy)
{
//   system = copy.system;
//   trajectory_list = copy.trajectory_list;
  wavevectors = copy.wavevectors;
  first_wavenumber_index = copy.first_wavenumber_index;
  last_wavenumber_index = copy.last_wavenumber_index;
  n_wavenumbers = copy.n_wavenumbers;
  n_times = copy.n_times;
  n_atoms_looped = copy.n_atoms_looped;
  atomset = copy.atomset;
  sf = copy.sf;
  currenttime = copy.currenttime;

  density = new complex<double> **[n_times];
  n_atoms = new int[n_times];

  for(int tii=0;tii<n_times;tii++)
  {
    n_atoms[tii] = copy.n_atoms[tii];
    density[tii] = new complex<double> * [n_wavenumbers];
    for(int wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(int vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        density[tii][wavenumberii-first_wavenumber_index][vectorii]=copy.density[tii][wavenumberii-first_wavenumber_index][vectorii];
      }
    }
  }
}

Wave_Density Wave_Density::operator = (const Wave_Density & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    trajectory_list = copy.trajectory_list;

    first_wavenumber_index = copy.first_wavenumber_index;
    last_wavenumber_index = copy.last_wavenumber_index;
    n_wavenumbers = copy.n_wavenumbers;
    n_times = copy.n_times;
    n_atoms_looped = copy.n_atoms_looped;
    atomset = copy.atomset;
    sf = copy.sf;
    currenttime = copy.currenttime;

    for(int tii=0;tii<n_times;tii++)
    {
      for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
      {
	delete [] density[tii][wavenumberii];
      }
      delete [] density[tii];
    }
    delete [] density;
    delete [] n_atoms;
    
    wavevectors = copy.wavevectors;

    density = new complex<double> **[n_times];
    n_atoms = new int[n_times];

    for(int tii=0;tii<n_times;tii++)
    {
      n_atoms[tii] = copy.n_atoms[tii];
      density[tii] = new complex<double> * [n_wavenumbers];
      for(int wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
      {
	density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
	for(int vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
	{
	  density[tii][wavenumberii-first_wavenumber_index][vectorii]=copy.density[tii][wavenumberii-first_wavenumber_index][vectorii];
	}
      }
    }
  }
  return *this;
}


Wave_Density::~Wave_Density()
{
  for(int tii=0;tii<n_times;tii++)
  {
    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      delete [] density[tii][wavenumberii];
    }
    delete [] density[tii];
  }
  delete [] density;
  delete [] n_atoms;
}

void Wave_Density::clear_memory()
{
  int tii, wavenumberii;
  for(tii=0;tii<n_times;tii++)
  {
    for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      delete [] density[tii][wavenumberii-first_wavenumber_index];
    }
    delete [] density[tii];
  }
  delete [] density;
  delete [] n_atoms;
}


void Wave_Density::set(System * sys, const Wave_Vectors * wv)
{

  int tii, wavenumberii, vectorii;

  clear_memory();

  system = sys; //set system for which wave density will be calculated
  wavevectors = wv;

  first_wavenumber_index = 0;
  n_wavenumbers = wavevectors->show_n_wavenumbers();
  last_wavenumber_index = n_wavenumbers-1;

  n_times = system->show_n_timesteps();		//get number of times in trajectory from system object

  /*allocate memory and zero out densities*/
  density = new complex<double> ** [n_times];
  n_atoms = new int[n_times];
  for(tii=0;tii<n_times;tii++)
  {
    n_atoms[tii] = 0;
    density[tii] = new complex<double> * [n_wavenumbers];
    for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].real(0));
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].imag(0));
      }
    }
  }

  /*zero out number of atoms looped over*/
  n_atoms_looped = 0;
}



void Wave_Density::set(System * sys, const Wave_Vectors * wv, int inner, int outer)
{
  int tii, wavenumberii, vectorii;

  clear_memory();

  system = sys;		//set system for which wave density will be calculated
  wavevectors = wv;

  /*define range of wavevector 'bin' indices to pull from Wave_Vector object; for a spherical geometry this corresponds to the indices of the desired wavenumbers*/
  first_wavenumber_index = inner;
  last_wavenumber_index = outer;
  n_wavenumbers = last_wavenumber_index - first_wavenumber_index + 1;

  n_times = system->show_n_timesteps();		//get number of times in trajectory from system object

  /*allocate memory and zero out densities*/
  density = new complex<double> ** [n_times];
  n_atoms = new int[n_times];
  for(tii=0;tii<n_times;tii++)
  {
    n_atoms[tii]=0;
    density[tii] = new complex<double> * [n_wavenumbers];
    for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)
    {
      density[tii][wavenumberii-first_wavenumber_index]=new complex<double> [wavevectors->vectorcount(wavenumberii)];
      for(vectorii=0;vectorii<wavevectors->vectorcount(wavenumberii);vectorii++)
      {
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].real(0));
        (density[tii][wavenumberii-first_wavenumber_index][vectorii].imag(0));
      }
    }
  }
  /*zero out number of atoms looped over*/
  n_atoms_looped = 0;
}




void Wave_Density::analyze(Trajectory_List * t_list)
{
  int tii;
  if(sf)
  {
    trajectory_list=t_list;
    for(tii=0;tii<system->show_n_exponentials();tii++)
    {
      currenttime=system->show_n_exponential_steps()*tii;
      n_atoms[currenttime]+=trajectory_list->show_n_trajectories(currenttime);
      trajectory_list->listloop(this,currenttime);
    }
  }
  else
  {
    trajectory_list=t_list;
    for(tii=0;tii<n_times;tii++)
    {
      currenttime=tii;
      n_atoms[currenttime]+=trajectory_list->show_n_trajectories(currenttime);
      trajectory_list->listloop(this,currenttime);
    }
  }
}


void Wave_Density::listkernel(Trajectory* current_trajectory)
{
  int wavenumberii, vectorii;
  vector<Coordinate>vectorlist;
  Coordinate coordinate;
  int vectorcount;
  float k_dot_r;
  complex<double> tempcomplex;

  coordinate = current_trajectory->show_coordinate(currenttime);	//get atom coordinate at give time

  for(wavenumberii=first_wavenumber_index;wavenumberii<=last_wavenumber_index;wavenumberii++)		//loop over wavenumbers for which wavedensity will be calculated
  {
    vectorlist = wavevectors->vectorlist(wavenumberii);	//call up wave vector list for this wavenumber
    vectorcount = wavevectors->vectorcount(wavenumberii);	//get count of wave vectors for this wavenumber
    for(vectorii=0;vectorii<vectorcount;vectorii++)		//loop over wavevectors for this wavenumber
    {
      k_dot_r = vectorlist[vectorii]&coordinate;		//calculate dot product of wave vector and present atomecoordinate
      tempcomplex.real(cos(k_dot_r));
      tempcomplex.imag(sin(k_dot_r));
      (density[currenttime][wavenumberii-first_wavenumber_index][vectorii])+=tempcomplex;
      //density[currenttime][wavenumberii-first_wavenumber_index][vectorii].real() += cos(k_dot_r);	//add contribution to real part of wave density
      //density[currenttime][wavenumberii-first_wavenumber_index][vectorii].imag() += sin(k_dot_r);	//add contribution to imaginary parto f wave density
//        if(tii==1&&wavenumberii==25&&molecule_index==0&&atom_index==0){cout<<k_dot_r<<"\n";}
    }
  }
}


