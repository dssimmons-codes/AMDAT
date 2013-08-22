/*Methods fo Structure_Factor class - calculates structure factor*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons*/


#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include "version.h"
#include "structure_factor.h"
#include <omp.h>


using namespace std;


Structure_Factor::Structure_Factor()
{
  int vectorcount;
  trajlist1=0;
  trajlist2=0;
  n_wavenumbers=0;
  current_wavedensity=0;
  wavevectors=0;
  n_wavenumbers = 0;
  n_atoms=0;
  currenttime=0;

  structure_factor = new float [n_wavenumbers];
  wavedensity1 = new complex<float> * [n_wavenumbers];
  wavedensity2 = new complex<float> * [n_wavenumbers];
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    vectorcount = wavevectors->vectorcount(wavenumberii);
    wavedensity1[wavenumberii]= new complex<float> [vectorcount];
    wavedensity2[wavenumberii]= new complex<float> [vectorcount];
    for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
    {
      wavedensity1[wavenumberii][wavevectorii].real() = 0;
      wavedensity1[wavenumberii][wavevectorii].imag() = 0;
      wavedensity2[wavenumberii][wavevectorii].real() = 0;
      wavedensity2[wavenumberii][wavevectorii].imag() = 0;
    }
  }
}


/*---------------------------------------------------------------------------------*/
// Standard Constructor
Structure_Factor::Structure_Factor(System* sys, const Wave_Vectors * wv, bool fblock)
{
  system = sys;
  int vectorcount;
  trajlist1=0;
  trajlist2=0;
  current_wavedensity=0;
  wavevectors = wv;
  fullblock=fblock;
  n_wavenumbers = wavevectors->show_n_wavenumbers();
  structure_factor = new float [n_wavenumbers];
  wavedensity1 = new complex<float> * [n_wavenumbers];
  wavedensity2 = new complex<float> * [n_wavenumbers];

  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    vectorcount = wavevectors->vectorcount(wavenumberii);
    wavedensity1[wavenumberii]= new complex<float> [vectorcount];
    wavedensity2[wavenumberii]= new complex<float> [vectorcount];
    for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
    {
      wavedensity1[wavenumberii][wavevectorii].real() = 0;
      wavedensity1[wavenumberii][wavevectorii].imag() = 0;
      wavedensity2[wavenumberii][wavevectorii].real() = 0;
      wavedensity2[wavenumberii][wavevectorii].imag() = 0;
    }
  }
}

Structure_Factor::~Structure_Factor()
{
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
   {
     delete [] wavedensity1[wavenumberii];
     delete [] wavedensity2[wavenumberii];
  }
  delete [] structure_factor;
  delete [] wavedensity1;
  delete [] wavedensity2;
}
/*---------------------------------------------------------------------------------------------*/



Structure_Factor::Structure_Factor(const Structure_Factor & copy)
:Analysis(copy)
{
  int vectorcount;
    fullblock = copy.fullblock;
    system = copy.system;
    trajlist1=copy.trajlist1;
    trajlist2=copy.trajlist2;
    n_wavenumbers=copy.n_wavenumbers;
    wavevectors=copy.wavevectors;
    n_atoms=copy.n_atoms;

    structure_factor = new float [n_wavenumbers];
    wavedensity1 = new complex<float> * [n_wavenumbers];
    wavedensity2 = new complex<float> * [n_wavenumbers];

    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      wavedensity1[wavenumberii]=new complex<float> [vectorcount];
      wavedensity2[wavenumberii]=new complex<float> [vectorcount];
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
	wavedensity1[wavenumberii][wavevectorii] = copy.wavedensity1[wavenumberii][wavevectorii];
	wavedensity2[wavenumberii][wavevectorii] = copy.wavedensity2[wavenumberii][wavevectorii];
      }
    }
}


Structure_Factor Structure_Factor::operator=(const Structure_Factor & copy)
{
  if(this!=&copy)
  {
    int vectorcount;
    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      delete [] wavedensity1[wavenumberii];
      delete [] wavedensity2[wavenumberii];
    }
    delete [] structure_factor;
    delete [] wavedensity1;
    delete [] wavedensity2;

    fullblock = copy.fullblock;
    system = copy.system;
    trajlist1=copy.trajlist1;
    trajlist2=copy.trajlist2;
    n_wavenumbers=copy.n_wavenumbers;
    wavevectors=copy.wavevectors;
    n_atoms=copy.n_atoms;

    structure_factor = new float [n_wavenumbers];
    wavedensity1 = new complex<float> * [n_wavenumbers];
    wavedensity2 = new complex<float> * [n_wavenumbers];

    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      wavedensity1[wavenumberii]=new complex<float> [vectorcount];
      wavedensity2[wavenumberii]=new complex<float> [vectorcount];
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
	wavedensity1[wavenumberii][wavevectorii] = copy.wavedensity1[wavenumberii][wavevectorii];
	wavedensity2[wavenumberii][wavevectorii] = copy.wavedensity2[wavenumberii][wavevectorii];
      }
    }
  }
  return *this;
}

/*calculates asymmetric structure factor*/
void Structure_Factor::analyze(Trajectory_List * t_list1, Trajectory_List * t_list2)
{
  int timeii, wavenumberii,  wavevectorii;
  int timeincrement, n_times;
  int vectorcount;

  trajlist1=t_list1;
  trajlist2=t_list2;

  //maxtime = system->show_n_timesteps();

  n_times = system->show_n_timesteps();

  n_atoms=0;

  if(fullblock)
  {
    timeincrement = 1;
  }
  else
  {
    timeincrement = system->show_n_exponential_steps();
  }


  /*allocate memory for wavedensity and structure factor and zero out structure factor*/
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    structure_factor[wavenumberii] = 0;
    wavedensity1[wavenumberii]=new complex<float> [wavevectors->vectorcount(wavenumberii)];
    wavedensity2[wavenumberii]=new complex<float> [wavevectors->vectorcount(wavenumberii)];
  }


  for(timeii=0;timeii<n_times;timeii+=timeincrement)
  {
    //Zero out wave densities
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
      {
        wavedensity1[wavenumberii][wavevectorii].real() = 0;
        wavedensity1[wavenumberii][wavevectorii].imag() = 0;
        wavedensity2[wavenumberii][wavevectorii].real() = 0;
        wavedensity2[wavenumberii][wavevectorii].imag() = 0;
      }
    }

    /*calculate wave densities*/
    currenttime=timeii;
    current_wavedensity = wavedensity1;
    analyze_wave_density(trajlist1);
    current_wavedensity = wavedensity2;
    analyze_wave_density(trajlist2);

    //Calculate structure factor
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
        structure_factor[wavenumberii] += (wavedensity1[wavenumberii][wavevectorii].real()*wavedensity2[wavenumberii][wavevectorii].real() + wavedensity1[wavenumberii][wavevectorii].imag()*wavedensity2[wavenumberii][wavevectorii].imag())/(float(vectorcount));
      }
    }

    n_atoms+=trajlist1->show_n_trajectories(timeii);
  }

  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    structure_factor[wavenumberii]/=float(n_atoms);
  }
}

/*---------------------------------------------------------------------------------*/

/*calculates symmetric structure factor*/
void Structure_Factor::analyze(Trajectory_List * t_list1)
{
  int timeii, wavenumberii,  wavevectorii;
  int timeincrement, n_times;
  int vectorcount;

  trajlist1=t_list1;
  trajlist2=0;

  n_times = system->show_n_timesteps();

  n_atoms=0;

  if(fullblock)
  {
    timeincrement = 1;
  }
  else
  {
    timeincrement = system->show_n_exponential_steps();
  }


  /*allocate memory for wavedensity and structure factor and zero out structure factor*/
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    structure_factor[wavenumberii] = 0;
    wavedensity1[wavenumberii]=new complex<float> [wavevectors->vectorcount(wavenumberii)];
  }


//cout <<endl << timeincrement << endl;cout.flush();


  for(timeii=0;timeii<n_times;timeii+=timeincrement) //TODO: This is probably the ideal spot, as it will call analyze_wave_density. Alternatively you could make listloop parallel and see how that works out
  {
    //Zero out wave densities
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
      {
        wavedensity1[wavenumberii][wavevectorii].real() = 0;
        wavedensity1[wavenumberii][wavevectorii].imag() = 0;
      }
    }

    /*calculate wave densities*/
    currenttime=timeii;
    current_wavedensity = wavedensity1;
    analyze_wave_density(trajlist1); 	//TODO: Examine how much time this call takes, the listloop may need to be parallel(and may be the majority of the time spent)
										//Results: takes 3-4 seconds for each iteration of the loop. Seems to be 99% of the total time. Parallelization in listloop would be good
    //Calculate structure factor
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors->vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
        structure_factor[wavenumberii] += (wavedensity1[wavenumberii][wavevectorii].real()*wavedensity1[wavenumberii][wavevectorii].real() + wavedensity1[wavenumberii][wavevectorii].imag()*wavedensity1[wavenumberii][wavevectorii].imag())/(float(vectorcount));
      }
    }
    n_atoms+=trajlist1->show_n_trajectories(timeii);
  }

  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
      if (float(n_atoms)!=0)
      {
         structure_factor[wavenumberii]/=float(n_atoms);
      }
    else
    {
        structure_factor[wavenumberii]=0;
    }
  }
}

void Structure_Factor::analyze_wave_density(Trajectory_List * t_list)
{
  t_list->listloop(this,currenttime);
}

void Structure_Factor::listkernel(Trajectory* current_trajectory)
{
  int wavenumberii, vectorii;
  Coordinate const * vectorlist;
  Coordinate coordinate;
  int vectorcount;
  float k_dot_r;

  coordinate = current_trajectory->show_coordinate(currenttime);	//get atom coordinate at give time // the currenttime reference is... actually probably ok.
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)		//loop over wavenumbers for which wavedensity will be calculated
  {
    vectorlist = wavevectors->vectorlist(wavenumberii);	//call up wave vector list for this wavenumber //The wavevectors reference may be what's breaking parallelization
    vectorcount = wavevectors->vectorcount(wavenumberii);	//get count of wave vectors for this wavenumber
    for(vectorii=0;vectorii<vectorcount;vectorii++)		//loop over wavevectors for this wavenumber
    {
      k_dot_r = vectorlist[vectorii]&coordinate;		//calculate dot product of wave vector and present atomecoordinate
	  #pragma omp flush(current_wavedensity)
      current_wavedensity[wavenumberii][vectorii].real() += cos(k_dot_r);	//add contribution to real part of wave density // This COULD break paralleliation. Might need to omp flush the variable?
//	  #pragma omp flush(current_wavedensity[wavenumerii][vectorii])
//	  #pragma omp flush(current_wavedensity)
      current_wavedensity[wavenumberii][vectorii].imag() += sin(k_dot_r);	//add contribution to imaginary part of wave density
	  #pragma omp flush(current_wavedensity)
    }
  }
}


/*Write correlation object to file*/
void Structure_Factor::write(string filename)const
{
  int binii;

  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Structure factor data created by AMDAT v." << VERSION << "\n";
  for(binii=0;binii<n_wavenumbers;binii++)
  {
    output << wavevectors->show_approx_wavenumber(binii) << "\t" << structure_factor[binii]  <<  endl;
  }
  output << "\n";

  output.close();
}
