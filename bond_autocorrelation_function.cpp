/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "bond_autocorrelation_function.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "multibody_list.h"
#include <omp.h>
#include "system.h"


using namespace std;


Bond_Autocorrelation_Function::Bond_Autocorrelation_Function()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  atomcount = 0;
}


Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(const Bond_Autocorrelation_Function & copy)
{
  int timeii;

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
}

Bond_Autocorrelation_Function::~Bond_Autocorrelation_Function()
{
  delete [] baf;
  delete [] weighting;
  delete [] timetable;
}


/** **/
Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(System*sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;

}




Bond_Autocorrelation_Function Bond_Autocorrelation_Function::operator = (const Bond_Autocorrelation_Function & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Bond_Autocorrelation_Function::initialize(System* sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
}



/*Methods to do analysis using trajectory list*/

void Bond_Autocorrelation_Function::analyze(Multibody_List * t_list)
{
  multibody_list=t_list;
  system->displacement_list(this);
  postprocess_list();
}

void Bond_Autocorrelation_Function::list_displacementkernel(int timegapii,int thisii, int nextii)
{
  weighting[timegapii]+=multibody_list->show_n_multibodies(thisii);
  multibody_list->listloop(this,timegapii, thisii, nextii);
}



void Bond_Autocorrelation_Function::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  float dotproduct=  (((*current_multibody)(1)->show_unwrapped(thisii)-(*current_multibody)(0)->show_unwrapped(thisii)).unit_vector())&(((*current_multibody)(1)->show_unwrapped(nextii)-(*current_multibody)(0)->show_unwrapped(nextii)).unit_vector());	//compute dot product between unit vectors at initial and later times
  baf[timegapii]+=0.5*(3.0*dotproduct*dotproduct - 1.0);	//increment baf by second legendre polynomial of dot product above
}


void Bond_Autocorrelation_Function::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        baf[timeii] /= float(weighting[timeii]);

  }
}



/*Method to write MSD data to file*/

void Bond_Autocorrelation_Function::write(string filename)const
{
  int timeii;

  cout << "\nWriting baf to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Bond autocorrelation function data created bys AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}

void Bond_Autocorrelation_Function::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting baf to file.";

  output << "Bond autocorrelation function data created bys AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}

#ifdef NEVER
void Bond_Autocorrelation_Function::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  multibody_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Bond_Autocorrelation_Function::postprocess_bins()
{
  postprocess_list();
}

#endif