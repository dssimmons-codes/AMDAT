/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to calculate mean-square-displacement*/
/*Written by David S. Simmons*/

#include "orientational_correlation.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "multibody_list.h"
#include <omp.h>
#include "system.h"


using namespace std;


Orientational_Correlation::Orientational_Correlation()
{
  n_times=0;

  correlation = new float [n_times];
  weighting = new int [n_times];
  overall_correlation = 0;
  weighting = 0;
  correlated_vector.set(0,0,0);
  multibody_list=0;
}


Orientational_Correlation::Orientational_Correlation(const Orientational_Correlation & copy)
{
  int timeii;
  system = copy.system;
  n_times = copy.n_times;
  multibody_list = copy.multibody_list;

  correlation = new float [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = copy.correlation[timeii];
    weighting[timeii] = copy.weighting[timeii];
  }
  overall_correlation=copy.overall_correlation;
  correlated_vector=copy.correlated_vector;
}


/** **/
Orientational_Correlation::Orientational_Correlation(System*sys)
{
  int timeii;
  system = sys;
  n_times = system->show_n_timesteps();
  correlation = new float [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = 0;
    weighting[timeii] = 0;
  }
  overall_correlation=0;
  correlated_vector.set(0,0,0);
}


Orientational_Correlation::Orientational_Correlation(System*sys, Coordinate vec)
{
  int timeii;
  system = sys;
  n_times = system->show_n_timesteps();
  correlation = new float [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = 0;
    weighting[timeii] = 0;
  }
  overall_correlation=0;
  correlated_vector.set(0,0,0);
  correlated_vector=vec;
}


Orientational_Correlation Orientational_Correlation::operator = (const Orientational_Correlation & copy)
{
  int timeii;

  if(this!=&copy)
  {

    system = copy.system;
  n_times = copy.n_times;
  multibody_list = copy.multibody_list;

  correlation = new float [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii] = copy.correlation[timeii];
    weighting[timeii] = copy.weighting[timeii];
  }
  overall_correlation=copy.overall_correlation;
  correlated_vector=copy.correlated_vector;

  }

  return *this;

}



void Orientational_Correlation::analyze(Multibody_List * t_list)
{
  multibody_list=t_list;
  for(int timeii=0; timeii<system->show_n_timesteps(); timeii++)
  {
    weighting[timeii]+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0, timeii, 0);
  }
  postprocess_list();
}


void Orientational_Correlation::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  float dotproduct=  (((*current_multibody)(1)->show_unwrapped(thisii)-(*current_multibody)(0)->show_unwrapped(thisii)).unit_vector())&(correlated_vector).unit_vector());	//compute dot product between unit vectors at initial and later times
  correlation[thisii]+=0.5*(3.0*dotproduct*dotproduct - 1.0);	//increment baf by second legendre polynomial of dot product above
}


void Orientational_Correlation::postprocess_list()
{
  float cumulative_weighting=0;
  for(int timeii=0;timeii<n_times;timeii++)
  {
        overall_correlation+=correlation[timeii];
        correlation[timeii] /= float(weighting[timeii]);
        cumulative_weighting+=weighting[timeii];
  }
  overall_correlation/=cumulative_weighting;


}




void Orientational_Correlation::write(string filename)const
{
  int timeii;

  cout << "\nWriting orientational_correlation to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Orientation correlation function data created by AMDAT v." << VERSION << "\n";
  output << "Mean correlation with vector " << correlated_vector.show_x() << " " <<correlated_vector.show_y() << correlated_vector.show_z()<<"\n";

    output << "overall\t"<<overall_correlatoin[timeii]<<"\n";

    for(timeii=0;timeii<n_times;timeii++)
  {
    output << system->show_time(timeii)<<"\t"<<correlation[timeii]<<"\n";
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
