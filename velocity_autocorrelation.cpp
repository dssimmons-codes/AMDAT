/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate velocity autocorrelation function*/
/*Written by David S. Simmons*/

#include "velocity_autocorrelation.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#ifndef TACC
#include <fftw3.h>
#endif

using namespace std;




/*Constructor*/
Velocity_Autocorrelation::Velocity_Autocorrelation(Mean_Square_Displacement* mean_square_displacement)
{
  initialize(mean_square_displacement);
}


/*Constructor*/
Velocity_Autocorrelation::Velocity_Autocorrelation(Mean_Square_Displacement* mean_square_displacement, int maxtimestep)
{
  initialize(mean_square_displacement, maxtimestep);
}



/*Method to initialize object*/
void Velocity_Autocorrelation::initialize(Mean_Square_Displacement* mean_square_displacement)
{
  int tii;

  msd = mean_square_displacement;
  system = msd->show_system();
  n_times = system->show_n_timegaps();
  timetable = system->displacement_times();
  vac = new float [n_times];
  linearity_error = 0;
  time_unit = system->show_time_unit();
  fourier_space_vac_real = new double [n_times];
  fourier_space_vac_imag = new double [n_times];
  fouriersize = 0;
  // Initialize velocity autocorrelation data to zero
  for(tii=0;tii<n_times;tii++)
  {
    vac[tii]=0;
    fourier_space_vac_real[tii] = 0;
    fourier_space_vac_imag[tii] = 0;
  }
  
  calculate();
  
}


/*Method to initialize object*/
void Velocity_Autocorrelation::initialize(Mean_Square_Displacement* mean_square_displacement, int ntimes)
{
  int tii;

  msd = mean_square_displacement;
  system = msd->show_system();
  
  if(ntimes>=system->show_n_timegaps())
  {
    n_times = system->show_n_timegaps();
  }
  else
  {
    n_times = ntimes;
  }
 
  time_unit = system->show_time_unit();
  timetable = system->displacement_times();
  vac = new float [n_times];
  linearity_error = 0;
  fourier_space_vac_real = new double [n_times];
  fourier_space_vac_imag = new double [n_times];
  fouriersize=0;
  // Initialize velocity autocorrelation data to zero
  for(tii=0;tii<n_times;tii++)
  {
    vac[tii]=0;
    fourier_space_vac_real[tii]=0;
    fourier_space_vac_imag[tii]=0;
  }
  
  calculate();
}


/*Method to calculate velocity autocorrelation function via central finite difference*/
void Velocity_Autocorrelation::calculate()
{
  int tii;

  /*derivative for zeroth time step*/
  vac[0]=(2*msd->show(0)-5*msd->show(1)+4*msd->show(2)-msd->show(3))/(time_unit*time_unit);

  /*loop over main time steps, taking second derivative of msd using central difference*/
  for(tii=1;tii<n_times-1;tii++)
  {
    vac[tii]=(msd->show(tii-1)-2*msd->show(tii)+msd->show(tii+1))/(time_unit*time_unit);
  }
  
  /*derivative for last time step*/
  vac[n_times-1]=(2*msd->show(n_times-1)-5*msd->show(n_times-2)+4*msd->show(n_times-3)-msd->show(n_times-4))/(time_unit*time_unit);
  
}
#ifndef TACC
void Velocity_Autocorrelation::fourier_transform()
{

  fftw_plan plan;
  int transformsize = n_times;
  double * fft_in;
  fftw_complex *fft_out;
  int tii;
  
  
  fouriersize=floor(transformsize/2)+1;
  fft_in = new double [n_times];
  for(tii=0;tii<n_times;tii++)
  {
    fft_in[tii]=double(vac[tii]);
  }
  fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fouriersize);
  plan = fftw_plan_dft_r2c_1d(transformsize, fft_in, fft_out,FFTW_ESTIMATE);
  fftw_execute(plan);
  for(tii=0;tii<fouriersize;tii++)
  {
    fourier_space_vac_real[tii] = fft_out[tii][0];
    fourier_space_vac_imag[tii] = fft_out[tii][1];
  }
  
}
#endif
void Velocity_Autocorrelation::write(string filename)const
{
  int timeii;
  
  cout << "\nWriting velocity autocorrelation function to file.";
  ofstream output(filename.c_str());
  
  output << "Velocity autocorrelation function data created by MDAT v." << VERSION << "\n"; 
  output << "\nThis algorithm assumes linearly spaced timesteps with spacing given by timestep in input file.\n Exponential or block time spacing schemes will yield incorrect results; be careful.\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
      output << timetable[timeii]<<"\t"<<vac[timeii]<<"\n";
  }
  output.close();
}


void Velocity_Autocorrelation::write_fourier(string filename)const
{
int tii;
  
  cout << "\nWriting fourier transform velocity autocorrelation function to file.";
  ofstream output(filename.c_str());
  
  output << "Fourier transform of velocity autocorrelation function data created by MDAT v." << VERSION << "\n"; 
  output << "\nThis algorithm assumes linearly spaced timesteps with spacing given by timestep in input file.\n Exponential or block time spacing schemes will yield incorrect results; be careful.\n";
  for(tii=0;tii<fouriersize;tii++)
  {
      output<<(tii+1)/timetable[n_times-1] <<"\t"<<fourier_space_vac_real[tii]<<"\t"<<fourier_space_vac_imag[tii]<<"\n";
  }
  output.close();
}
