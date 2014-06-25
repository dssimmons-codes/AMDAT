/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for Space-Time_Correlation_Function class: a general class to hold a time correlation function*/
/*Written by David S. Simmons*/


#include "time_correlation_function.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#ifndef TACC
#include <fftw3.h>
#endif
#include "version.h"

using namespace std;

#include "progress.h"

#ifndef PI
#define PI 3.141592653589793238462643383280
#endif




Space-Time_Correlation_Function::~Space-Time_Correlation_Function()
{
  //clear_memory();
}



void Space-Time_Correlation_Function::clear_memory()
{
 int ii;
 for(ii=0;ii<n_times;ii++)
 {
   delete [] correlation[ii];
 }
 delete [] correlation;
 delete [] weighting;
 delete [] timetable;
}



Space-Time_Correlation_Function Space-Time_Correlation_Function::operator+ (const Space-Time_Correlation_Function & increment)const
{
  int timeii, binii;
  Space-Time_Correlation_Function temp;

  if(n_bins==increment.n_bins)
  {temp.n_bins=n_bins;}
  else
  {
  cout << "\nNumber of bins do not match!\n";
  exit(1);
  }
  
  if(n_times==increment.n_times)
  {temp.n_times=n_times;}
  else
  {
  cout << "\nNumber of times do not match!\n";
  exit(1);
  }
  
  if(bin_size==increment.bin_size)
  {temp.bin_size=bin_size;}
  else
  {
  cout << "\nSize of bins do not match!\n";
  exit(1);
  }
  
  if(max_value==increment.max_value)
  {temp.max_value=max_value;}
  else
  {
  cout << "\nBin maxima do not match!\n";
  }

  if(min_value==increment.min_value)
  {temp.min_value=min_value;}
  else
  {
  cout << "\nBin minima do not match!\n";
  }

  temp.system = system;		//naively copy system point from present object; the user must simply be smart about what this means and when it is appropriate.  In general, it is acceptable if systems have the same timestep properties and merely have different trajectories.  Could put a warning in here triggered if time properties are different.
  temp.correlation = new int*[n_times];
  temp.weighting = new int[n_times];
  temp.timetable = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    temp.correlation[timeii]=new int[n_bins];
    temp.weighting[timeii]=weighting[timeii]+increment.weighting[timeii];
    for(binii=0;binii<n_bins;binii++)
    {
      temp.correlation[timeii][binii] = (correlation[timeii][binii]*float(weighting[timeii]) + increment.correlation[timeii][binii]*float(increment.weighting[timeii]))/float(temp.weighting[timeii]);
    }
    if(timetable[timeii]==increment.timetable[timeii])
    {
      temp.timetable[timeii] = timetable[timeii];
    }
    else
    {
      cout << "\nError in correlation addition: correlation data time schemes not consistent!\n";
      exit(1);
    }
  }
  
  return temp;
}



/*------------------------------------------------------------------------------*/



// Space-Time_Correlation_Function Space-Time_Correlation_Function::operator= (Space-Time_Correlation_Function increment)
// {
//   int timeii, binii;
//   Space-Time_Correlation_Function temp;
//   temp.n_bins=n_bins;
//   temp.n_times=n_times;
//   temp.bin_size=bin_size;
//   temp.max_value=max_value;
//   temp.min_value=min_value;
//   temp.system = system;
  	
  	  
//   temp.correlation = new int*[n_times];
//   temp.weighting = new int[n_times];
//   for(timeii=0;timeii<n_times;timeii++)
//   {
//     temp.correlation[timeii]=new int[n_bins];
//     for(binii=0;binii<n_bins;binii++)
//     {
//       temp.correlation[timeii][binii] = increment.correlation[timeii][binii];
//     }
//     temp.weighting[timeii] =increment.weighting[timeii];
//   }
  
//   return temp;
// }

/*------------------------------------------------------------------------------*/


#ifndef TACC
/*calculates and writes to file the spherically symmetric spacial fourier transform of the correlation data*/
float** Space-Time_Correlation_Function::spherical_fourier(string filename, int minbin)
{
  #ifdef MANUAL
  float ** normal;				//define variable to hold normalized correlation data
  normal = normalized();			//calculate normalized correlation data
  float ** fourier;				//define pointer to array of fourier transformed data
  int n_wavenumbers;				//number of wavenumbers
  float * wavenumber;				//array of wavenumbers
  int kii;					//index over wavenumbers
  int rii;					//index over radius
  float r;					//mean radius of present shell
  int timeii;					//index over time
  float rho;					//mean density
  
  n_wavenumbers = 100;				//number of wavenumbers; this may eventually be an input value
  
  rho = system->show_rho();
  /*allocate memory for wavenumbers and array of fourier transformed data*/
  wavenumber = new float [n_wavenumbers];
  fourier = new float* [n_times];		
  for(timeii=0;timeii<n_times;timeii++)
  {
    fourier[timeii] = new float [n_wavenumbers];
  }
  
  /*Calculate wavenumbers*/
  for(kii=0;kii<n_wavenumbers;kii++)
  {
    wavenumber[kii] = PI/max_value*(2.*float(kii)+.5);		//calculate wavenumber corresponding to this index over wavenumbers
  }
  
  /*calculate spherically symmetric spacial fourier transform*/
  for(timeii=0;timeii<n_times;timeii++)
  {
    for(kii=0;kii<n_wavenumbers;kii++)
    {
      fourier[timeii][kii]=1/wavenumber[kii];			//initialize value of fourier transform at this wavenumber to zero
      
      for(rii=0;rii<n_bins;rii++)		//sum over radius
      {
        r = min_value+(rii+0.5)*bin_size;	//calculate mean radius of current shell
	fourier[timeii][kii] += r * sin(wavenumber[kii] * r) * (normal[timeii][rii]/rho-1.0) * bin_size;	
      }
      fourier[timeii][kii] *= 4*PI*rho/wavenumber[kii];
    }
  }
  
    ofstream output (filename.c_str());		//open correlation file
  
  /*Write wavenumbers to file*/
  output << "\t";
  for(kii=0;kii<n_wavenumbers-1;kii++)
  {
    output << wavenumber[kii] << "\t";
  }
  output << "\n";
  
  /*Write fourier transform data to file, with each line corresponding to a time-spacing and led off by that time-spacing */
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii] << "\t";
    for(kii=0;kii<n_wavenumbers;kii++)
    {
      output << fourier[timeii][kii] << "\t";
    }
    output << "\n";
  }
  output.close();
  
  return fourier;
  #endif
  
  float ** normal;
  normal = normalized();
  double * fft_in;
  fftw_complex *fft_out;
  fftw_plan plan;
  float ** transform;
  int binii;
  int tii;
  int transformsize = n_bins-minbin;
  double r;
  double rmax;
  
  rmax = min_value + (n_bins+0.5)*bin_size;
  
  transform = new float * [n_times];
  fft_in = new double[transformsize];

  fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (floor(transformsize/2)+1));
  
  plan = fftw_plan_dft_r2c_1d(transformsize, fft_in, fft_out,FFTW_MEASURE);
  
  for(tii=0;tii<n_times;tii++)
  {
    for(binii=0;binii<transformsize;binii++)
    {
      r = min_value+(binii+0.5)*bin_size;	//calculate mean radius of current shell
      fft_in[binii] = double(r*normal[tii][binii+minbin]);
    }
    
    fftw_execute(plan);
    
    transform[tii]=new float [int(transformsize/2)-1];
    
    for(binii=0;binii<(int(transformsize/2)-1);binii++)
    {
      transform[tii][binii] = -4*PI*fft_out[binii][1]/((binii+1)*(2*PI/rmax)); 
    }
  }
  //fftw_destroy_plan(plan);
  
  ofstream output (filename.c_str());		//open correlation file
   
  output << "\t";
  for(binii=0;binii<(int(transformsize/2)-1);binii++)
  {
    output << (binii+1)*(2*PI/rmax) << "\t";
  }
  output << "\n";
  
  for(tii=0;tii<n_times;tii++)
  {
    output << timetable[tii] << "\t";
    for(binii=0;binii<(int(transformsize/2)-1);binii++)
    {
      output << transform[tii][binii] << "\t";		//write bins at this time to file
    }
  output << "\n";
  }
  output.close();
  
  return transform;
}
#endif


/*------------------------------------------------------------------------------*/





void Space-Time_Correlation_Function::bin(int timestep, float distance)
 {
  int binindex;
  binindex = int((distance-min_value)/bin_size);
    
  if(binindex>=0)
  {
    if(binindex<n_bins)
    {
      (correlation[timestep][binindex])++;
    }
  }
  //(weighting[timestep])++;  //putting the weighting out here (so that the weighting considers even particles falling outside the bin range) ensures that every timestep is normalized in the same way (maybe).
  //if(timestep!=1){cout << timestep << "\t";}
  
}



/*------------------------------------------------------------------------------*/




void Space-Time_Correlation_Function::write(string filename)const
{
  int timeii;
  int binii;
  
  float ** normal;				//temporarily store normalized correlation
  
  normal = normalized();
  
  ofstream output (filename.c_str());		//open correlation file
  
  output << "Correlation data created by MDAT v." << VERSION << "\n"; 
  output << n_bins << " bins\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting correlation function to file " <<filename<<"." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << min_value+bin_size/2+float(binii)*bin_size << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << normal[timeii][binii] << "\t";		//write bins at this time to file
//      output << correlation[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << normal[timeii][n_bins-1] << "\n";		//write last bin at this time to file
//    output << correlation[timeii][n_bins-1] << "\n";		//write last bin at this time to file
    delete [] (normal[timeii]);				//deallocate temporary memory for normalized correlation at this time
  }
  delete [] normal;					//deallocate temporary memory for normalized correlation
  output.close();					//close file  exit(0);
}



void Space-Time_Correlation_Function::postprocess_list()
{
  int timeii;
  for (timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii]/=float(weighting[timeii]);
  }
  
}
