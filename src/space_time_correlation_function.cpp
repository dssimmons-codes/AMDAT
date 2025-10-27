// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Space_Time_Correlation_Function class: a general class to hold a time correlation function*/
/*Also calculates and stores inverse-space */
/*Written by David S. Simmons*/


#include "space_time_correlation_function.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "generated/version.h"

using namespace std;

#include "progress.h"

#ifndef PI
#define PI 3.141592653589793238462643383280
#endif




Space_Time_Correlation_Function::~Space_Time_Correlation_Function()
{
  //clear_memory();
}



void Space_Time_Correlation_Function::clear_memory()
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



Space_Time_Correlation_Function Space_Time_Correlation_Function::operator+ (const Space_Time_Correlation_Function & increment)const
{
  int timeii, binii;
  Space_Time_Correlation_Function temp;

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

  temp.system = system;		//naively copy system point from present object; the user must simply be smart about what this means and when it is appropriate.  In general, it is acceptable if systems have the same timestep properties and merely have different trajectories.  Could put a warning in here triggered if time properties are different.
  temp.correlation = new float*[n_times];
  temp.weighting = new int[n_times];
  temp.timetable = new float[n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    temp.correlation[timeii]=new float[n_bins];
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



// Space_Time_Correlation_Function Space_Time_Correlation_Function::operator= (Space_Time_Correlation_Function increment)
// {
//   int timeii, binii;
//   Space_Time_Correlation_Function temp;
//   temp.n_bins=n_bins;
//   temp.n_times=n_times;
//   temp.bin_size=bin_size;
//   temp.max_value=max_value;
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

/*sets the spatial inverse correlation equal to the spherical fourier transform of the correlation function*/

void Space_Time_Correlation_Function::spherical_fourier()
{
  //float ** fourier;				//define pointer to array of fourier transformed data
  
  float * wavenumber;				//array of wavenumbers
  int kii;					//index over wavenumbers
  int rii;					//index over radius
  float r;					//mean radius of present shell
  int timeii;					//index over time
  float rho;					//mean density
  
  //n_wavenumbers = 100;				//number of wavenumbers; this may eventually be an input value
  
  //rho = system->show_rho();
  /*allocate memory for wavenumbers and array of fourier transformed data*/
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] spatial_inverse[timeii];
  }
  delete [] spatial_inverse;
  
  wavenumber = new float [n_wavenumbers];
  spatial_inverse = new float* [n_times];		
  for(timeii=0;timeii<n_times;timeii++)
  {
    spatial_inverse[timeii] = new float [n_wavenumbers];
  }
  
  /*Calculate wavenumbers*/
  for(kii=0;kii<n_wavenumbers;kii++)
  {
    wavenumber[kii] = 2.*PI*float(kii+1)/max_value;		//calculate wavenumber corresponding to this index over wavenumbers
   // wavenumber[kii] = PI/max_value*(2.*float(kii)+.5);		//calculate wavenumber corresponding to this index over wavenumbers
  }
  
  /*calculate spherically symmetric spacial fourier transform*/
  for(timeii=0;timeii<n_times;timeii++)
  {
    for(kii=0;kii<n_wavenumbers;kii++)
    {
      spatial_inverse[timeii][kii]=0;			//initialize value of fourier transform at this wavenumber to zero
      
      for(rii=0;rii<n_bins;rii++)		//sum over radius
      {
        r = (rii+0.5)*bin_size;	//calculate mean radius of current shell
        spatial_inverse[timeii][kii] += r * correlation[timeii][rii]* sin(wavenumber[kii] * r);
	//fourier[timeii][kii] += r * sin(wavenumber[kii] * r) * (normal[timeii][rii]/rho-1.0) * bin_size;	
      }
      spatial_inverse[timeii][kii] *= 4*PI*rho/wavenumber[kii];
    }
  }
  
}


/*default method to set the spatial inverse to be the fourier transform of the real-space correlation function*/
void Space_Time_Correlation_Function::calculate_spatial_inverse(int n_wavenums)
{
  n_wavenumbers = n_wavenums;
  spherical_fourier();
  spatial_inverse_calculated=true;
}


/*------------------------------------------------------------------------------*/





void Space_Time_Correlation_Function::bin(int timestep, float distance)
 {
  int binindex;
  binindex = int((distance)/bin_size);
    
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




void Space_Time_Correlation_Function::write(string filename)const
{
  int timeii;
  int binii;
  
  ofstream output (filename.c_str());		//open correlation file
  
  output << "Correlation data created by AMDAT v." << amdat::build::SEMVER << "\n"; 
  output << n_bins << " bins\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting correlation function to file " <<filename<<"." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << bin_size/2+float(binii)*bin_size << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << correlation[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << correlation[timeii][n_bins-1] << "\n";		//write last bin at this time to file
  }
  output.close();					//close file  exit(0);
}


void Space_Time_Correlation_Function::write(ofstream& output)const
{
  int timeii;
  int binii;
  
  output << "Correlation data created by AMDAT v." << amdat::build::SEMVER << "\n"; 
  output << n_bins << " bins\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting correlation function to file." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << bin_size/2+float(binii)*bin_size << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << correlation[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << correlation[timeii][n_bins-1] << "\n";		//write last bin at this time to file
  }
}


/*------------------------------------------------------------------------------*/




void Space_Time_Correlation_Function::write_spatial_inverse(string filename)const
{
  int timeii;
  int binii;
  
  ofstream output (filename.c_str());		//open correlation file
  
  output << "Inverse-space correlation function data created by AMDAT v." << amdat::build::SEMVER << "\n"; 
  output << n_wavenumbers << " wavenumbers\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting inverse-space correlation function to file " <<filename<<"." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << 2*PI*(binii+1)/max_value << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << spatial_inverse[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << spatial_inverse[timeii][n_bins-1] << "\n";		//write last bin at this time to file
  }
  output.close();					//close file  exit(0);
}



void Space_Time_Correlation_Function::postprocess_list()
{
  int timeii,binii;
  float rshell, shellvolume;
  
  for(binii=0;binii<n_bins;binii++)
  {
    rshell = binii*bin_size;						//determine inner radius of bin
    shellvolume = (4.0/3.0)*PI*(pow(rshell+bin_size,3.0)-pow(rshell,3.0));		//calculate volume of bin
    for (timeii=0;timeii<n_times;timeii++)
    {
      correlation[timeii][binii]/=(float(weighting[timeii])*shellvolume);
    }
  }
  
}
