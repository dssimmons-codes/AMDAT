/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to compare self Van Hove function to the so-called Gaussian approximation for the Van Hove, establishing boundaries and fractions for slow and fast particle categories*/
/*Written by David S. Simmons*/


#include "gaussian_comparison.h"
#include "math.h"
#include "mean_square_displacement.h"
#include <iostream>
#include <stdlib.h>
#include "version.h"

#define PI 3.14159265

using namespace std;

Gaussian_Comparison::Gaussian_Comparison()
{
  int binii;

  system = 0;
  self_van_hove = 0;
  non_gaussian_parameter = 0;
  
  n_bins = 0;
  bin_size = 0;
  
  time_index = 0;
  mean_square_displacement = 0;
  
  gaussian_approx = new float [n_bins];
  van_hove = new float [n_bins];
  
  for(binii=0;binii<n_bins;binii++)
  {
    gaussian_approx[binii] = van_hove[binii] = 0;
  }
  
  slowboundary_index = 0;
  fastboundary_index = 0;
  slowfraction = fastfraction = gaussian_slowfraction = gaussian_fastfraction = 0;
  
  errorstate=0;
}


Gaussian_Comparison::Gaussian_Comparison(System* sys, const Non_Gaussian_Parameter* ngp, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd)
{
  int binii;
  system = sys;
  self_van_hove = vhs;
  non_gaussian_parameter = ngp;
  
  n_bins = self_van_hove->show_n_bins();
  bin_size = self_van_hove->show_bin_size();
  
  time_index = non_gaussian_parameter->max();
  
  mean_square_displacement = msd->show(time_index);

  gaussian_approx = new float [n_bins];
  van_hove = new float [n_bins];
  
  for(binii=0;binii<n_bins;binii++)
  {
    gaussian_approx[binii] = van_hove[binii] = 0;
  }
  
  errorstate=0;
  
  slowboundary_index = 0;
  fastboundary_index = 0;
  slowfraction = fastfraction = gaussian_slowfraction = gaussian_fastfraction = 0;
  
  scale_van_hoves();
  find_crossovers();
  determine_fractions();
}


Gaussian_Comparison::Gaussian_Comparison(System* sys, int time_in, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd)
{
	int binii;
	system = sys;
	self_van_hove = vhs;
  
	n_bins = self_van_hove->show_n_bins();
	bin_size = self_van_hove->show_bin_size();
  
	time_index = time_in;
  
	mean_square_displacement = msd->show(time_index);

	gaussian_approx = new float [n_bins];
	van_hove = new float [n_bins];
  
	for(binii=0;binii<n_bins;binii++)
	{
		gaussian_approx[binii] = van_hove[binii] = 0;
	}
  
	errorstate=0;
  
	slowboundary_index = 0;
	fastboundary_index = 0;
	slowfraction = fastfraction = gaussian_slowfraction = gaussian_fastfraction = 0;
  
	scale_van_hoves();
	find_crossovers();
	determine_fractions();
}


void Gaussian_Comparison::clear()
{
  delete [] gaussian_approx;
  delete [] van_hove;
}



void Gaussian_Comparison::set(System* sys, const Non_Gaussian_Parameter* ngp, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd)
{
  int binii;
  
  clear();
  
  system = sys;
  self_van_hove = vhs;
  non_gaussian_parameter = ngp;
  
  n_bins = self_van_hove->show_n_bins();
  bin_size = self_van_hove->show_bin_size();
  
  time_index = non_gaussian_parameter->max();
  
  mean_square_displacement = msd->show(time_index);

  gaussian_approx = new float [n_bins];
  van_hove = new float [n_bins];
  
  for(binii=0;binii<n_bins;binii++)
  {
    gaussian_approx[binii] = van_hove[binii] = 0;
  }
   errorstate=0;
    
  slowboundary_index = 0;
  fastboundary_index = 0;
  slowfraction = fastfraction = gaussian_slowfraction = gaussian_fastfraction = 0;
  
  scale_van_hoves();
  find_crossovers();
  determine_fractions();
}


/*Method to intital gaussian comparison with time selected by user*/
void Gaussian_Comparison::set(System* sys, int time_in, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd)
{
	int binii;
  
	clear();
  
	system = sys;
	self_van_hove = vhs;
  
	n_bins = self_van_hove->show_n_bins();
	bin_size = self_van_hove->show_bin_size();
  
	time_index = time_in;
  
	mean_square_displacement = msd->show(time_index);

	gaussian_approx = new float [n_bins];
	van_hove = new float [n_bins];
  
	for(binii=0;binii<n_bins;binii++)
	{
		gaussian_approx[binii] = van_hove[binii] = 0;
	}
	errorstate=0;
    
	slowboundary_index = 0;
	fastboundary_index = 0;
	slowfraction = fastfraction = gaussian_slowfraction = gaussian_fastfraction = 0;
  
	scale_van_hoves();
	find_crossovers();
	determine_fractions();
}


/*Method to generate van hove and gaussian van hove (self) that are scaled by shell volume*/
void Gaussian_Comparison::scale_van_hoves()
{
  float alpha;		//scaling factor of gaussian
  int binii;
  float shellvolume;

  alpha = 3.0/(2.0*mean_square_displacement);
  for(binii=0;binii<n_bins;binii++)
  {
    shellvolume = 4.0/3.0*PI*(pow((float(binii)+1.0)*bin_size,3.0)-pow(float(binii)*bin_size,3.0));		//calculate shell volume
    van_hove[binii] = self_van_hove->show(time_index,binii)*shellvolume; 	//scale self van hove by shell volume
    gaussian_approx[binii] = pow(alpha/PI,3.0/2.0)*exp(-alpha*pow((float(binii)+.5)*bin_size,2.0))*shellvolume;  //generate array of gaussian approximation van hove, scaled by shell volume
  }

}



/*Method to find points where self van hove and gaussian self van hove cross*/
void Gaussian_Comparison::find_crossovers()
{
  int binii;

  for(binii=1;binii<n_bins;binii++)
  {
    if(gaussian_approx[binii]>=van_hove[binii])
    {
      if(binii==0)
        {slowboundary_index = 0;}
      else
        {slowboundary_index = binii-1;}
      break;
    }
  }
  for(binii=slowboundary_index+1;binii<n_bins;binii++)
  {
    if(gaussian_approx[binii]<=van_hove[binii])
    {
      fastboundary_index = binii;
      break;
    }
  }
  if(fastboundary_index==0)
  {
    cout << "\nWarning: program cannot identify two crossing points between self Van Hove and Gaussian approximation self Van Hove.  Fractions of slow and fast particles will be incorrect.\n";
    errorstate = bool(1);
  }
  
  slowboundary = pow((float(slowboundary_index)+.5)*bin_size,2);
  fastboundary = pow((float(fastboundary_index)-.5)*bin_size,2);
  
}





/*Method to determine fractions of slow and fast particles*/
/*presently this uses the points outside the gaussian approximation curve as the boundary.  A better method, possibly to be implemented later, could do a linear interpolation to find a more accurate crossing point and weight the edge bins based on this determination.*/
void Gaussian_Comparison::determine_fractions()
{
  int binii;

  for(binii=1;binii<=slowboundary_index;binii++)
  {
    slowfraction += van_hove[binii];
    gaussian_slowfraction += gaussian_approx[binii];
  }
  
  fastfraction = 1.0-slowfraction;
  gaussian_fastfraction = 1.0-gaussian_slowfraction;
  for(binii=slowboundary_index+1;binii<fastboundary_index;binii++)
  {
    fastfraction  -= van_hove[binii];
    gaussian_fastfraction -= gaussian_approx[binii];
  }
  
}



void Gaussian_Comparison::write(string filename)const
{
  int binii;

  cout << "\nWriting fractions of slow and fast particles and comparison between gaussian and actual self van hove to file " << filename <<".";
  
  ofstream output(filename.c_str());
  if(errorstate){output << "Warning: program cannot identify two crossing points between self Van Hove and Gaussian approximation self Van Hove.  Fractions of slow and fast particles will be incorrect.\n";}
  
  output << "Gaussian comparison data created by MDAT v." << VERSION << "\n"; 
  output << "t*\n" << system->displacement_times()[time_index] << "\n";
  output << "Slow_Cutoff: " << pow(slowboundary,0.5) << "\n";
  output << "Fast_Cutoff: " << pow(fastboundary,0.5) << "\n";
  output << "\t" << "Slow_Fraction"<<"\t"<<"Fast_Fraction\n";
  output << "Gaussian" << "\t" << gaussian_slowfraction << "\t" << gaussian_fastfraction<<"\n";
  output << "Actual" << "\t" << slowfraction << "\t" << fastfraction << "\n";
  output << "Excess" << "\t" << slowfraction - gaussian_slowfraction << "\t" << fastfraction - gaussian_fastfraction << "\n\n";
  
  output << "Self_Van_Hove*shellvolume\n";
  
  output << "radius\tgaussian\tactual\n";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << (binii+.5)*bin_size << "\t" << gaussian_approx[binii] << "\t" << van_hove[binii] << "\n";
  }
  output.close();
}

void Gaussian_Comparison::write(ofstream& output)const
{
  int binii;

  cout << "\nWriting fractions of slow and fast particles and comparison between gaussian and actual self van hove to file.";
  
  if(errorstate){output << "Warning: program cannot identify two crossing points between self Van Hove and Gaussian approximation self Van Hove.  Fractions of slow and fast particles will be incorrect.\n";}
  
  output << "Gaussian comparison data created by MDAT v." << VERSION << "\n"; 
  output << "t*\n" << system->displacement_times()[time_index] << "\n";
  output << "Slow_Cutoff: " << pow(slowboundary,0.5) << "\n";
  output << "Fast_Cutoff: " << pow(fastboundary,0.5) << "\n";
  output << "\t" << "Slow_Fraction"<<"\t"<<"Fast_Fraction\n";
  output << "Gaussian" << "\t" << gaussian_slowfraction << "\t" << gaussian_fastfraction<<"\n";
  output << "Actual" << "\t" << slowfraction << "\t" << fastfraction << "\n";
  output << "Excess" << "\t" << slowfraction - gaussian_slowfraction << "\t" << fastfraction - gaussian_fastfraction << "\n\n";
  
  output << "Self_Van_Hove*shellvolume\n";
  
  output << "radius\tgaussian\tactual\n";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << (binii+.5)*bin_size << "\t" << gaussian_approx[binii] << "\t" << van_hove[binii] << "\n";
  }
}