/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Trajectory_RgTensor class: calculates time-dependent Rg Tensor for a trajectory*/
/*Written by David S. Simmons*/

/*Uses matrix math methods from the Template Numerical Toolkit (TNT)*/
/*TNT can be found at math.nist.gov/tnt/ */

/*need to incorporate center of mass part*/

#include "rgtensor.h"
#include "math.h"
#include <jama_eig.h>
#include <stdlib.h>
#include "version.h"

using namespace std;

RgTensor::RgTensor(System* sys,Trajectory*traj)
{
  int timeii;
  
  system=sys;
  trajectory=traj;
  
  n_blocks = system->show_n_exponentials();
  
  rgtensor = new sixfloat * [n_blocks];
  maxtimes = system->show_n_timegaps();
  
  blocks_per_time = new int [maxtimes];
  mean_eigenvalues = new threefloat [maxtimes];
  
  blocksize = new int [n_blocks];
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    blocks_per_time[timeii]=0;
    mean_eigenvalues[timeii][0]=0;
    mean_eigenvalues[timeii][1]=0;
    mean_eigenvalues[timeii][2]=0;
  }
  
  
  calc_tensor();
  calc_eigenvalues();
}


void RgTensor::calc_tensor()
{
  int blockii;
  int * timelist;
  float * times;
  
  sixfloat * temptensor;
  
  Coordinate * coordinatelist;
  
  for(blockii=0;blockii<n_blocks;blockii++)
  {
    
    blocksize[blockii] = system->big_blocksize(blockii);
    timelist = new int [blocksize[blockii]];
    times = new float [blocksize[blockii]];
    system->big_block(blockii,timelist);	//list of times for block (counting allowed inter-block timegaps)
    system->show_times(blocksize[blockii],timelist,times);
    
    temptensor = new sixfloat [blocksize[blockii]];
    coordinatelist=trajectory->show_coordinates(timelist,blocksize[blockii]);	//get list of coordinates	
    gyration_tensor(coordinatelist,times,blocksize[blockii],temptensor);		//calculate gyration tensor as function of time
    rgtensor[blockii]=temptensor;			//copy gyration tensor into long-term memory
  }
  
  
}

void RgTensor::calc_eigenvalues()
{
  int blockii,timeii;
  TNT::Array2D<float> rgarray(3,3,0.0);
  eigenvalues = new threefloat * [n_blocks];
  for(blockii=0;blockii<n_blocks;blockii++)
  {
    eigenvalues[blockii] = new threefloat [blocksize[blockii]];
    for(timeii=0;timeii<blocksize[blockii];timeii++)
    {
      /*create TNT Array2D corresponding to current gyration tensor*/
      rgarray[0][0]=rgtensor[blockii][timeii][0];
      rgarray[1][0]=rgtensor[blockii][timeii][1];
      rgarray[2][0]=rgtensor[blockii][timeii][2];
      rgarray[0][1]=rgtensor[blockii][timeii][1];
      rgarray[1][1]=rgtensor[blockii][timeii][3];
      rgarray[2][1]=rgtensor[blockii][timeii][4];
      rgarray[0][2]=rgtensor[blockii][timeii][2];
      rgarray[1][2]=rgtensor[blockii][timeii][4];
      rgarray[2][2]=rgtensor[blockii][timeii][5];
      
      eigkernel(rgarray,blockii,timeii);
      
      mean_eigenvalues[timeii][0]+=eigenvalues[blockii][timeii][0];
      mean_eigenvalues[timeii][1]+=eigenvalues[blockii][timeii][1];
      mean_eigenvalues[timeii][2]+=eigenvalues[blockii][timeii][2];
      
      blocks_per_time[timeii]++;
    }
  }
  
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    mean_eigenvalues[timeii][0]/=blocks_per_time[timeii];
    mean_eigenvalues[timeii][1]/=blocks_per_time[timeii];
    mean_eigenvalues[timeii][2]/=blocks_per_time[timeii];
  }
}

void RgTensor::eigkernel(TNT::Array2D<float> rgarray, int blockii, int timeii)
{
  TNT::Array1D<float> eigvals;
  
  JAMA::Eigenvalue <float> eigmat(rgarray);
  eigmat.getRealEigenvalues(eigvals);
      
  eigenvalues[blockii][timeii][0]=eigvals[0];
  eigenvalues[blockii][timeii][1]=eigvals[1];
  eigenvalues[blockii][timeii][2]=eigvals[2];
  
}


void RgTensor::write(string filename)
{
  int timeii;
  float * times;
  
  ofstream output(filename.c_str());
  output << "Single particle gyration tensor data created by MDAT v." << VERSION << "\n"; 
  
  times = system->displacement_times();
  
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    //output << times[timeii] << "\t" << mean_eigenvalues[timeii][0] << "\t" << mean_eigenvalues[timeii][1] << "\t" << mean_eigenvalues[timeii][2] << "\n";
    output << times[timeii] << "\t" << eigenvalues[0][timeii][0] << "\t" << eigenvalues[0][timeii][1] << "\t" << eigenvalues[0][timeii][2] << "\n";
  }
}