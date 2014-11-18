/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Trajectory_RgTensor class: calculates time-dependent Rg Tensor for a trajectory*/
/*Written by David S. Simmons*/

#include "system.h"
#include "trajmath.h"
//#include <tnt.h>
#include "tntjama/tnt.h"

namespace std {

typedef float threefloat [3];
  
class RgTensor 
{    
  private:
    
    System* system;
    Trajectory * trajectory;		//trajectory for which the rg tensor is calculated
    sixfloat ** rgtensor;		//Rg tensor of trajectory as function of time, determined separately for the trajectory starting from the beginning of each block: rgtensor[blockii][timeii][index]
    int n_blocks; 
    int * blocksize;
    threefloat ** eigenvalues; 		//Eigenvalues of Rg tensor as a function of time, determined separately for the trajectory starting from the beginning of each block: eigenvalues[blockii][index];
    threefloat * mean_eigenvalues;
    int maxtimes;
    int * blocks_per_time;
    
    void calc_tensor();
    void calc_eigenvalues();
    
    void eigkernel(TNT::Array2D<float> rgarray, int blockii,int timeii);
  public:
    RgTensor(System* sys, Trajectory* traj);
    
    void reset(Trajectory*traj);
    
    void show_eigenvalues(int timeii,threefloat eigout){eigout=mean_eigenvalues[timeii];};
    float show_eigenvalues(int timeii,int ii){return mean_eigenvalues[timeii][ii];};
    void write(string);

    
};

}