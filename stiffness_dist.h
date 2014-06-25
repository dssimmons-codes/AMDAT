
#ifndef STIFFNESS_DIST
#define STIFFNESS_DIST

#include "analysis.h"
#include <string>

namespace std {

class Stiffness_Dist: public Analysis
{
    float * distribution;	//stores distribution in bins
    int n_bins;			//number of distribution bins
    float binsize;		//size of distribution bins
    int time_index;		//index of time-separation for which to calculate msd
    float time;			//actual time-separation corresponding to msd
    int atomcount;		//number of atoms used to produce this distribution (used for normalization)
    int timeweighting;		//normalization factor reflecting the number of timegaps used to produce distribution
    float mean;
    float variance;
    float square_term;
    float log_mean;
    float log_variance;
    float log_square_term;

   
    //calculation variables
    int currenttime, nexttime;

   
  public:
    Stiffness_Dist(System* sys, int bins, float maxvalue, float t);
    

    void write(string)const;			//write distribution to file
    
    void analyze(Trajectory_List* t_list);
    void list_displacementkernel(int timegapii, int thisii, int nextii);
    void listkernel(Trajectory * current_trajectory);
    void postprocess_list();
};

}

#endif