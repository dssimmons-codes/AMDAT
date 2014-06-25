/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate distribution of particle displacements to user specified power*/
/*Written by David S. Simmons*/

#ifndef DISPLACEMENT_DISTRIBUTION
#define DISPLACEMENT_DISTRIBUTION

#include <string>
#include <sstream>
#include "analysis.h"

namespace std {

class Displacement_Distribution: public Analysis
{
    float * distribution;	//stores distribution in bins
    int n_bins;			//number of distribution bins
    float binsize;		//size of distribution bins
    int time_index;		//index of time-separation for which to calculate msd
    float time;			//actual time-separation corresponding to msd
    int atomcount;		//number of atoms used to produce this distribution (used for normalization)
    int timeweighting;		//normalization factor reflecting the number of timegaps used to produce distribution
    float power;

    float mean;
    float variance;
    float square_term;

    int currenttime, nexttime;
   
   
  public:
    Displacement_Distribution(System* sys, float pow ,int bins, float maxvalue, float t);
    
    Analysis_Type what_are_you(){Analysis_Type type = displacement_dist; return type;};		//virtual method to report the type of analysis
    
    void write(string)const;			//write distribution to file
    
    void analyze(Trajectory_List* t_list);
    void list_displacementkernel(int timegapii, int thisii, int nextii);
    void listkernel(Trajectory * current_trajectory);
    void postprocess_list();
    
    void atomkernel(int,int, int, int){cout<<"Error: displacement_distribution does not work with system loops";};
};

}
#endif