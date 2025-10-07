// Copyright (C) 2010-2025 David S. Simmons

#ifndef DISPLACEMENT_MAP
#define DISPLACEMENT_MAP

#include "analysis.h"
#include <string>
#include "value_list.h"

namespace std {

class Displacement_Map: public Value_List <float>, public Analysis
{
    int time_index;		//index of time-separation for which to calculate displacement
    
    float time;			//actual time-separation corresponding to msd
   
    void postprocess(){};

    int firstblock;
    int lastblock;

    int currenttime, nexttime, currentblock;
   
    float maxdisplacement;	//maximum displacement value to write out
   
  public:
    Displacement_Map(System* sys, int timespacing, int blockstart, int blockend, float maxdisp = 0);
    
    Analysis_Type what_are_you(){Analysis_Type type = displacement_map; return type;};		//virtual method to report the type of analysis
    
    void write(string)const;			//write distribution to file
    
    void analyze(Trajectory_List* t_list);
    void list_displacementkernel(int timegapii, int thisii, int nextii);
    void listkernel(Trajectory * current_trajectory);
    void postprocess_list(){};
    
};

}
#endif