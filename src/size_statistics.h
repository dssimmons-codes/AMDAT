// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class Size_Statistics: computes size distributions and statistics for multibodies*/
/*Written by David S. Simmons*/


#ifndef SIZE_STATISTICS
#define SIZE_STATISTICS

#include "multibody_analysis.h"
#include "system.h"

namespace std{

class Size_Statistics: public Multibody_Analysis
{
    vector<float> size_count;
    int weighting;
    float * moments;
    int n_moments;
    

  public:
    Size_Statistics();
    Size_Statistics(const Size_Statistics &);
    Size_Statistics operator=(const Size_Statistics &);
    ~Size_Statistics();
    
    Size_Statistics(System*,int);
    
    void set(System*, int);
    
    void analyze(Multibody_List * mblist);
    void listkernel(Multibody *, int, int, int);
    void postprocess();
    
    void write(string)const;
    void write(ofstream&)const;
  
};
}

#endif