


#ifndef SIZE_STATISTICS
#define SIZE_STATISTICS

#include "multibody_analysis.h"

namespace std{

class Size_Statistics: public Multibody_Analysis
{
    vector<float> size_count;
    int weighing;
    float * moments;
    int maxmoment;
    

  public:
    Size_Statistics();
    Size_Statistics(const Size_Statistics &);
    Size_Statistics operator=(const Size_Statistics &);
    ~Size_Statistics(){};
    
    Size_Statistics(System*);
    
    void set(System*)
    
    void analyze(Multibody_List * mblist);
    void listkernel(Multibody *, int, int, int);
    void postprocess();
    
    void write(string);
    void write(ofstream&)const;
  
};
}