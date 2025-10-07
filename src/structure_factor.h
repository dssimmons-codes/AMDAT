// Copyright (C) 2010-2025 David S. Simmons

/*Written by David S. Simmons*/

#include <complex>
#include <string>
#include "analysis.h"
#include "wave_vectors.h"
#include "system.h"
#include "analysis_onetime.h"

namespace std{

class Structure_Factor:public Analysis_Onetime
{
    //Trajectory_List * trajlist1;
    //Trajectory_List * trajlist2;
    float * structure_factor;
    int n_wavenumbers;
    
    int n_atoms;
    int currenttime;
    
    complex<float> ** wavedensity1;
    complex<float> ** wavedensity2;
    
    complex<float> ** current_wavedensity;
    
    void atomkernel(int,int,int,int){};
    
    Wave_Vectors const * wavevectors;
    
    void preprocess();
    void preprocess2();
    void timekernel (int timeii);	//method that is looped over by analyze with single trajectory list
    void timekernel2 (int timeii);	//method that is looped over by analyze with two trajectory lists
    void postprocess_list();
    
    
  public:
    Structure_Factor();
    ~Structure_Factor();
    Structure_Factor(const Structure_Factor &);
    Structure_Factor operator = (const Structure_Factor &);
    Structure_Factor(System * sys, const Wave_Vectors * wv, int timescheme = -1);
    
    void analyze_wave_density(Trajectory_List * t_list);
    void listkernel(Trajectory* current_trajectory);
    void write(string filename)const;
    void write(ofstream& output)const;
    //bool isThreadSafe(){return true;};    
};

}
