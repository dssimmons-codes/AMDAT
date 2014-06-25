/*Written by David S. Simmons*/

#include <complex>
#include <string>
#include "analysis.h"
#include "wave_vectors.h"
#include "system.h"


namespace std{

class Structure_Factor:public Analysis
{
    Trajectory_List * trajlist1;
    Trajectory_List * trajlist2;
    float * structure_factor;
    int n_wavenumbers;
    bool fullblock;
    
    int n_atoms;
    int currenttime;
    
    complex<float> ** wavedensity1;
    complex<float> ** wavedensity2;
    
    complex<float> ** current_wavedensity;
    
    void atomkernel(int,int,int,int){};
    
    Wave_Vectors const * wavevectors;
    
    
  public:
    Structure_Factor();
    ~Structure_Factor();
    Structure_Factor(const Structure_Factor &);
    Structure_Factor operator = (const Structure_Factor &);
    Structure_Factor(System * sys, const Wave_Vectors * wv, bool fullblock=1);
    
    void analyze(Trajectory_List * t_list1);
    void analyze(Trajectory_List * t_list1, Trajectory_List * t_list2);	
    
    void analyze_wave_density(Trajectory_List * t_list);
    void listkernel(Trajectory* current_trajectory);
    void write(string filename)const;
    bool isThreadSafe(){return true;};    
};

}
