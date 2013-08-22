/*Written by David S. Simmons*/

#include "correlation_2d.h"
#include "wave_vectors.h"
#include "particle_list.h"

namespace std{

class Incoherent_Scattering_Function:public Correlation_2D
{
    /*computational members - just used in calculations; no useful value later on*/
    int currenttime, nexttime, timegap;
    
    void preprocess();
    
    bool fullblock;
    
  public:
    Incoherent_Scattering_Function();		//default constructor
    Incoherent_Scattering_Function(const Incoherent_Scattering_Function &);		//copy constructor
    Incoherent_Scattering_Function operator=(const Incoherent_Scattering_Function &);
    //~Incoherent_Scattering_Function();
    
    Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, bool fblock=0);
    Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, int inner, int outer, bool fblock = 0);
    
    void displacementkernel(int timegapii, int thisii, int nextii, Trajectory * traj);
    void atomkernel(Trajectory * traj);
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
	void listkernel(Trajectory *, int, int, int);
    
    void bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii);
    void postprocess_bins();
	bool isThreadSafe(){return true;};
};

}
