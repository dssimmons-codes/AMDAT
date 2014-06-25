/*Written by David S. Simmons*/

#include "correlation_2d.h"
#include "wave_vectors.h"

namespace std{

class Incoherent_Scattering_Function:public Correlation_2D
{
    /*computational members - just used in calculations; no useful value later on*/
    int currenttime, nexttime, timegap;

    bool fullblock;
    
  public:
    Incoherent_Scattering_Function();		//default constructor
    Incoherent_Scattering_Function(const Incoherent_Scattering_Function &);		//copy constructor
    Incoherent_Scattering_Function operator=(const Incoherent_Scattering_Function &);
    //~Incoherent_Scattering_Function();
    
    Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, bool fblock=0);
    Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, int inner, int outer, bool fblock = 0);
    
    Analysis_Type what_are_you(){Analysis_Type type = incoherent_scattering_function; return type;};		//virtual method to report the type of analysis
    
    
     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void listkernel(Trajectory *, int, int, int);
    
    void bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii);
    void postprocess_bins();
    //bool isThreadSafe(){return true;};
};

}
