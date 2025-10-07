/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_Region class - stores list of multibody objects that are within a given region*/
/*Written by David S. Simmons*/


#ifndef MULTIBODY_BEAD_REGION
#define MULTIBODY_BEAD_REGION

#include "coordinate.h"
#include "multibody_list.h"
#include "multibody_analysis.h"

namespace std{

class Multibody_Bead_Region: public Multibody_List, public Multibody_Analysis
{
  private:

    Coordinate upperbound;
    Coordinate lowerbound;
    int threshold;

  public:


    Multibody_Bead_Region();
    Multibody_Bead_Region(System*syst);
    Multibody_Bead_Region(System*syst, Coordinate low, Coordinate high, int thresh);
    Multibody_Bead_Region(const Multibody_Bead_Region &);
    Multibody_Bead_Region operator=(const Multibody_Bead_Region &);

    void analyze(Multibody_List * mblist);
    void list_displacementkerneyl(int,int,int){};
    void listkernel(Multibody*,int,int,int);
    void postprocess_list();

    void write(string) const;

};

}


#endif
