/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Header for Multibody_Set class - Class to store sets of user-defined multibodies*/
/*Written by David S. Simmons*/

#ifndef MULTIBODY_SET
#define MULTIBODY_SET

#include <vector>

#include "multibody.h"
#include "trajectory.h"

namespace std
{
class Multibody_Set
{
    int n_multibodies;
    vector <Multibody> multibodies;

  public:
    Multibody_Set();
    Multibody_Set(const Multibody_Set &);
    ~Multibody_Set();
    Multibody_Set operator =(const Multibody_Set &);
    Multibody_Set(int multibody_count);

    void add_multibody(const Multibody & multibody);

    /*Methods to define a multibody in the array of multibodies*/
    void set_multibody(int multibody_index, int n_bodies, Trajectory** bodies){multibodies[multibody_index].set(n_bodies,bodies);};
    void set_multibody(int multibody_index, const Multibody & multibody){multibodies[multibody_index]=multibody;};
    Multibody * show_multibody(int index);
    int show_n_multibodies(){return n_multibodies;};




};

}

#endif
