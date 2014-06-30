/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_Set class - Class to store sets of user-defined multibodies
/*Written by David S. Simmons*/

#include "multibody_set.h"

using namespace std;

/*Default constructor*/
Multibody_Set::Multibody_Set()
{
  n_multibodies=0;
  multibodies = new Multibody [n_multibodies];
}


/*Copy constructor*/
Multibody_Set::Multibody_Set(const Multibody_Set & copy)
{
  int multibodyii;

  n_multibodies=copy.n_multibodies;
  multibodies = new Multibody [n_multibodies];

  for(multibodyii=0;multibodyii<n_multibodies;multibodyii++)
  {
    multibodies[multibodyii]=copy.multibodies[multibodyii];
  }
}



/*Destructor*/
Multibody_Set::~Multibody_Set()
{

}



/*Equality operator*/
Multibody_Set Multibody_Set::operator=(const Multibody_Set & copy)
{
  int multibodyii;

  if (this!=&copy)
  {
    int multibodyii;
    delete [] multibodies;
    n_multibodies=copy.n_multibodies;
    multibodies = new Multibody [n_multibodies];

    for(multibodyii=0;multibodyii<n_multibodies;multibodyii++)
    {
      multibodies[multibodyii]=copy.multibodies[multibodyii];
    }
  }
  return *this;
}


/*Constructor that sets number of multibodies*/
Multibody_Set::Multibody_Set(int multibody_count)
{
  n_multibodies=multibody_count;
  multibodies = new Multibody [n_multibodies];
}




/*reset number of multibodies and reinitialize array of multibodies*/
void Multibody_Set::set(int multibody_count)
{
  delete [] multibodies;
  n_multibodies=multibody_count;
  multibodies = new Multibody [n_multibodies];
}


Multibody* Multibody_Set::show_multibody(int index)
{
    if(index<n_multibodies)
    {
        return &multibodies[index];
    }
    else
    {
        return 0;
    };
}
