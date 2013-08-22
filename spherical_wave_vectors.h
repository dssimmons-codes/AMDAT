/*Spherical_Wave_Vectors class.  This class bins a grind of wave vectors into spherical shells of wavenumber.*/
/*Written by David S. Simmons*/

#include "wave_vectors.h"

namespace std{

class Spherical_Wave_Vectors:public Wave_Vectors
{
  private:
    void bin (int xii, int yii, int zii);
  public:
    Spherical_Wave_Vectors(System* sys,int shellcount=300){calculate(sys,shellcount);};
};

}