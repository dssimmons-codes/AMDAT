/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_List class*/
/*Written by David S. Simmons*/

#include <fstream>
#include <stdlib.h>
#include "multibody_list.h"
#include "analysis.h"
#include "multibody_analysis.h"
#include "system.h"

using namespace std;

Multibody_List::Multibody_List();
{
  sys=0;

  n_times=0;
  n_bodies=-2;

  time_conversion.resize(n_times);
  multibodies.resize(n_times);

}





Multibody_List::Multibody_List(const Multibody_List & copy);
{
  sys=copy.sys;

  n_times=copy.n_times;
  n_bodies=copy.n_bodies;

  time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
  multibodies=copy.multibodies;

}



Multibody_List Multibody_List::operator=(const Multibody_List & copy);
{
  if(this!=&copy)
  {

    delete [] time_conversion;
    sys=copy.sys;

    n_times=copy.n_times;
    n_bodies=copy.n_bodies;

     time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
    multibodies=copy.multibodies;
  }

  return *this;
}



Multibody_List::~Multibody_List()
{
    delete [] time_conversion;
}


Multibody_List::Multibody_List(System sys,int timecount);
{
  int timeii;

  sys=copy.sys;
  n_times=0;
  n_bodies=-2;

  time_conversion=new int[sys->show_n_timesteps()];
  multibodies.resize(n_times);

  for(timeii=0;timeii<n_times;timeii++)
  {
    multibodies[timeii].resize(0);
  }
}


Multibody_List::Multibody_List(System syst, Multibody_Set multibodyset)
{
    int multibodyii;

    sys = syst;
    n_times = 1;
    multibodies.resize(n_times);

    time_conversion = new int[sys->show_n_timesteps()]
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=0;
    }

    multibodies[0].resize(multibodyset.show_n_multibodies());

    for(multibodyii=0;multibodyii<multibodyset.show_n_multibodies();multibodyii++)
    {
        multibodies[0][multibodyii]=show_multibody(multibodyii);
    }

    check_n_bodies();
}


Multibody_List Multibody_List::operator+(const Multibody_List)
{

}


void Multibody_List::check_n_bodies()
{
    int timeii, multibodyii;

    n_bodies = multibodies[0][0]->show_n_bodies();


    for(timeii=0;timeii<multibodies.size();timeii++)
    {
        for(multibodyii=0;multibodyii<multibodies[timeii].size();multibodyii++)
        {
            if(n_bodies!=multibodies[timeii][multibodyii]->show_n_bodies())
            {
                n_bodies = -1;
            }
        }
    }
}
