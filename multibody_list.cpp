/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_List class*/
/*Written by David S. Simmons*/

#include <fstream>
#include <stdlib.h>
#include "multibody_list.h"
//include "multibody_analysis.h"
#include "system.h"

using namespace std;


/*Default constructor*/
Multibody_List::Multibody_List()
{
  sys=0;

  n_times=0;
  n_bodies=-2;

  time_conversion=new int [0];
  multibodies.resize(n_times);

}




/*Copy constructor*/
Multibody_List::Multibody_List(const Multibody_List & copy)
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


/*Equality operator*/
Multibody_List Multibody_List::operator=(const Multibody_List & copy)
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


/*Destructor*/
Multibody_List::~Multibody_List()
{
    delete [] time_conversion;
}

/*Constructor that provides system and number of times*/
Multibody_List::Multibody_List(System * syst,int timecount)
{
  int timeii;

  sys=syst;
  n_times=timecount;
  n_bodies=-2;

  time_conversion=new int[sys->show_n_timesteps()];
  multibodies.resize(n_times);

  for(timeii=0;timeii<n_times;timeii++)
  {
    multibodies[timeii].resize(0);
  }
}

/*Constructor that provides system and a pointer to a Multibody_Set - fully creates list with multibodies in it*/
/*This is comparable to a static trajectory list in that it defines only one list of multibodies that is the same for all times in the system*/
Multibody_List::Multibody_List(System * syst, Multibody_Set * multibodyset)
{
    int multibodyii;

    sys = syst;
    n_times = 1;
    multibodies.resize(n_times);

    time_conversion = new int[sys->show_n_timesteps()];
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=0;
    }

    multibodies[0].resize(multibodyset->show_n_multibodies());

    for(multibodyii=0;multibodyii<multibodyset->show_n_multibodies();multibodyii++)
    {
        multibodies[0][multibodyii]=multibodyset->show_multibody(multibodyii);
    }

    check_n_bodies();
}



/*Method that resets object with system and a pointer to a Multibody_Set - fully creates list with multibodies in it*/
/*This is comparable to a static trajectory list in that it defines only one list of multibodies that is the same for all times in the system*/
void Multibody_List::set(System * syst, Multibody_Set * multibodyset)
{
    int multibodyii;
    
    delete [] time_conversion;
    
    sys = syst;
    n_times = 1;
    multibodies.resize(n_times);

    time_conversion = new int[sys->show_n_timesteps()];
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=0;
    }

    multibodies[0].resize(multibodyset->show_n_multibodies());

    for(multibodyii=0;multibodyii<multibodyset->show_n_multibodies();multibodyii++)
    {
        multibodies[0][multibodyii]=multibodyset->show_multibody(multibodyii);
    }

    check_n_bodies();
}



Multibody_List Multibody_List::operator+(const Multibody_List &)
{

}

/*Determines whether all multibodies in the list have the same number of bodies.*/
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


void Multibody_List::listloop(Multibody_Analysis* analysis, int timegap, int currenttime, int nexttime)
{
	int internal_time;

	current_time=convert_time(current_time);
	for(int trajectoryii=0;trajectoryii<multibodies[internal_time].size();trajectoryii++) 	
	{
		analysis->listkernel(multibodies[internal_time][trajectoryii], timegap, currentime, nexttime);
	}
}
