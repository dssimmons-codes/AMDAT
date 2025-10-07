// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_List class*/
/*Written by David S. Simmons*/

#include <fstream>
#include <stdlib.h>
#include "multibody_list.h"
#include "multibody_analysis.h"
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

/*Method that resets object with system and a pointer to a list of Multibody_Set onjects - fully creates list with multibodies in it*/
Multibody_List::Multibody_List(System * syst, int timecount, Multibody_Set ** multibodysets, int*time_con)
{
    int multibodyii;
    
    sys = syst;
    n_times = timecount;
    multibodies.clear();
    multibodies.resize(n_times);
    time_conversion = new int[sys->show_n_timesteps()];
    
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=time_con[timeii];
    }
    
    for(int timeii=0;timeii<n_times;timeii++)
    {
      multibodies[timeii].resize(multibodysets[timeii]->show_n_multibodies());
      for(multibodyii=0;multibodyii<multibodysets[timeii]->show_n_multibodies();multibodyii++)
      {
	multibodies[timeii][multibodyii]=multibodysets[timeii]->show_multibody(multibodyii);
      }
      
    }
    check_n_bodies();
}


/*Bounds are inclusive*/
Multibody_List::Multibody_List(const Multibody_List& copy, int bound, bool greater)
{
  int timeii, mbodyii;
  
  Multibody*multibody;
  
  sys=copy.sys;
  
  int system_times=sys->show_n_timesteps();
  n_times=copy.n_times;
  
  
  multibodies.resize(n_times);
  
  time_conversion = new int [system_times];
  
  for(timeii=0;timeii<system_times;timeii++)
  {
    time_conversion[timeii]=copy.time_conversion[timeii];
  }
  
  if(greater)
  {
    for(timeii=0;timeii<n_times;timeii++)
    {
      for(mbodyii=0;mbodyii<copy.multibodies[timeii].size();mbodyii++)
      {
	multibody=copy.multibodies[timeii][mbodyii];
	if(multibody->show_n_bodies()>=bound)
	{
	  multibodies[timeii].push_back(multibody);
	}
      }
    }
  }
  else
  {
    for(timeii=0;timeii<n_times;timeii++)
    {
      for(mbodyii=0;mbodyii<copy.multibodies[timeii].size();mbodyii++)
      {
	multibody=copy.multibodies[timeii][mbodyii];
	if(multibody->show_n_bodies()<=bound)
	{
	  multibodies[timeii].push_back(multibody);
	}
      }
    }
  }
}



/*Bounds are inclusive*/
Multibody_List::Multibody_List(const Multibody_List& copy, int lowerbound, int upperbound)
{
  int timeii, mbodyii;
  int system_times=sys->show_n_timesteps();
  Multibody*multibody;
  
  sys=copy.sys;
  n_times=copy.n_times;
  
  multibodies.resize(n_times);
  
  time_conversion = new int [system_times];
  
  for(timeii=0;timeii<system_times;timeii++)
  {
    time_conversion[timeii]=copy.time_conversion[timeii];
  }
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    for(mbodyii=0;mbodyii<copy.multibodies[timeii].size();mbodyii++)
    {
      multibody=copy.multibodies[timeii][mbodyii];
      if(multibody->show_n_bodies()<=upperbound&&multibody->show_n_bodies()>=lowerbound)
      {
	multibodies[timeii].push_back(multibody);
      }
    }
  }
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



/*Method that resets object with system and a pointer to a list of Multibody_Set onjects - fully creates list with multibodies in it*/
void Multibody_List::set(System * syst, vector<Multibody_Set*> multibodysets, int*time_con)
{
    int multibodyii;
    int timecount = multibodysets.size();
    
    delete [] time_conversion;
    
    sys = syst;
    n_times = timecount;
    multibodies.clear();
    multibodies.resize(n_times);
    time_conversion = new int[sys->show_n_timesteps()];
    
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=time_con[timeii];
    }
    
    for(int timeii=0;timeii<n_times;timeii++)
    {
      multibodies[timeii].resize(multibodysets[timeii]->show_n_multibodies());
      for(multibodyii=0;multibodyii<multibodysets[timeii]->show_n_multibodies();multibodyii++)
      {
	multibodies[timeii][multibodyii]=multibodysets[timeii]->show_multibody(multibodyii);
      }
      
    }
    check_n_bodies();
}


Multibody_List Multibody_List::operator+(const Multibody_List & comparison) const
{
    int internal_time1, internal_time2;		//variables to hold internal Trajectory_List time indexes
    int new_internal_time = 0;			//variable to hold internal Trajectory_List time index for new Trajectory_List
    int last_internal_time1 = -1;			//variable to hold previous value of internal Trajectory_List time index
    int last_internal_time2 = -1;			//variable to hold previous value of internal Trajectory_List time index
    int n_internal_times;				//number of internal times for new Trajectory_List object
    int multibodyii;


    /*check if system objects are the same for the two Multibody_Lists; if not, return error*/
    if(sys!=comparison.sys)
    {
        cout << "Error: list systems do not match.\n";
        exit(1);
    }

    /*loop to count the number of independent internal times the new Multibody_List must have*/
    for(int system_timeii=0;system_timeii<sys->show_n_timesteps();system_timeii++)
    {
        internal_time1 = convert_time(system_timeii);
        internal_time2 = comparison.convert_time(system_timeii);

        if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
        {
            new_internal_time++;
        }

        last_internal_time1=internal_time1;
        last_internal_time2=internal_time2;
    }
    n_internal_times = new_internal_time;


    new_internal_time = -1;
    last_internal_time1 = -1;
    last_internal_time2 = -1;
    
    Multibody_List new_list(sys, n_internal_times);
    

    /*calculate new time_conversion table and new inclusion list based on intersection of two trajectories*/
    for(int system_timeii=0;system_timeii<sys->show_n_timesteps();system_timeii++)
    {
        internal_time1 = convert_time(system_timeii);
        internal_time2 = comparison.convert_time(system_timeii);

        if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
        {
            new_internal_time++;
	    new_list.multibodies[new_internal_time]=multibodies[internal_time1];
	    for(multibodyii=0;multibodyii<comparison.multibodies[internal_time2].size();multibodyii++)
	    {
	      (new_list.multibodies[new_internal_time]).push_back(comparison.multibodies[internal_time2][multibodyii]);
	    }
        }

        new_list.time_conversion[system_timeii]=new_internal_time;//set value in time conversion table to current new_internal_time for this system_timeii

        last_internal_time1=internal_time1;
        last_internal_time2=internal_time2;
    }

    return new_list;
}

/*Determines whether all multibodies in the list have the same number of bodies.*/
void Multibody_List::check_n_bodies()
{
    int timeii, multibodyii;

    // If nothing has been pushed yet, report "not determined"
    if (multibodies.empty() || multibodies[0].empty())
      {
        n_bodies = -2;
        return;
      }

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

	internal_time=convert_time(currenttime);
	for(int trajectoryii=0;trajectoryii<multibodies[internal_time].size();trajectoryii++) 	
	{
		analysis->listkernel(multibodies[internal_time][trajectoryii], timegap, currenttime, nexttime);
	}
}


/*Returns maximum bodies per multibody.*/
int Multibody_List::maxsize()const
{
    int timeii, multibodyii, max_bodies;

    max_bodies = 0;


    for(timeii=0;timeii<multibodies.size();timeii++)
    {
        for(multibodyii=0;multibodyii<multibodies[timeii].size();multibodyii++)
        {
            if(max_bodies<multibodies[timeii][multibodyii]->show_n_bodies())
            {
                max_bodies=multibodies[timeii][multibodyii]->show_n_bodies();
            }
        }
    }
    
    return max_bodies;
}
