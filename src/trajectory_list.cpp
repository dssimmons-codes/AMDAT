/*Methods for trajectories class - stores a list of trajectories*/
/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Written by David S. Simmons*/

#include "trajectory_list.h"
#include "analysis.h"
#include <fstream>
#include "system.h"
#include <stdlib.h>
#include <omp.h>

using namespace std;

Trajectory_List::Trajectory_List(int timecount, int cap)
{

	capacity=cap;
	n_times = timecount;
    n_atomtypes = 0;
	//Allocate memory for arrays
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	time_conversion = new int [0];
	included = new Boolean_List [n_times];
	sys = 0;
	for(int timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
		for(int trajii=0;trajii<capacity;trajii++)
		{
			trajectories[timeii][trajii]=0;
		}
	}


}

Trajectory_List::Trajectory_List(System* syst, int timecount, int cap, Boolean_List * boollist, int*time_conv)
{
  sys=syst;
  capacity=cap;

    n_atomtypes = sys->show_n_atomtypes();
    n_times = timecount;

      //Allocate memory for arrays
    trajectories = new Trajectory ** [n_times];
    n_trajectories = new int [n_times];
    included = new Boolean_List [n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
        trajectories[timeii] = new Trajectory * [capacity];
        n_trajectories[timeii]=0;
        included[timeii] = boollist[timeii];
        for(int trajii=0; trajii<capacity; trajii++)
        {
            trajectories[timeii][trajii]=0;
        }
    }

    time_conversion = new int [sys->show_n_timesteps()];

    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=time_conv[timeii];
    }

    trajlist_from_boollist();

}

Trajectory_List Trajectory_List::operator = (const Trajectory_List & copy)
{
  if(this!=&copy)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {
      delete [] trajectories[timeii];
    }
    delete [] trajectories;
    delete [] n_trajectories;
    delete [] included;
    delete [] time_conversion;

    sys = copy.sys;
    n_times = copy.n_times;
    capacity = copy.capacity;
    n_atomtypes = copy.n_atomtypes;
    cout<<"\n"<<capacity;
    
    if(capacity==0&&sys!=0){capacity=sys->show_n_trajectories();}

    /*END NEW LINES*/

    trajectories = new Trajectory**[n_times];
    n_trajectories = new int[n_times];
    included = new Boolean_List[n_times];

    for(int timeii=0; timeii<n_times; timeii++)
    {
      trajectories[timeii]=new Trajectory* [capacity];
      n_trajectories[timeii]=copy.n_trajectories[timeii];
      cout<<"\n"<<n_trajectories[timeii];
      included[timeii]=copy.included[timeii];
      for(int trajii=0; trajii<n_trajectories[timeii]; trajii++)
      {
        trajectories[timeii][trajii]=copy.trajectories[timeii][trajii];
      }
      for(int trajii=n_trajectories[timeii]; trajii<capacity; trajii++)
      {
        trajectories[timeii][trajii]=0;
      }
    }

    time_conversion = new int[sys->show_n_timesteps()];
    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
  }
  return *this;
}

Trajectory_List::Trajectory_List(const Trajectory_List & copy)
{
  sys=copy.sys;
  n_times=copy.n_times;
  capacity=copy.capacity;
  n_atomtypes = copy.n_atomtypes;
  trajectories = new Trajectory**[n_times];
  n_trajectories = new int[n_times];
  included = new Boolean_List[n_times];

  for(int timeii=0; timeii<n_times; timeii++)
  {
    trajectories[timeii]=new Trajectory* [capacity];
    n_trajectories[timeii]=copy.n_trajectories[timeii];
    included[timeii]=copy.included[timeii];
    for(int trajii=0; trajii<capacity; trajii++)
    {
      trajectories[timeii][trajii]=copy.trajectories[timeii][trajii];
    }
  }

  time_conversion = new int[sys->show_n_timesteps()];
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=copy.time_conversion[timeii];
  }
}

Trajectory_List::~Trajectory_List()
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] n_trajectories;
  delete [] time_conversion;
  delete [] included;
}


Trajectory* Trajectory_List::operator()(int trajii)
{
	Trajectory * nullp;
	nullp=0;

	if(trajii>=n_trajectories[0]) return nullp;
	else return trajectories[0][trajii];
}


Trajectory* Trajectory_List::operator()(int timeii,int trajii)
{
	Trajectory * nullp;
	nullp=0;

	if(convert_time(timeii)>=n_times) return nullp;
	else if(trajii>=n_trajectories[convert_time(timeii)]) return nullp;
	else return trajectories[convert_time(timeii)][trajii];
}


void Trajectory_List::listloop(Analysis* analysis, int time)
{
	int current_time;
	current_time=convert_time(time);
	#pragma omp parallel for schedule(dynamic) if(analysis->isThreadSafe()) // This is needed to make Structure Factor parallel 
	for(int trajectoryii=0;trajectoryii<n_trajectories[current_time];trajectoryii++)
	{
		analysis->listkernel(trajectories[current_time][trajectoryii]);
	}
}

void Trajectory_List::listloop(Analysis* analysis, int timegap, int curTime, int nextTime)
{
	int current_time;

	current_time=convert_time(curTime);
	for(int trajectoryii=0;trajectoryii<n_trajectories[current_time];trajectoryii++) 	
	{
		analysis->listkernel(trajectories[current_time][trajectoryii], timegap, curTime, nextTime);
	}
}

void Trajectory_List::listloop2(Analysis* analysis, Trajectory* traj1, int timegap, int curTime, int nextTime)
{
	int current_time;

	current_time=convert_time(curTime);
	for(int trajectoryii=0;trajectoryii<n_trajectories[current_time];trajectoryii++) 	
	{
		analysis->listkernel2(traj1, trajectories[current_time][trajectoryii], timegap, curTime, nextTime);
	}
}


int Trajectory_List::n_system_trajectories()const
{
  return sys->show_n_trajectories();
}


void Trajectory_List::write_count(string filename)const
{
	double avg_trajectories=0;


	/*calculate average number of atoms*/
	for(int timeii=0;timeii<n_times;timeii++)
	{
		avg_trajectories += double(n_trajectories[timeii])/double(n_times);
	}

	ofstream output(filename.c_str());
	output << "Average_trajectories:\t"<< avg_trajectories<<"\n";
	output << "Count_List\n";
	for(int timeii=0;timeii<n_times;timeii++)
	{
		output << "\t" << n_trajectories[timeii] << "\n";
	}
}


Trajectory_List Trajectory_List::time_intersection(int t1, int t2)const
{
  int time1,time2;
  Boolean_List * temp_included;
  int * temp_time_conversion;

  temp_time_conversion = new int [sys->show_n_trajectories()];
  Trajectory_List new_list;
  temp_included = new Boolean_List;

  for(int timeii=0;timeii<sys->show_n_trajectories();timeii++)
  {
    temp_time_conversion[timeii]=0;
  }

  time1 = convert_time(t1);
  time2 = convert_time(t2);

  temp_included[0] = included[time1]&&included[time2];

  new_list.set(sys,1,capacity,temp_included,temp_time_conversion);

  delete temp_included;
  delete [] temp_time_conversion;

  return new_list;
}




Trajectory_List Trajectory_List::time_union(int t1, int t2)const
{
  int time1,time2;
  Boolean_List * temp_included;
  int * temp_time_conversion;
  temp_time_conversion = new int [sys->show_n_trajectories()];
  Trajectory_List new_list;
  temp_included = new Boolean_List;

  for(int timeii=0;timeii<sys->show_n_trajectories();timeii++)
  {
    temp_time_conversion[timeii]=0;
  }

  Trajectory_List newlist;

  time1 = convert_time(t1);
  time2 = convert_time(t2);

  temp_included[0] = included[time1]||included[time2];

  new_list.set(sys,1,capacity,temp_included,temp_time_conversion);

  delete temp_included;
  delete [] temp_time_conversion;

  return new_list;
}



/*Method to build array of trajectories from previously provided Boolean_List array*/
void Trajectory_List::trajlist_from_boollist()
{
  int * trajid_list;
  int trajcount;

  for(int timeii=0;timeii<n_times;timeii++)
  {
    trajcount = included[timeii].show_n_included();
    trajid_list = new int [trajcount];
    n_trajectories[timeii]=included[timeii].show_trajectory_ids(trajcount,trajid_list);

    for(int trajii=0;trajii<n_trajectories[timeii];trajii++)
    {
      trajectories[timeii][trajii]=sys->show_trajectory(trajid_list[trajii]);
    }
    delete [] trajid_list;
  }

}

/*Method to set trajectory_list based on a Boolean_List array*/

void Trajectory_List::set(System* syst, int timecount, int cap, Boolean_List * boollist, int*time_conv)
{

  sys=syst;
  capacity=cap;

    n_atomtypes = sys->show_n_atomtypes();




    delete [] time_conversion;

    //deallocate previous allocations
    for(int timeii=0; timeii<n_times; timeii++)
    {
        delete [] trajectories[timeii];
    }
    delete [] trajectories;
    delete [] n_trajectories;
    delete [] included;

    n_times = timecount;

      //Allocate memory for arrays
    trajectories = new Trajectory ** [n_times];
    n_trajectories = new int [n_times];
    included = new Boolean_List [n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
        trajectories[timeii] = new Trajectory * [capacity];
        n_trajectories[timeii]=0;
        included[timeii] = boollist[timeii];
        for(int trajii=0; trajii<capacity; trajii++)
        {
            trajectories[timeii][trajii]=0;
        }
    }

    time_conversion = new int [sys->show_n_timesteps()];

    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=time_conv[timeii];
    }

    trajlist_from_boollist();

}


//Method to initialize trajectory list based on vector of trajectory sets
void Trajectory_List::set(System * syst, vector<Trajectory_Set*> trajectory_sets, int*time_conv)
{
  int trajii, timeii;
  
  sys=syst;
  //system=const_cast<System*>(sys);
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
	delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] n_trajectories;
  delete [] time_conversion;
  delete [] included;
  
  n_atomtypes = sys->show_n_atomtypes();
  n_times=trajectory_sets.size();
  
  capacity=sys->show_n_trajectories();
  
  
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  time_conversion=new int [sys->show_n_timesteps()];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_trajectories [timeii] = trajectory_sets[timeii]->show_n_trajectories();
    included[timeii].set(sys);
    trajectories[timeii] = new Trajectory * [n_trajectories[timeii]];
    for(trajii=0;trajii<n_trajectories[timeii];trajii++)
    {
      trajectories[timeii][trajii]=trajectory_sets[timeii]->show_trajectory(trajii);
      (included[timeii])(trajectory_sets[timeii]->show_trajectory(trajii)->show_trajectory_ID(),1);
    }
  }
  
  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=time_conv[timeii];
  }
}



void Trajectory_List::flatten_multibodies(const Multibody_List& mblist)
{
  int timeii, mbodyii, bodyii;
  int system_times;
  int n_mbodies, n_bodies;
  int max_trajectories;
  
  vector<Trajectory*> bodies;
  
  sys=mblist.sys;
  
  delete [] time_conversion;

  //deallocate previous allocations
  for(int timeii=0; timeii<n_times; timeii++)
  {
      delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] n_trajectories;
  delete [] included;
  
  n_times = mblist.n_times;
  
  system_times = sys->show_n_timesteps();
  n_atomtypes=sys->show_n_atomtypes();
  
  //copy time conversion table
  time_conversion = new int [system_times]; 
  for(timeii=0;timeii<system_times;timeii++)
  {
    time_conversion[timeii]=mblist.time_conversion[timeii];
  }
  
  //allocate memory for trajectory arrays
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  max_trajectories=0;
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(sys);
    n_trajectories[timeii]=0;
    n_mbodies=mblist.multibodies[timeii].size();
    for(mbodyii=0;mbodyii<n_mbodies;mbodyii++)
    {
      n_bodies=mblist.multibodies[timeii][mbodyii]->show_n_bodies();
      n_trajectories[timeii]+=n_bodies;
      bodies=mblist.multibodies[timeii][mbodyii]->show_bodies();
      for(bodyii=0;bodyii<n_bodies;bodyii++)
      {
	included[timeii](bodies[bodyii]->show_trajectory_ID(),1);
      }
    }
    if(n_trajectories[timeii]>max_trajectories)
    {
      max_trajectories=n_trajectories[timeii];
    }
    
  }
  
  capacity=max_trajectories;
  for(int timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [n_trajectories[timeii]];
  }
  
  trajlist_from_boollist();
}



void Trajectory_List::addtrajectory(int time, Trajectory* traj)
{
    trajectories[time][n_trajectories[time]]=traj;
    included[time](trajectories[time][n_trajectories[time]]->show_trajectory_ID(),bool(1));
    n_trajectories[time]++;
}

/*Method to return intersection of two Trajectory_Lists*/

Trajectory_List Trajectory_List::operator&& (const Trajectory_List & comparison) const
{

    int internal_time1, internal_time2;		//variables to hold internal Trajectory_List time indexes
    int new_internal_time = 0;			//variable to hold internal Trajectory_List time index for new Trajectory_List
    int last_internal_time1 = -1;			//variable to hold previous value of internal Trajectory_List time index
    int last_internal_time2 = -1;			//variable to hold previous value of internal Trajectory_List time index
    int n_internal_times;				//number of internal times for new Trajectory_List object
    int * temp_time_conversion;			//table of system frame index to internal time index conversions for new Trajectory_List
    int newcapacity;

    Boolean_List * temp_included;			//trajectory inclusion table for new Trajectory_List
    Trajectory_List new_list;

    temp_time_conversion = new int [sys->show_n_timesteps()];		//allocate memory for new time conversion table

    /*select new capacity as larger of the capacity of the two Trajectory_Lists*/
    if(capacity>=comparison.capacity)
    {
        newcapacity=capacity;
    }
    else
    {
        newcapacity=comparison.capacity;
    }


    /*check if system objects are the same for the two Trajectory_Lists; if not, return error*/
    if(sys!=comparison.sys)
    {
        cout << "Error: list systems do not match.\n";
        exit(1);
    }

    /*loop to count the number of independent internal times the new Trajectory_List must have*/
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

    temp_included = new Boolean_List [n_internal_times];		//allocate memory for new inclusion array

    new_internal_time = -1;
    last_internal_time1 = -1;
    last_internal_time2 = -1;

    /*calculate new time_conversion table and new inclusion list based on intersection of two trajectories*/
    for(int system_timeii=0;system_timeii<sys->show_n_timesteps();system_timeii++)
    {
        internal_time1 = convert_time(system_timeii);
        internal_time2 = comparison.convert_time(system_timeii);

        if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
        {
            new_internal_time++;
            temp_included[new_internal_time] = included[internal_time1]&&comparison.included[internal_time2];
        }

        temp_time_conversion[system_timeii]=new_internal_time;	//set value in time conversion table to current new_internal_time for this system_timeii

        last_internal_time1=internal_time1;
        last_internal_time2=internal_time2;
    }

    /*create new trajectory list and return*/
    new_list.set(sys,n_internal_times,newcapacity,temp_included,temp_time_conversion);

    delete [] temp_included;									//MODIFIED DSS
    delete [] temp_time_conversion;

    return new_list;
}



/*Method to return intersection of two Trajectory_Lists that have the same time conversion table*/

Trajectory_List Trajectory_List::same_timetable_intersection(const Trajectory_List & comparison) const
{
  Boolean_List * temp_included;			//trajectory inclusion table for new Trajectory_List
  Trajectory_List new_list;

  /*check if system objects are the same for the two Trajectory_Lists; if not, return error*/
  if(sys!=comparison.sys)
  {
    cout << "Error: list systems do not match.\n";
    exit(1);
  }

  temp_included = new Boolean_List [n_times];		//allocate memory for new inclusion array


  for(int timeii=0;timeii<n_times;timeii++)
  {
    temp_included[timeii]=included[timeii]&&comparison.included[timeii];
  }


  /*create new trajectory list and return*/
   new_list.set(sys,n_times,capacity,temp_included,time_conversion);

   delete [] temp_included;									//MODIFIED DSS

   return new_list;

}


/*Method to return union of two Trajectory_Lists*/


Trajectory_List Trajectory_List::operator|| (const Trajectory_List & comparison) const
{
  int internal_time1, internal_time2;		//variables to hold internal Trajectory_List time indexes
  int new_internal_time = 0;			//variable to hold internal Trajectory_List time index for new Trajectory_List
  int last_internal_time1 = -1;			//variable to hold previous value of internal Trajectory_List time index
  int last_internal_time2 = -1;			//variable to hold previous value of internal Trajectory_List time index
  int n_internal_times;				//number of internal times for new Trajectory_List object
  int * temp_time_conversion;			//table of system frame index to internal time index conversions for new Trajectory_List
  int newcapacity;

  Boolean_List * temp_included;			//trajectory inclusion table for new Trajectory_List
  Trajectory_List new_list;

  temp_time_conversion = new int [sys->show_n_timesteps()];		//allocate memory for new time conversion table

  /*select new capacity as larger of the capacity of the two Trajectory_Lists*/
  if(capacity>=comparison.capacity)
  {
    newcapacity=capacity;
  }
  else
  {
    newcapacity=comparison.capacity;
  }


  /*check if system objects are the same for the two Trajectory_Lists; if not, return error*/
  if(sys!=comparison.sys)
  {
    cout << "Error: list systems do not match.\n";
    exit(1);
  }

  /*loop to count the number of independent internal times the new Trajectory_List must have*/
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

  temp_included = new Boolean_List [n_internal_times];		//allocate memory for new inclusion array

  new_internal_time = -1;
  last_internal_time1 = -1;
  last_internal_time2 = -1;

  /*calculate new time_conversion table and new inclusion list based on intersection of two trajectories*/
  for(int system_timeii=0;system_timeii<sys->show_n_timesteps();system_timeii++)
  {
    internal_time1 = convert_time(system_timeii);
    internal_time2 = comparison.convert_time(system_timeii);
    if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
    {
      new_internal_time++;
      temp_included[new_internal_time] = included[internal_time1]||comparison.included[internal_time2];
    }

    temp_time_conversion[system_timeii]=new_internal_time;	//set value in time conversion table to current new_internal_time for this system_timeii

    last_internal_time1=internal_time1;
    last_internal_time2=internal_time2;
  }

  /*create new trajectory list and return*/
   new_list.set(sys,n_internal_times,newcapacity,temp_included,temp_time_conversion);

   delete [] temp_included;
   delete [] temp_time_conversion;

   return new_list;

}



void Trajectory_List::write_xyz(string filename)const
{
  ofstream output(filename.c_str());
  int internaltime;
  Coordinate coord;

  for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
  {
    internaltime = convert_time(timeii);
    output << n_trajectories[internaltime]<<"\n";
    output << "Atoms\n";

    for(int trajii=0;trajii<n_trajectories[internaltime];trajii++)
    {
      coord = trajectories[internaltime][trajii]->show_coordinate(timeii);
      output << trajectories[internaltime][trajii]-> show_type() << "\t" << coord.show_x() << "\t" << coord.show_y() << "\t" << coord.show_z() << "\n";
    }
  }
}

void Trajectory_List::write_full_xyz(string filename)const
{
    ofstream output(filename.c_str());
    int internaltime;
    Coordinate coord;

    for(int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {
        internaltime = convert_time(timeii);
        output << sys->show_n_atoms()<<"\n";
        output << "Atoms\n";
        int inc_trajii=0;

        for(int trajii=0;trajii< sys->show_n_atoms();trajii++)
        {
            if (included[internaltime](trajii))
            {
                coord = trajectories[internaltime][inc_trajii]->show_coordinate(timeii);
                output << trajectories[internaltime][inc_trajii]-> show_type() << "\t" << coord.show_x() << "\t" << coord.show_y() << "\t" << coord.show_z() << "\n";
                inc_trajii++;
            }
            else
            {
                output << sys->show_trajectory(trajii)->show_type() << "\t" << 1000 << "\t" << 1000 << "\t" << 1000 << "\n";
            }
        }
    }
}

bool Trajectory_List::is_included(int timeii, int trajii)
{
    int current_time;
    current_time = convert_time(timeii);
    return included[current_time](trajii);
}

void Trajectory_List::inversion(Trajectory_List* t_list,Trajectory_List* original)
{
    Boolean_List* new_included;
    int current_time;
    int original_time;
    new_included = new Boolean_List [n_times];

    for (int timeii=0;timeii<sys->show_n_timesteps();timeii++)
    {   current_time = convert_time(timeii);
        original_time = original->convert_time(timeii);
        new_included[current_time].set(sys);
        for (int trajii=0;trajii<sys->show_n_trajectories();trajii++)
        {
//            cout<< "time = "<<current_time<<endl;cout.flush();
//            cout << "traj = "<< trajii<<  endl;cout.flush();
            if(!(included[current_time](trajii)) && original->included[original_time](trajii))
            {
                new_included[current_time](trajii,1);
            }
            else
            {
                new_included[current_time](trajii,0);
            }
        }
    }
   t_list->set(sys, n_times, n_system_trajectories(), new_included, time_conversion);

}


