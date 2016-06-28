/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Provisional_Multibodies class methods - provides parent class for any analysis method that needs to create multibodies*/
/*Written by David S. Simmons*/



#include "provisional_multibodies.h"
#include "system.h"
#include "control.h"
#include "multibody_list.h"
#include "trajectory_list.h"
#include <stdlib.h>
#include <string>


using namespace std;

Provisional_Multibodies::Provisional_Multibodies()
{
  time_conversion = new int [1];
  n_times=1;
  
}


Provisional_Multibodies::Provisional_Multibodies(const Provisional_Multibodies& copy)
{
  multibodies=copy.multibodies;
  set_pointers=copy.set_pointers;
  basename=copy.basename;
  set_names=copy.set_names;
  n_times=copy.n_times;
  time_conversion=new int [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=copy.time_conversion[timeii];
  }
  
}


Provisional_Multibodies Provisional_Multibodies::operator=(const Provisional_Multibodies& copy)
{
  if(this!=&copy)
  {
    multibodies=copy.multibodies;
    set_pointers=copy.set_pointers;
    basename=copy.basename;
    set_names=copy.set_names;
    n_times=copy.n_times;
    delete [] time_conversion;
    time_conversion=new int [n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }   
  }
  
  return *this;
}


Provisional_Multibodies::~Provisional_Multibodies()
{
  delete [] time_conversion;
}



void Provisional_Multibodies::convert(System* syst, Control* control, string setname, string traj_typename, bool centertype)
{
  create_multibody_sets();
  add_sets_to_system(syst, setname, traj_typename, centertype);
  add_lists_to_control(syst, control);
}



void Provisional_Multibodies::create_multibody_sets()
{
  Multibody_Set * mbset;
  int timeii;
  
  set_pointers.resize(multibodies.size());
  
  for(timeii=0;timeii<multibodies.size();timeii++)
  {
    mbset = new Multibody_Set;
    mbset->set(multibodies[timeii]);   
    set_pointers[timeii]=mbset;
  }
  
}



void Provisional_Multibodies::add_sets_to_system(System* syst, string setname, string traj_typename, bool centertype)
{
  int setii;
  basename = setname;
  
  for(setii=0;setii<set_pointers.size();setii++)
  {
    set_names.push_back(basename+to_string(static_cast<long long>(setii)));
    syst->add_multibody_set(set_names[setii],set_pointers[setii]);
    syst->create_trajectory_set(set_names[setii], set_names[setii], traj_typename, centertype);
    
  }
}



void Provisional_Multibodies::add_lists_to_control(System* syst, Control* control)
{
  vector<Trajectory_Set*> trajset;
  
  Multibody_List* new_multibody_list;
  Trajectory_List*new_trajectory_list;
  
  new_multibody_list = new Multibody_List;
  new_trajectory_list = new Trajectory_List;	//how is this going to work? Need way to build trajectory list from multiple multibody lists
  
  cout<<"\n"<<set_pointers.size()<<endl;cout.flush();
  
  new_multibody_list->set(syst,set_pointers,time_conversion);
  control->add_multibody_list(new_multibody_list,basename);
  
  for(int setii=0;setii<set_pointers.size();setii++)
  {
    trajset.push_back((Trajectory_Set*)(set_pointers[setii]));
  }
  
  new_trajectory_list->set(syst,trajset,time_conversion);
  control->add_trajectorylist(new_trajectory_list, basename); 
}

