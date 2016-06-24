/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Provisional_Multibodies class methods - provides parent class for any analysis method that needs to create multibodies*/
/*Written by David S. Simmons*/

using namespace std;

#include "provisional_multibodies.h"
#include "system.h"
#include "control.h"
#include "multibody_list.h"
#include "trajectory_list.h"







Provisional_Multibodies::Provisional_Multibodies()
{
  time_conversion = new[1];
  
}


Provisional_Multibodies::Provisional_Multibodies(const Provisional_Multibodies&)
{
  
}


Provisional_Multibodies::Provisional_Multibodies operator=(const Provisional_Multibodies&)
{
  
}



void Provisional_Multibodies::create_multibody_sets()
{
  Multibody_Set * mbset;
  int timeii;
  
  for(timeii=0;timeii<multibodies.size();timeii++)
  {
    mbset = new Multibody_Set;
    mbset.set(multibodies[timeii]);   
    set_pointers[timeii]=mbset;
  }
  
}



void Provisional_Multibodies::add_sets_to_system(System* syst, string setname, string traj_typename, bool centertype)
{
  int setii;
  basename = setname;
  
  for(setii=0;setii<set_pointers.size();setii++)
  {
    set_names[setii]= = basename+string(itoa(setii));
    syst->add_multibody_set(set_names[setii],set_pointers[setii]);
    syst->create_trajectory_set(set_names[setii], set_names[setii], traj_typename, centertype);
    
  }
}



void Provisional_Multibodies::add_lists_to_control(System* syst, Control* control)
{
  Multibody_List* new_multibody_list;
  Trajectory_List*new_trajectory_list;
  
  new_multibody_list = new Multibody_List;
  new_trajectory_list = new Trajectory_List;	//how is this going to work? Need way to build trajectory list from multiple multibody lists
  
  new_multibody_list->set(syst,set_pointers,time_conversion);
  control->add_multibody_list(new_multibody_list,basename);
  
  new_trajectory_list->set(syst,set_pointers,time_conversion);
  control->add_trajectorylist(new_trajectory_list, basename); 
}

void Provisional_Multibodies::convert(System* syst, Control* control, string setname, string traj_typename, bool centertype)
{
  create_multibody_sets();
  add_sets_to_sytem( syst, setname, traj_typename, centertype);
  add_lists_to_control(syst, control);
}