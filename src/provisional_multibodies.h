/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Provisional_Multibodies class - provides parent class for any analysis method that needs to create multibodies*/
/*Written by David S. Simmons*/


#ifndef PROVISIONAL_MULTIBODIES
#define PROVISIONAL_MULTIBODIES

#include "multibody.h"
#include "multibody_set.h"
#include "multibody_list.h"
#include <vector>
#include <sstream>
#include <stdlib.h>

namespace std
{

class Control;
class System;
  
class Provisional_Multibodies
{
  protected:
    vector<vector<Multibody>> multibodies;
    vector<Multibody_Set*> set_pointers;
    string basename;
    vector<string> set_names;
    int*time_conversion;
    int n_times;
   
    void create_multibody_sets();
    void add_sets_to_system(System* syst, string setname, string traj_typename, bool centertype);	//need to code functionality to add sets to system data structures
    void add_lists_to_control(System * syst, Control* control);
  public:
    
    Provisional_Multibodies();
    Provisional_Multibodies(const Provisional_Multibodies&);
    Provisional_Multibodies operator=(const Provisional_Multibodies&);
    ~Provisional_Multibodies();
    
    Multibody_List* temporary_multibodies(System* syst);
    void delete_sets();
    
    void convert(System* syst, Control* control, string setname, string traj_typename, bool centertype);
  
};


}

#endif