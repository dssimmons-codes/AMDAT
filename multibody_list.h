/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_List class - stores list of multibody objects*/
/*Written by David S. Simmons*/

#include "multibody.h"

namespace std{
  
  class System;
 
  class Multibody_List
  {
    System sys;
    
    vector<vector<Multibody*>> multibodies
    
    //Multibody *** multibodies;
    vector<int> time_conversion;
    
    int n_bodies;	//number of bodies in each multibody if all multibodies have same number of bodies; -1 otherwise; -2 if not checked
    
    int n_times;
    
    int convert_time(int timeii)const{return time_conversion[timeii];};	//convert requested time (Where the index is the time index from the system object) to internal time index
    
    public:
      Multibody_List();
      Multibody_List(const Multibody_List &);
      Multibody_List operator = (const Multibody_list &)
      
      void listloop(Multibody_Analysis* analysis, int time);
      
      int show_n_bodies(){return n_bodies};
  }
  