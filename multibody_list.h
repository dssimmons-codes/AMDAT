/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_List class - stores list of multibody objects*/
/*Written by David S. Simmons*/

#ifndef MULTIBODY_LIST
#define MULTIBODY_LIST

#include "multibody.h"
#include "multibody_set.h"

namespace std{
  
  class Multibody_Analysis;

  class System;

  class Multibody_List
  {
    System * sys;

    vector<vector<Multibody*>> multibodies;

    //Multibody *** multibodies;
    int * time_conversion;

    int n_bodies;	//number of bodies in each multibody if all multibodies have same number of bodies; -1 otherwise; -2 if not checked
    void check_n_bodies();  //determine n_bodies;

    int n_times;

    int convert_time(int timeii)const{return time_conversion[timeii];};	//convert requested time (Where the index is the time index from the system object) to internal time index

    public:
      Multibody_List();
      Multibody_List(const Multibody_List &);
      Multibody_List operator = (const Multibody_List &);
      ~Multibody_List();

      Multibody_List(System * sys, int timecount);    //construct Multibody_List, setting the total number of times in the list.
      Multibody_List (System * sys, Multibody_Set * multibodyset);


      Multibody_List operator + (const Multibody_List &); //combine two Multibody_Lists
      void set(System * sys, Multibody_Set * multibodyset);

      int show_n_bodies(){return n_bodies;}; //return number of bodies in multibodies in list if they all contain the same number of bodies. Return -1 if they are not all the same or -2 if it has not been determined.
      int show_n_multibodies(int timeii){return multibodies[convert_time(timeii)].size();};   //return number of multibodies in list at a given time
      
      void listloop(Multibody_Analysis* analysis, int timegap, int currenttime, int nexttime);
      
  };
}

#endif