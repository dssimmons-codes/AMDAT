/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Multibody_List class - stores list of multibody objects*/
/*Written by David S. Simmons*/


namespace std{
  
  class System;
 
  class Multibody_List
  {
    System const* sys;
    Multibody *** multibodies;
    int * time_conversion;
    int * n_multibodies;
    
    int n_bodies;	//number of bodies in each multibody - must be the same for all multibodies
    
    int n_times;
    
    int convert_time(int timeii)const{return time_conversion[timeii];};	//convert requested time (Where the index is the time index from the system object) to internal time index
    
    public:
      Multibody_List();
      Multibody_List(const Multibody_List &);
      Multibody_List operator = (const Multibody_list &)
      
      void listloop(Multibody_Analysis* analysis, int time);
      
      int show_n_bodies(){return n_bodies};
  }
  