/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Header for class Boolean_List, which stores an array of booleans - one for each trajectory in the system - denoting inclusion in a trajectory list*/
/*Written by David S. Simmons*/

#include <vector>

#ifndef BOOLEAN_LIST
#define BOOLEAN_LIST

namespace std
{

class System;

class Boolean_List
{
  System * system;
  //bool * included;		//array holding a boolean for each trajectory in the system
  vector <bool> included;	//array holding a boolean for each trajectory in the system

  public:
  Boolean_List(System * sys);
  Boolean_List(System * sys, int * inc, int n_included);
  Boolean_List();
  Boolean_List(const Boolean_List &); //MEM - copy constructor
  ~Boolean_List();
  Boolean_List operator= (const Boolean_List &);


  void set(System * sys);
  void set(System * sys, int * inc, int n_included);

  bool operator() (int)const;
  void operator() (int,bool);
  Boolean_List operator&& (const Boolean_List &);
  Boolean_List operator|| (const Boolean_List &);
  bool operator== (const Boolean_List &);
  bool operator!= (const Boolean_List &);

  int show_n_included()const;
  int show_trajectory_ids(int, int *)const;
  int first_included()const;

  void grow_list(){included.push_back(0);};		//add additional '0' to end of boolean vector - use whenever a new trajectory is added to system
  void grow_list(int growsize); //add additional '0's to end of boolean vector - use whenever multiple new trajectories are added to system

};
}
#endif
