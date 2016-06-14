/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Header for class Boolean_List, which stores an array of booleans - one for each trajectory in the system - denoting inclusion in a trajectory list*/
/*Written by David S. Simmons*/

#include "boolean_list.h"
#include "system.h"
#include <stdlib.h>

using namespace std;


Boolean_List::Boolean_List()
{
  system = 0;
}

Boolean_List::Boolean_List(System * sys)
{
  int trajii;

  system = sys;

  included.resize(system->show_n_trajectories()); 
}


Boolean_List::Boolean_List(System * sys, int * inc, int n_included)
{
   /**Create boolean list from list of trajectory IDs of included trajectories
    * @author David Simmons
    * @date 3/1/2013
    */
  int trajii;
  
  system = sys;

  //included = new bool [system->show_n_trajectories()];
  included.resize(system->show_n_trajectories());

  //for(trajii=0;trajii<system->show_n_trajectories();trajii++)
  //{
  //  included[trajii]=0;
  //}

  for(trajii=0;trajii<n_included;trajii++)
  {
    included[inc[trajii]]=1;
  }

}



/*New MEM - copy constructor*/
Boolean_List::Boolean_List(const Boolean_List & copy)
{
  system = copy.system;
  //included = new bool [system->show_n_trajectories()];
  included.resize(system->show_n_trajectories());
  for(int trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    included[trajii]=copy.included[trajii];
  }

}
/*End New MEM*/



/*NEW SECTION DSS*/
Boolean_List::~Boolean_List()
{
}
/*END NEW SECTION*/

int Boolean_List::show_size()const
{
  return included.size();
}


Boolean_List Boolean_List::operator= (const Boolean_List & copy)
{

   if(this!=&copy)
   {
     int max_trajectories;
     system = copy.system;
     //delete [] included;
     //included = new bool [system->show_n_trajectories()];
     update_size();
     max_trajectories = system->show_n_trajectories();
     if(copy.show_size()<max_trajectories)
     {
       max_trajectories=copy.show_size();
     }
     for(int trajii=0;trajii<max_trajectories;trajii++)
     {
       included[trajii]=copy.included[trajii];
     }
   }
  return *this;				//NEW LINE DSS
}


void Boolean_List::set(System * sys)
{
  //delete [] included;

  system = sys;

  //included = new bool [system->show_n_trajectories()];
  included.resize(system->show_n_trajectories());
  //for(int trajii=0;trajii<system->show_n_trajectories();trajii++)
  //{
 //   included[trajii]=0;
  //}
}


void Boolean_List::set(System * sys, int * inc, int n_included)
{
  int trajii=0;
//  delete [] included;

  
  system = sys;

  //included = new bool [system->show_n_trajectories()];
  included.resize(system->show_n_trajectories());
//  for(trajii=0;trajii<system->show_n_trajectories();trajii++)
//  {
//    included[trajii]=0;
//  }

  for(trajii=0;trajii<n_included;trajii++)
  {
    included[inc[trajii]]=1;
  }

}



bool Boolean_List::operator()(int index)
{
  update_size();
  if(index>=system->show_n_trajectories())
  {
    cout << "Error: index greater than number of trajectories in system.\n";
    exit(1);
  }
  else
  {
    return included [index];
  }
}


void Boolean_List::operator()(int index, bool inc)
{
  update_size();
  if(index>=system->show_n_trajectories())
  {
    cout << "Error: index greater than number of trajectories in system.\n";
    exit(1);
  }
  else
  {
    included [index] = inc;
  }

}



Boolean_List Boolean_List::operator&&(Boolean_List& comparison)
{
  Boolean_List intersection(system);
  int trajii;

  if(system!=comparison.system)
  {
    cout << "Error: cannot combine boolean lists referencing different system objects.\n";
    exit(1);
  }
  update_size();
  comparison.update_size();

  for(trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    intersection.included[trajii] = included[trajii] && comparison.included[trajii];
  }

  return intersection;
}



Boolean_List Boolean_List::operator||(Boolean_List& comparison)
{
  Boolean_List intersection(system);
  int trajii;

  if(system!=comparison.system)
  {
    cout << "Error: cannot combine boolean lists referencing different system objects.\n";
    exit(1);
  }

  update_size();
  comparison.update_size();
  
  for(trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    intersection.included[trajii] = included[trajii] || comparison.included[trajii];
  }

  return intersection;
}





bool Boolean_List::operator==(Boolean_List& comparison)
{
  bool same;
  if(system!=comparison.system)
  {
    cout << "Error: cannot combine boolean lists referencing different system objects.\n";
    exit(1);
  }
  update_size();
  comparison.update_size();
  for(int trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    if ((included[trajii] != comparison.included[trajii]))
        {
            same = 0;
            break;
        }
        else
        {
            same = 1;
        }
  }

  return same;
}

bool Boolean_List::operator!=(Boolean_List& comparison)
{
  return !(*this == comparison);
}




int Boolean_List::show_n_included()const
{
  int trajii;
  int n_included=0;

  for(trajii=0;trajii<included.size();trajii++)
  {
    n_included+=int(included[trajii]);
  }

  return n_included;
}



int Boolean_List::show_trajectory_ids(int alloc_size, int * idlist)
{
  int trajii, includedii;
  int n_included;
  n_included = show_n_included();
  includedii=0;


  update_size();
  for(trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    if(included[trajii])
    {
      idlist[includedii]=trajii;
      includedii++;
      if(includedii>=alloc_size){return n_included;}
    }
  }

  return n_included;
}


int Boolean_List::first_included()const
{
  for(int trajii=0;trajii<system->show_n_trajectories();trajii++)
  {
    if(included[trajii]){return trajii;}
  }

  return int(-1);
}


void Boolean_List::update_size()
{
  if(included.size()<system->show_n_trajectories())
  {
    included.resize(system->show_n_trajectories());
  }
}