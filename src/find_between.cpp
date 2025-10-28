// Copyright (C) 2010-2025 David S. Simmons

//Add standard methods up here. Contructors, equality operator, etc.
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

#include "find_between.h"

using namespace std;

/*Default constructor*/
Find_Between::Find_Between()
{
  int timeii;

  n_times = 0;

  capacity = 0;
  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  time_conversion = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
  }

}

Find_Between::Find_Between(System * syst)
{
  int timeii;

  sys=syst;
  system = const_cast<System*>(sys);

  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.
  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  time_conversion = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
    time_conversion[timeii]=timeii;
  }
}

Find_Between::Find_Between(System * syst, float d_cutoff, float t_cutoff, bool only_diff)
{
  int timeii;

  sys=syst;
  system = const_cast<System*>(sys);

  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.
  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  time_conversion = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
    time_conversion[timeii]=timeii;
  }
  only_diff_molecule=only_diff;
  dist_cutoff=d_cutoff;
  dist_cutoffx2=d_cutoff*2;
  costheta_cutoff=t_cutoff;
}


Find_Between::Find_Between(const Find_Between & copy):Trajectory_List(copy), Analysis_Onetime(copy)
{
//  Trajectory_List::Trajectory_List(copy);
//  Analysis::Analysis(copy);
  only_diff_molecule=copy.only_diff_molecule;
  dist_cutoff=copy.dist_cutoff;
  dist_cutoffx2=copy.dist_cutoff*2;
  costheta_cutoff=copy.costheta_cutoff;

}

Find_Between::~Find_Between()
{
////  Trajectory_List::~Trajectory_List();
////  Analysis::~Analysis();
}

void Find_Between::set(System * syst, float d_cutoff, float t_cutoff, bool only_diff)
{
  int timeii;
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] trajectories[timeii];
  }
  delete [] trajectories;
  delete [] included;
  delete [] n_trajectories;
  delete [] time_conversion;
  
  sys=syst;
  system = const_cast<System*>(sys);

  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.
  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  trajectories = new Trajectory ** [n_times];
  n_trajectories = new int [n_times];
  included = new Boolean_List [n_times];
  time_conversion = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectories[timeii] = new Trajectory * [capacity];
    n_trajectories[timeii]=0;
    included[timeii].set(sys);
    time_conversion[timeii]=timeii;
  }
  only_diff_molecule=only_diff;
  dist_cutoff=d_cutoff;
  dist_cutoffx2=d_cutoff*2;
  costheta_cutoff=t_cutoff;

//    for(int timeii=0;timeii<n_times;timeii++)
//    {
//        time_conversion[timeii]=time_conv[timeii];
//    }

}

Find_Between Find_Between::operator = (const Find_Between & copy)
{
  sys=copy.system;
  system = const_cast<System*>(sys);
  trajectory_list = copy.trajectory_list;
  n_trajectories = copy.n_trajectories;
  included = copy.included;
//  Trajectory_List::operator=(copy);
//  Analysis::operator=(copy);

  capacity=system->show_n_atoms()+system->show_n_molecules();   //This sets how much memory will be allocated for trajectory list at each time. This is just an estimate.

  n_times = system->show_n_timesteps();     //sets number of times for memory allocation to equal number of times in system

  only_diff_molecule=copy.only_diff_molecule;
  dist_cutoff=copy.dist_cutoff;
  dist_cutoffx2=copy.dist_cutoff*2;
  costheta_cutoff=copy.costheta_cutoff;

}

void Find_Between::timekernel2(int timeii)
{
   trajectory_list->listloop(this,0, timeii, 0);
}


void Find_Between::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Find_Between::listkernel2(Trajectory* traj1, Trajectory* traj2, int timegapii, int thisii, int nextii) //need to check back on which of these time indices is the right one.
{
  if (is_included(thisii, traj1->show_trajectory_ID())) //skip if traj1 is already included
    return;

  //check if distance between traj1 and traj2 is under dist_cutoff before proceeding over tertiary loop
  Coordinate dist_vector_1 = (traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).vector_unwrapped(system->size(thisii)); //calculate shortest distance between two coordinates, taking into account periodic boundaries
  float dist_1 = dist_vector_1.length();
  if (dist_1 > dist_cutoffx2)
    return;

  int n_trajs=trajectory_list2->show_n_trajectories(thisii);  //record how many particles are in the second trajectory list at this time

  Trajectory* traj3; //holds second reference bead
  Atom_Trajectory* atom_traj2 = dynamic_cast<Atom_Trajectory*>(traj2); //to compare atomIDs with second reference bead

  for(int trajii=0;trajii<n_trajs;trajii++)  //here is the second loop over particles
  {

    traj3 = (*trajectory_list2)(thisii,trajii);  //this is the syntax to pull the second particle in set 2

    //skip if traj3 is the same bead as traj2
    Atom_Trajectory* atom_traj3 = dynamic_cast<Atom_Trajectory*>(traj3);
    if (atom_traj2->show_atomID() == atom_traj3->show_atomID())
      continue;

    //if only_diff_molecule==1, check that traj2 and traj3 are part of a different cluster
    if (only_diff_molecule && atom_traj2->show_moleculeID() == atom_traj3->show_moleculeID())
      continue;
    // TODO: we may want to also code in a version that only computes for intra-molecular/intra-cluster

    //check if distance between traj1 and traj3 is under dist_cutoff
    Coordinate dist_vector_2 = (traj3->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).vector_unwrapped(system->size(thisii));
    float dist_2 = dist_vector_2.length();
    if (dist_2 > dist_cutoffx2)
      continue;

    // check dist_cutoff criteria
    float distx2 = dist_1+dist_2;
    if (distx2 > dist_cutoffx2)
      continue;

    // check costheta_cutoff criteria
    float cos_theta = dist_vector_1&dist_vector_2/dist_1/dist_2;
    if (cos_theta > costheta_cutoff)
      continue;

    addtrajectory(thisii,traj1); //this line will add the trajectory to the trajectory list
    break; //break to not repeat adding traj1
//    //print stuff to checkk
//    cout << "Bead: "    << traj1->show_coordinate(thisii).show_x() << " " << traj1->show_coordinate(thisii).show_y() << " " <<traj1->show_coordinate(thisii).show_z() << "\n";
//    cout << "Center1: " << atom_traj2->show_moleculeID() << " " << atom_traj2->show_atomID() << " " << traj2->show_coordinate(thisii).show_x() << " " << traj2->show_coordinate(thisii).show_y() << " " <<traj2->show_coordinate(thisii).show_z() << "\n";
//    cout << "Center2: " << atom_traj3->show_moleculeID() << " " << atom_traj3->show_atomID() << " " << traj3->show_coordinate(thisii).show_x() << " " << traj3->show_coordinate(thisii).show_y() << " " <<traj3->show_coordinate(thisii).show_z() << "\n";
//    cout << "Diff1: "   << dist_vector_1.show_x() << " " << dist_vector_1.show_y() << " " <<dist_vector_1.show_z() << "\n";
//    cout << "Diff2: "   << dist_vector_2.show_x() << " " << dist_vector_2.show_y() << " " <<dist_vector_2.show_z() << "\n";
//    cout << dist_1 << " " << dist_2 << " " << cos_theta << "\n";
  }
}

void Find_Between::preprocess()
{
}

void Find_Between::preprocess2()
{
}

void Find_Between::bin_hook(Trajectory_List *,int,int,int)
{
}

void Find_Between::postprocess_bins()
{
  //I don't think any postprocessing is needed in this case
}
