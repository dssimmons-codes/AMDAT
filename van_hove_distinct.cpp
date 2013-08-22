/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Van_Hove_Distinct: a class for distinct part of Van Hove correlation function.*/
/*Written by David S. Simmons*/

#include "van_hove_distinct.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;


/*----------------------------------------------------------------------------------------------*/

/*Method to intialize object*/
void Van_Hove_Distinct::initialize(System*sys, int bin_count, const Spacial_Decomposition* cells, float min_size, float value_max, float value_min)
{
  int timeii, binii;

  system = sys;
  cellobject = cells;
  n_bins = bin_count;
  min_value = value_min;

  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;

  if(min_size==0) {min_cell_size=value_max-value_min;}
  else {min_cell_size = min_size;}

  n_atoms_looped = 0;


  max_cell_range = ceil(max_value/min_cell_size);

  if(max_cell_range > cellobject->show_overhang())
  {
    cout << "Error: insufficient number of image cells in spacially decomposed atom data!  Increase overhang!";
    exit(1);
  }

  bin_size = (max_value-min_value)/float(n_bins);

  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  //allocate memory for van hove distinct correlation function and weighting and initialize to zero
  correlation = new int * [n_times];
  n_atoms = new double [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii]=0;
    correlation[timeii]=new int [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }

  weighting = system->timegap_weighting();	//get timegap weighting from system

}


/*----------------------------------------------------------------------------------------------*/

/*Constructor*/

Van_Hove_Distinct::Van_Hove_Distinct()
{
  system = 0;
  cellobject = 0;
  n_bins = 0;
  min_value = 0;
  max_value = 0;
  max_cell_range = 0;
  n_times = 0;
  bin_size = 0;
  timetable = 0;
  correlation = new int * [1];
  weighting = new int [0];
  correlation[0]=new int [0];
}


/*----------------------------------------------------------------------------------------------*/


/*Constructor*/
Van_Hove_Distinct::Van_Hove_Distinct(System*sys, int bin_count, const Spacial_Decomposition* cells, float min_size, float value_max, float value_min)
{
  initialize(sys, bin_count, cells, min_size, value_max, value_min);
}



/*----------------------------------------------------------------------------------------------*/


/*atom kernel method*/
void Van_Hove_Distinct::atomkernel(Trajectory * traj)
{
  system->displacement_loop(this, traj);
  n_atoms_looped++;		//increment number of atoms for which correlation has been binned.
}



/*----------------------------------------------------------------------------------------------*/



/*displacement kernel method*/
void Van_Hove_Distinct::displacementkernel(int timegap,int thisii, int nextii,Trajectory * traj)
{
  int deltaxii, deltayii, deltazii;
  Coordinate coordinate;
  int atomID;

  coordinate = traj->show_coordinate(thisii);
  atomID = traj->show_trajectory_ID();
  for(deltaxii=-max_cell_range;deltaxii<=max_cell_range;deltaxii++)
  {
    for(deltayii=-max_cell_range;deltayii<=max_cell_range;deltayii++)
    {
      for(deltazii=-max_cell_range;deltazii<=max_cell_range;deltazii++)
      {
        cellobject->loop_cell(this, timegap, nextii, coordinate, atomID, deltaxii, deltayii, deltazii);
      }
    }
  }
}


/*----------------------------------------------------------------------------------------------*/


/*cellkernel method*/
void Van_Hove_Distinct::cellkernel(int timegap, Coordinate coordinate1, Coordinate coordinate2, int atom_type)
{
  float distance = (coordinate2-coordinate1).length();

  bin(timegap,distance);
}



/*----------------------------------------------------------------------------------------------*/


void Van_Hove_Distinct::postprocess()
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
		n_atoms[timeii] = n_atoms_looped;
	}

}



/*----------------------------------------------------------------------------------------------*/


void Van_Hove_Distinct::list_displacementkernel(int timegapii, int thisii, int nextii, Particle_List* particle_list)
{
  nexttime = nextii;
  timegap = timegapii;

  particle_list->atomloop(this,thisii);
}




/*----------------------------------------------------------------------------------------------*/


void Van_Hove_Distinct::atomlist_kernel(int thisii, int species_index, int molecule_index, int atom_type, int atom_index,Particle_List* particle_list)
{

  int deltaxii, deltayii, deltazii;
  Atom_Trajectory * atom;
  Coordinate coordinate;
  int atomID;

  atom = ((system->show_molecule(species_index, molecule_index))->show_atom_trajectory(atom_type, atom_index));
  coordinate = atom->show_coordinate(thisii);
  atomID = atom->show_trajectory_ID();
  for(deltaxii=-max_cell_range;deltaxii<=max_cell_range;deltaxii++)
  {
    for(deltayii=-max_cell_range;deltayii<=max_cell_range;deltayii++)
    {
      for(deltazii=-max_cell_range;deltazii<=max_cell_range;deltazii++)
      {
        cellobject->loop_cell(this, timegap, nexttime, coordinate, atomID, deltaxii, deltayii, deltazii);
      }
    }
  }

  n_atoms[timegap]+=double(particle_list->show_n_atoms(thisii))/double(weighting[timegap]);
}

#ifdef NO

/*----------------------------------------------------------------------------------------------*/


void Van_Hove_Distinct::list_displacementkernel(int timegapii, int thisii, int nextii)
{
  nexttime = nextii;
  timegap = timegapii;

  trajectory_list->listloop(this,thisii);
}


/*----------------------------------------------------------------------------------------------*/


void Van_Hove_Distinct::listkernel(Trajectory* traj)
{

  int deltaxii, deltayii, deltazii;
  Atom_Trajectory * atom;
  Coordinate coordinate;
  int atomID;

  coordinate = traj->show_coordinate(thisii);
  atomID = traj->show_trajectory_ID();
  for(deltaxii=-max_cell_range;deltaxii<=max_cell_range;deltaxii++)
  {
    for(deltayii=-max_cell_range;deltayii<=max_cell_range;deltayii++)
    {
      for(deltazii=-max_cell_range;deltazii<=max_cell_range;deltazii++)
      {
        cellobject->loop_cell(this, timegap, nexttime, coordinate, atomID, deltaxii, deltayii, deltazii);
      }
    }
  }

  n_atoms[timegap]+=double(particle_list->show_n_atoms(thisii))/double(weighting[timegap]);
}

#endif
