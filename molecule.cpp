/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for molecules class*/
/*Written by David S. Simmons*/

#include <stdio.h>
#include <stdlib.h>
#include "molecule.h"
#include <iostream>
#include <math.h>

using namespace std;



/*Overloaded constructor to set number of atom types*/
Molecule::Molecule(int types, int ts)
{
  trajtype=molecule;
  n_atomtypes = types;	
  n_atoms = new int[types];			//allocated memory for array of atom counts
  atoms = new Atom_Trajectory*[0];		//allocate zero memory for trajectory array
  unwrapped = new Coordinate [0];		//allocate zero memory for center of mass trajectory array
  coordinates = new Coordinate[0];		//allocate zero memory for wrapped center of mass trajectory array
  mass = 0;
  n_timesteps = ts;
  atoms_unwrapped=0;
  is_unwrapped=0;
  is_wrapped=0;
}



/*Overloaded constructor to set atom counts*/
Molecule::Molecule(int types, int * natoms, int ts=0)
{
  int typeii;				//index over types
  trajtype=molecule;
  n_atomtypes = types;	
  n_atoms = new int[types];		//instantiate array of atom numbers
  n_atoms = natoms;			//copy array into member
  total_atoms = 0;
  atoms=new Atom_Trajectory*[n_atomtypes];				//allocate memory for row of trajectory array columns
  
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    atoms[typeii] = new Atom_Trajectory[n_atoms[typeii]];		//allocate memory for each column of trajectory array
    total_atoms += n_atoms[typeii];
  }
  
  unwrapped = new Coordinate [0];	//allocate zero memory for center of mass trajectory array
  coordinates = new Coordinate[0];	//allocate zero memory for wrapped center of mass trajectory array
  mass = 0;
  
  n_timesteps = ts;
  atoms_unwrapped = 0;
  is_unwrapped=0;
  is_wrapped=0;
}




Molecule::~Molecule()
{
	clear_memory();
}



/*Clear memory allocated to arrays*/
void Molecule::clear_memory()
{
  int typeii;						//index over types
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    delete [] (*atoms);			//delete column of atoms of each type
  }
  delete [] atoms;				//delete row of columns	
  delete [] n_atoms;				//delete list of atom counts
}




/*Method to calculate mass of molecule by summing up mass of constituent atoms*/
void Molecule::calculate_mass()
{
  int typeii, atomii;
  mass = 0;
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    for(atomii=0;atomii<n_atoms[typeii];atomii++)
    {
      mass += atoms[typeii][atomii].show_mass();
    }
  }
}




void Molecule::calculate_center_of_mass(const Coordinate * box_size, Coordinate ** box_boundaries)
{
  int typeii, atomii, timeii;
  calculate_mass();		//recalculate mass to make sure it is up to date
  Coordinate temp_coordinate1;
  Coordinate temp_coordinate2;
  Coordinate displacement;
  int atomcount;
  Coordinate com_displacement;
  int goodframe = 0;
  bool good;
  
  coordinates = new Coordinate [n_timesteps];
  unwrapped = new Coordinate [n_timesteps];
  
  /*find a frame at which center of mass can be determined unambigously*/
  for(timeii=0;timeii<n_timesteps;timeii++)	//loop over timesteps
  {
    good = 1;
    atomcount=0;
    for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
    {
      for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
      {
	if (atomcount==0){temp_coordinate1=atoms[typeii][atomii].show_coordinate(timeii);}
	temp_coordinate2=atoms[typeii][atomii].show_coordinate(timeii);
	displacement = temp_coordinate2 - temp_coordinate1;
	
	if((displacement.show_x()<box_size[timeii].show_x()*.55&&displacement.show_x()>box_size[timeii].show_x()*.45)||(displacement.show_y()<box_size[timeii].show_y()*.55&&displacement.show_y()>box_size[timeii].show_y()*.45)||(displacement.show_z()<box_size[timeii].show_z()*.55&&displacement.show_z()>box_size[timeii].show_z()*.45))  //check for excessively large real displacements
	{
	  good = 0;
	}
	atomcount++;
      }
    }
    if(good)
    {
      goodframe = timeii;
      break;
    }
  }
  
  if(!good)
  {
    cout << "\nWarning: Distance between atoms within molecule " << moleculeID << " at all times too large to reliably determine center of mass due to possible molecular box-crossing. A larger box may be necessary to resolve this ambiguity.";
    goodframe=0;
  }
  
  
  
  /*calculate center of mass at first time at which unambiguous unwrapping is possible*/
  
  coordinates[goodframe].set(0,0,0);		//initialize center of mass
  atomcount = 0;
  for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
  {
    for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
    {
      if (atomcount==0){temp_coordinate1=atoms[typeii][atomii].show_coordinate(goodframe);}
      temp_coordinate2=atoms[typeii][atomii].show_coordinate(goodframe);
      displacement = temp_coordinate2 - temp_coordinate1;
      
      displacement-=box_size[goodframe]*((displacement/(box_size[goodframe]*.55)).integer());
      temp_coordinate2=temp_coordinate1+displacement;
      coordinates[goodframe]+=temp_coordinate2*atoms[typeii][atomii].show_mass()/mass;
      atomcount++;
    }
  }
  coordinates[goodframe]-=box_size[goodframe]*((coordinates[goodframe]-box_boundaries[goodframe][0])/box_size[goodframe]).coord_floor();
  
  
  /*determine coordinates at initial time on the basis of the previously determined good time*/
  /*this step is not really necessary (could calculate everything directly based on the good time), but causes the unwrapped coordinates equal the wrapped coordinates at the initial time.*/
  coordinates[0].set(0,0,0);		//initialize center of mass
  atomcount = 0;
  com_displacement.set(0,0,0);
  for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
  {
    for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
    {
      com_displacement += (atoms[typeii][atomii].show_unwrapped(0)-atoms[typeii][atomii].show_unwrapped(goodframe))*atoms[typeii][atomii].show_mass()/mass;
    }
  }
  coordinates[0] = coordinates[goodframe] + com_displacement;
  coordinates[0]-=box_size[goodframe]*((coordinates[goodframe]-box_boundaries[goodframe][0])/box_size[goodframe]).coord_floor();
  unwrapped[0] = coordinates[0];
  
  
  /*calculate coordinates at later times*/
  for(timeii=1;timeii<n_timesteps;timeii++)	//loop over timesteps
  {
    coordinates[timeii].set(0,0,0);		//initialize center of mass
    atomcount = 0;
    com_displacement.set(0,0,0);
    for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
    {
      for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
      {
	com_displacement += (atoms[typeii][atomii].show_unwrapped(timeii)-atoms[typeii][atomii].show_unwrapped(0))*atoms[typeii][atomii].show_mass()/mass;;
      }
    }
    unwrapped[timeii] = coordinates[0] + com_displacement;
    coordinates[timeii]= unwrapped[timeii] - box_size[timeii]*((unwrapped[timeii]-box_boundaries[timeii][0])/box_size[timeii]).coord_floor();
  }
 
}




//works only for molecule that does not cross box boundaries
void Molecule::calculate_center_of_mass()
{
  int typeii, atomii, timeii;
  calculate_mass();		//recalculate mass to make sure it is up to date
  
  if(!atoms_unwrapped){cout<<"\nError: Coordinates not unwrapped.";exit(1);}
  
  unwrapped = new Coordinate [n_timesteps];
  
  for(timeii=0;timeii<n_timesteps;timeii++)	//loop over timesteps
  {
    unwrapped[timeii].set(0,0,0);		//initialize center of mass
    for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
    {
      for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
      {
        unwrapped[timeii] += atoms[typeii][atomii].show_unwrapped(timeii)*atoms[typeii][atomii].show_mass()/mass;
      }
    }
  }
}




//works only for molecule that does not cross box boundaries
float Molecule::gyration_radius()
{
  int timeii, typeii, atomii;
 
  float * radius;
  radius = new float [n_timesteps];
  float time_average_radius = 0;
  
  calculate_center_of_mass();
  
  if(!atoms_unwrapped)
  {
     cout << "\nError: Molecule not unwrapped!\n";exit(1);
  }
  
  calculate_mass(); //recalculate mass to make sure it is up to date
  
  for(timeii=0;timeii<n_timesteps;timeii++)	//loop over timesteps
  {
    radius[timeii] = 0;			//initialize gyration radius
    for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
    {
      for(atomii=0;atomii<n_atoms[typeii];atomii++)	//loop over atoms
      {
        radius[timeii] += (atoms[typeii][atomii].show_unwrapped(timeii)-unwrapped[timeii])&(atoms[typeii][atomii].show_unwrapped(timeii)-unwrapped[timeii])*atoms[typeii][atomii].show_mass()/mass;
      }
    }
    time_average_radius += radius[timeii];
    radius[timeii] = pow(radius[timeii],0.5);
    
  }
  time_average_radius = pow(time_average_radius/n_timesteps,0.5);
  
  return time_average_radius;
}




/*Set object values, defining array sizes and appropriately assigning memory.  Also initially define atom objects by passing them number of timesteps.*/
void Molecule::set(int types, int * natoms, int timesteps)
{
  int typeii;						//index over types
  clear_memory();
  n_atomtypes = types;	
  n_atoms = new int[types];		//instantiate array of atom numbers
  n_atoms = natoms;			//copy array into member
  total_atoms = 0;
  int atomii;		  
  
  n_timesteps=timesteps;
  
  atoms=new Atom_Trajectory*[n_atomtypes];				//allocate memory for row of atom array columns
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    atoms[typeii] = new Atom_Trajectory[n_atoms[typeii]];		//allocate memory for each column of atom array
    total_atoms += n_atoms[typeii];
    for(atomii=0;atomii<n_atoms[typeii];atomii++)
    {
      (atoms[typeii][atomii]).set(typeii+1,timesteps);
    }
  }

  is_unwrapped=0;
  is_wrapped=0;
  
}



int Molecule::typecount(int type)				//method to return number of atoms of a given type
{
  if(type<n_atomtypes)
  {
    return n_atoms[type];
  }
  else
  {
    cout << "Invalid atom type!";
    exit(1);
  }
}



int Molecule::atomcount()			//method to return total number of atoms
{
return total_atoms;
}



Atom_Trajectory Molecule::show_atom(int type, int index)
{
  if(type>n_atomtypes)
  {
    cout << "Atom type index greater than number of atom types!";
    exit(1);
  }
  if(index>n_atoms[type])
  {
    cout << "Atom index greater than number of atoms of this type!";
    exit(1);
  }
  
  return atoms[type][index];
}



void Molecule::set_coordinate(int type, int index, const Coordinate & coordinate, int timestep)
{
  if(type>n_atomtypes)
  {
    cout << "Atom type index greater than number of atom types!";
    exit(1);
  }
  if(index>n_atoms[type])
  {
    cout << "Atom index greater than number of atoms of this type!";
    exit(1);
  }
  (atoms[type][index]).set(coordinate, timestep);
}



void Molecule::set_unwrapped(int type, int index, const Coordinate & coordinate, int timestep)
{
  if(type>n_atomtypes)
  {
    cout << "Atom type index greater than number of atom types!";
    exit(1);
  }
  if(index>n_atoms[type])
  {
    cout << "Atom index greater than number of atoms of this type!";
    exit(1);
  }
  (atoms[type][index]).set_unwrapped(coordinate, timestep);
}



void Molecule::set_velocity(int type, int index, const Coordinate & coordinate, int timestep)
{
  if(type>n_atomtypes)
  {
    cout << "Atom type index greater than number of atom types!";
    exit(1);
  }
  if(index>n_atoms[type])
  {
    cout << "Atom index greater than number of atoms of this type!";
    exit(1);
  }
  (atoms[type][index]).set_velocity(coordinate, timestep);
}




int Molecule::unwrap_atoms(const Coordinate & box_size)
{
  int typeii;
  int atomii;

  for(typeii=0; typeii<n_atomtypes; typeii++)
  {
    for(atomii=0; atomii<n_atoms[typeii]; atomii++)
    {
      (atoms[typeii][atomii]).unwrap(box_size);
    }
  }
  
  atoms_unwrapped = 1;
  
  return total_atoms;
}



int Molecule::wrap_atoms(const Coordinate * box_size, Coordinate ** box_boundaries)
{
  int typeii;
  int atomii;

  for(typeii=0; typeii<n_atomtypes; typeii++)
  {
    for(atomii=0; atomii<n_atoms[typeii]; atomii++)
    {
      (atoms[typeii][atomii]).wrap(box_size, box_boundaries);
    }
  }
  
  atoms_unwrapped = 1;
  
  return total_atoms;
}


Coordinate Molecule::show_unwrapped(int atom_type, int atom_index, int timestep)const
{
  return atoms[atom_type][atom_index].show_unwrapped(timestep);
}



Atom_Trajectory* Molecule::show_atom_trajectory(int atom_type, int atom_index)
{
  return &atoms[atom_type][atom_index];
}


float Molecule::atom_distance(int atom_type, int atom_index, int start_time, int finish_time)
{
  return atoms[atom_type][atom_index].distance(start_time,finish_time);
}



/*Pass unique ID of molecule down to constituent atoms*/
void Molecule::ID_to_atoms()
{
  int typeii, atomii;
  
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    for(atomii=0;atomii<n_atoms[typeii];atomii++)
    {
      (atoms[typeii][atomii]).set_moleculeID(moleculeID);
    }
  }
}