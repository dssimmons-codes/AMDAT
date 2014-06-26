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
  int typeii;
  total_atoms=0;
  n_atomtypes = types;	
  n_atoms = new int[n_atomtypes];		//allocate memory for array of atom counts
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    n_atoms[typeii]=0;
  }
  atoms = new Atom_Trajectory*[0];		//allocate zero memory for trajectory array
  n_timesteps = ts;
  atoms_unwrapped=0;
  moleculeID=-1;
}



/*Overloaded constructor to set atom counts*/
Molecule::Molecule(int types, int * natoms, int ts=0)
{
  int typeii;				//index over types
  n_atomtypes = types;	
  n_atoms = new int[types];		//instantiate array of atom numbers
  
  total_atoms = 0;
  atoms=new Atom_Trajectory*[n_atomtypes];				//allocate memory for row of trajectory array columns
  
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    n_atoms[typeii] = natoms[typeii];			//copy array into member
    atoms[typeii] = new Atom_Trajectory[n_atoms[typeii]];		//allocate memory for each column of trajectory array
    total_atoms += n_atoms[typeii];
  }
  
  n_timesteps = ts;
  atoms_unwrapped = 0;
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
      (atoms[typeii][atomii]).set(typeii+1,timesteps,1);
    }
  }
  
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


/*return multibody consisting of all atoms in molecule*/
Multibody Molecule::create_multibody()const
{
  int typeii,atomii;
  int bodyii=0;
  
  Multibody multibody(total_atoms);
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    for(atomii=0;atomii<n_atoms[typeii];atomii++)
    {
      multibody.set(bodyii,&(atoms[typeii][atomii]));
      bodyii++;
    }
  }
  
  return multibody;
}



/*return multibody consisting of all atoms of type typeii in molecule*/
Multibody Molecule::create_multibody(int typeii)const
{
  int atomii;
  int bodyii=0;
  
  Multibody multibody(n_atoms[typeii]);
  for(atomii=0;atomii<n_atoms[typeii];atomii++)
  {
    multibody.set(bodyii,&(atoms[typeii][atomii]));
    bodyii++;
  }
  
  return multibody;
}





/*return multibody consisting of atoms specified by a list of types and indices*/
Multibody Molecule::create_multibody(int n_bodies, int * typeii, int * index)const
{
  int atomii;
  int bodyii=0;
  
  Multibody multibody(n_bodies);
  for(atomii=0;atomii<n_bodies;atomii++)
  {
    multibody.set(bodyii,&(atoms[typeii[atomii]][index[atomii]]));
    bodyii++;
  }
  
  return multibody;
}
