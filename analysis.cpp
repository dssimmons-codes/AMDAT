/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for parent Analysis class*/
/*Written by David S. Simmons*/

#include "system.h"
#include "particle_list.h"
#include <iostream>
#include <stdlib.h>

#define ARGMAX 100

using namespace std;

#include "tokenize.h"


/*----------------Constructors and assignment----------------------*/

/*Default constructor*/
Analysis::Analysis()
{
  system = 0;
  trajectory_list=0;
}


/*Copy constructor*/
Analysis::Analysis(const Analysis & copy)
{
  system = copy.system;
  trajectory_list = copy.trajectory_list;
}
  

/*Assignment operator*/
Analysis Analysis::operator =(const Analysis & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    trajectory_list = copy.trajectory_list;
  }
  
  return *this;
}


/*Methods for employing the loops implemented in class System over various sets of atoms*/

/*Do analysis on time-varying particle list, with different list of particles corresponding to first time element of each exponential block*/
void Analysis::atomlist(Particle_List * particle_list)
{preprocess_list(particle_list);particle_list->atomloop(this);postprocess_list(particle_list);}

void Analysis::all_moleculecom()
{system->loop_all_moleculecom(this);}

void Analysis::species_moleculecom(int species_index)
{preprocess();system->loop_species_moleculecom(this, species_index);postprocess();}

void Analysis::species_moleculecom(string species_name)
{
  int species_index;
  species_index = system->show_species_index(species_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  system->loop_species_moleculecom(this, species_index);
}

void Analysis::atom_species(int species_index, int atom_type, int atom_index)
{preprocess();system->loop_atom_species(this, species_index, atom_type, atom_index);postprocess();}

void Analysis::atom_species(string species_name, string atomtype_name, int atom_index)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_atom_species(this, species_index, atom_type, atom_index);
}

void Analysis::type_molecule(int species_index, int molecule_index, int atom_type)
{preprocess();system->loop_type_molecule(this, species_index, molecule_index, atom_type);postprocess();}

void Analysis::type_molecule(string species_name, int molecule_index, string atomtype_name)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_molecule(this, species_index, molecule_index, atom_type);
}
		
void Analysis::molecule(int species_index, int molecule_index)
{preprocess();system->loop_molecule(this, species_index, molecule_index);postprocess();}

void Analysis::molecule(string species_name, int molecule_index)
{
  int species_index;
  species_index = system->show_species_index(species_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  system->loop_molecule(this, species_index, molecule_index);
}

void Analysis::species(int species_index)
{preprocess();system->loop_species(this,species_index);postprocess();}


void Analysis::species(string species_name)
{
  int species_index;
  species_index = system->show_species_index(species_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  system->loop_species(this,species_index);
}


void Analysis::type_species(int species_index, int atom_type)
{preprocess();system->loop_type_species(this,species_index,atom_type);postprocess();}

void Analysis::type_species(string species_name, string atomtype_name)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_species(this,species_index,atom_type);
}


void Analysis::type_system(int atom_type)
{preprocess();system->loop_type_system(this,atom_type);postprocess();}


void Analysis::type_system(string atomtype_name)
{
  int atom_type;
  atom_type = system->show_atomtype_index(atomtype_name);
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_system(this,atom_type);
}


void Analysis::all()
{system->loop_system(this);}


void Analysis::analyze(string runline)
{
  int n_args;
  int args_used=0;
  int args_needed;
  string args [ARGMAX];
  string command;
  n_args = tokenize(runline,args);
  
  if(n_args==0){cout << "Error: No atom set command found.";exit(1);}
  
  preprocess();
  
  while(args_used<n_args)
  {
    command = args[args_used];
    args_used++;
    if(command == "atom_species")
    {
      args_needed = atom_species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      atom_species(args[args_used],args[args_used+1],atoi(args[args_used+2].c_str()));
      args_used+=args_needed;
    }
    else if(command == "type_molecule")
    {
      args_needed = type_molecule();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      
      type_molecule(args[args_used],atoi(args[args_used+1].c_str()),args[args_used+2]);
      args_used+=args_needed;
    }
    else if(command == "molecule")
    {
      args_needed = molecule();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      molecule(args[args_used],atoi(args[args_used+1].c_str()));
      args_used+=args_needed;
    }
    else if(command == "species")
    {
      args_needed = species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      species(args[args_used]);
      args_used+=args_needed;
    }
    else if(command == "type_species")
    {
      args_needed = type_species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      type_species(args[args_used],args[args_used+1]);
      args_used+=args_needed;
    }
    else if(command == "type_system")
    {
      args_needed = type_system();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      type_system(args[args_used]);
      args_used+=args_needed;
    }
    else if(command == "all")
    {
      args_needed = 0;
    
      all();
    }
    else if(command == "moleculecom_all")
    {
      args_needed = 0;
      
      all_moleculecom();
    }
    else if (command == "moleculecom_species")
    {
      args_needed = species_moleculecom();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      species_moleculecom(args[args_used]);
      args_used+=args_needed;
    }
    else
    {
      cout << "\nAtom set command "<<command<< " for analysis not understood.\n";
      exit(1);
    }
  }
  
  postprocess();
}
