/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*System class methods*/
/*Written by David S. Simmons*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "system.h"
#include <time.h>
#include <sstream>
#include <string>
#include "version.h"
#include "control.h"
#include <vector>
#include <float.h>
#include <omp.h>

#define CPLUSPLUS
#ifndef TACC
#include "xdrfile_xtc.h"
#endif

using namespace std;
#include "tokenize.h"
#include "progress.h"
#define ARGMAX 100


/*-------------------------------------------------------------------------------------*/
/*-----------------------------------CONSTRUCTORS -------------------------------------*/
/*-------------------------------------------------------------------------------------*/

System::System()
{
}
/*----------------------------------------------------------------------------------*/


/*Constructor that reads system information from input file*/
System::System(vector<string> file_in, bool ensemble)
{
  /** Create system object
  * @param file_in - user input file
  * @param ensemble - 0 if constant volume (default), 1 if non-constant volume
  * @author David S. Simmons
  **/

  string line, fileline;
  string args [ARGMAX];
  int n_args, argii, speciesii, moleculeii, typeii, atomii, atomcount, moleculecount;
  string trajectory_type;
  int extra_trajectories = 0;	//extra space to allocate for additional trajectories

  /*Initialize booleans specifying state of coordinates*/
  wrapped = 0;
  unwrapped = 0;
  boxified = 0;

  np = ensemble;	//determine whether system is constant volume varies (1 if yes)

  /*read type of trajectory file to be read in*/
//  getline(*file_in,line);
  line = Control::read_line();
//  line=Control::replace_constants(line);
  n_args = tokenize(line, args);
  trajectory_type=args[0];
  if(n_args==2)
  {
    extra_trajectories=atoi(args[1].c_str());
    cout << "\nAllocating space for " << extra_trajectories << " additional trajectories.";
  }

  /*read trajectory filename information*/
//  getline(*file_in,fileline);
  fileline = Control::read_line();
  /*read timestep information*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++){args[argii]="";}

//  getline(*file_in,line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  timetype = args[0];

  if(timetype=="linear")
  {
    n_timesteps = atoi(args[1].c_str());
    n_exponentials = n_timesteps-1;
    exp_base = 1.0;
    n_exponential_steps = 1;
    first_exponent=0;
    frt=0;
    time_unit = atof(args[2].c_str());
  }
  else if (timetype == "exponential")
  {
    n_exponentials = atoi(args[1].c_str());
    n_exponential_steps = atoi(args[2].c_str());
    exp_base = atof(args[3].c_str());
    frt = atoi(args[4].c_str());
    first_exponent = atoi(args[5].c_str());
    n_timesteps = n_exponential_steps*n_exponentials+1;
    time_unit = atof(args[6].c_str());
  }
  else if (timetype == "snapshot")
  {
    n_timesteps = 1;
    n_exponentials = 1;
    exp_base = 1.0;
    n_exponential_steps = 0;
    first_exponent = 0;
    frt = 0;
    time_unit = 1.0;
  }
  else
  {
    cout << "Error: time scheme " << timetype <<" not recognized!\n";
    exit(1);
  }

  /*allocate memory for box-size info*/
  box_size = new Coordinate[n_timesteps];
  box_boundary = new Coordinate * [n_timesteps];
  for (int timeii=0;timeii<n_timesteps;timeii++) box_boundary[timeii]=new Coordinate [2];
  rho = new float [n_timesteps];

  /*create table of times in trajectory*/
  create_timelist();

  read_trajectory(trajectory_type,file_in,fileline);
  cout << "\nTrajectory data read successfully.";

  /*initialize master lists of trajectory objects (foreign key arrays)*/
  moleculecount=0;
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
	  for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
	  {
		  moleculecount++;		//count up total number of molecules
	  }
  }
  total_molecules=moleculecount;	//set total number of molecules in system
  moleculecount=0;

  total_trajectories=total_molecules+n_atoms+extra_trajectories;
  atomlist = new Atom_Trajectory*[n_atoms];
  moleculelist = new Molecule*[total_molecules];
  trajectorylist = new Trajectory*[total_trajectories];


  atomcount=0;
  /*generate a unique ID for each atom in the system (serves as a foreign key for atoms)*/
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
    {
      for(typeii=0;typeii<n_atomtypes;typeii++)
      {
        for(atomii=0;atomii<molecules[speciesii][moleculeii].typecount(typeii);atomii++)
        {
          atomlist[atomcount]=molecules[speciesii][moleculeii].show_atom_trajectory(typeii,atomii);	//add atom to master atom list
	  trajectorylist[atomcount]=atomlist[atomcount];	//add atom to master trajectory list
          atomlist[atomcount]->set_atomID(atomcount);		//set atom id in atom_trajectory object to correspond to position in master atom list
	  atomlist[atomcount]->set_trajectory_ID(atomcount);	//set trajectory id in atom_trajectory object to correspond to position in master trajectory list
          atomcount++;
        }
      }
    }
  }

  /*generate a unique ID for each molecule in the system (serves as foreigh key for molecules)*/
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
	  for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
	  {
		moleculelist[moleculecount]=&molecules[speciesii][moleculeii];	//add molecule to master molecule list
		moleculelist[moleculecount]->set_type(n_atomtypes+speciesii);
		trajectorylist[moleculecount+atomcount]=moleculelist[moleculecount];	//add molecule to master trajectory list
		moleculelist[moleculecount]->set_moleculeID(moleculecount);		//set atom id in atom_trajectory object to correspond to position in master atom list
		moleculelist[moleculecount]->set_trajectory_ID(moleculecount+atomcount);	//set trajectory id in atom_trajectory object to correspond to position in master trajectory list
		moleculelist[moleculecount]->set_species(speciesii);		//set internal register of species of molecule
		moleculelist[moleculecount]->ID_to_atoms();	//Pass unique ID of molecule down to constituent atoms
		moleculecount++;

	  }
  }

  /*Calculate molecules centers of mass*/
  cout << "\nGenerating molecule center of mass trajectories.";
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
    {
      molecules[speciesii][moleculeii].calculate_center_of_mass(box_size,box_boundary);
    }
  }


  displacement_limit=0;

  /*determine number of timegaps*/
  n_timegaps = n_exponential_steps + n_exponentials -1 +int(!frt);
  cout << "\nSystem created.";

  /*determine system density*/
  for(int timeii=0; timeii< n_timesteps; timeii++)
  {
    rho[timeii] = n_atoms/(box_size[timeii].show_x()*box_size[timeii].show_y()*box_size[timeii].show_z());
  }
}




/*-------------------------------------------------------------------------------------*/
/*---------------------METHODS TO READ IN TRAJECTORY DATA -----------------------------*/
/*-------------------------------------------------------------------------------------*/


/*Method to determine what type of trajectory file to read*/
void System::read_trajectory(string trajectory_type, vector<string> file_in, string fileline)
{
  if(np)
  {
    if(trajectory_type=="xyz"||trajectory_type=="xyz_log"||trajectory_type=="xtc")
    {
      cout << "Error: AMDAT can currently only handle custom LAMMP trajectory files for non-constant volume (system_np) trajectories.\n";
      exit(0);
    }
    else if (trajectory_type=="custom")
    {
      custom_prep(file_in,fileline);
    }
    else
    {
      cout<<"Error: trajectory file type "<<trajectory_type<<" not recognized.";
      exit(1);
    }
  }
  else
  {
    if(trajectory_type=="xyz")
    {
      xyz_prep(file_in,fileline);
    }
    else if(trajectory_type=="xyz_log")
    {
      xyz_prep_withlog(file_in,fileline);
    }
    #ifndef TACC
    else if(trajectory_type=="xtc")
    {
      xtc_prep(file_in,fileline);
    }
    #endif
    else if (trajectory_type=="custom")
    {
      custom_prep(file_in,fileline);
    }
    else
    {
      cout<<"Error: trajectory file type "<<trajectory_type<<" not recognized.";
      exit(1);
    }
  }
}



/*------------------------Methods to read xyz trajectories --------------------- ------*/


/*----------------------------------------------------------------------------------*/
void System::xyz_prep(vector<string> file_in, string fileline)
{
  string line;
  string args [ARGMAX];
  int n_args;
  int ** natoms;						//array of number of atoms of each type in each species
  int speciesii;						//species type index
  int typeii, argii;
  int atomii;							//atom type index
  float L_x, L_y, L_z;						//temporary storage for box dimensions
  float xlo, xhi, ylo, yhi, zlo, zhi;
  string xyzfilename;
  bool trajtemplate_given=0;			//does user specify a trajectory template?
  string trajtemplate;

  n_args = tokenize(fileline, args);
  xyzfilename = args[0];
  if(n_args>1)
  {
    trajtemplate_given=1;
    trajtemplate = args[1];
  }

  ifstream ifile(xyzfilename.c_str());
  if (!ifile)
  {
	Control::throw_error("XYZ Trajectory file not found.", true);
  }

  /*read in species names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  n_species = n_args/2;
  atoms_per_species = new int [n_species];
  species_name = new string [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    species_name[speciesii] = args[speciesii*2];
    atoms_per_species[speciesii]=0;
  }
  n_molecules = new int[n_species];	//allocate memory for molecule-number array
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    n_molecules[speciesii] = atoi(args[speciesii*2+1].c_str());
  }

  natoms = new int*[n_species];			//allocate memory for atom-number-per-molecule array

  /*read in atom type names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_atomtypes = tokenize(line, args);
  atomtype_name = new string [n_atomtypes];
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    atomtype_name[typeii] = args[typeii];
  }

  /*read in atomtype counts for each species*/
  for(speciesii=0; speciesii<n_species; speciesii++)
  {
    natoms[speciesii] = new int [n_atomtypes];
    line = "";
    for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//    getline(*file_in,line);
    line = Control::read_line();
    n_args = tokenize(line, args);
    if(n_args!=n_atomtypes){cout << "Error: number of atom counts ("<<n_args<<") does not match number of atom types("<<n_atomtypes<<").\n";exit(1);}
    for(atomii=0; atomii<n_atomtypes; atomii++)
    {
      natoms[speciesii][atomii] = atoi(args[atomii].c_str());
      atoms_per_species += natoms[speciesii][atomii];
    }
  }


  /*read box dimensions*/
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  L_x = atof(args[0].c_str());
  L_y = atof(args[1].c_str());
  L_z = atof(args[2].c_str());

//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  xlo = atof(args[0].c_str());
  xhi = atof(args[1].c_str());
  ylo = atof(args[2].c_str());
  yhi = atof(args[3].c_str());
  zlo = atof(args[4].c_str());
  zhi = atof(args[5].c_str());

  for(int timeii=0;timeii<n_timesteps;timeii++)
  {
    box_size[timeii].set(L_x,L_y,L_z);
    box_boundary[timeii][0].set(xlo, ylo, zlo);
    box_boundary[timeii][1].set(xhi, yhi, zhi);
  }
  count_atoms(natoms);

  create_molecules(natoms);

  if(trajtemplate_given){read_xyz(xyzfilename,trajtemplate);}
  else {read_xyz(xyzfilename);}

  wrapped = 1;
  /*calculate unwrapped trajectories*/
  cout << "\nUnwrapping trajectories.";
  unwrap();

}


/*method to read xyz file with log file*/
/*format is <xyzfilename> <logfilename> <templatefilename (optional)>*/

/*----------------------------------------------------------------------------------*/
void System::xyz_prep_withlog(vector<string> file_in, string fileline)
{


  string line;
  string args [ARGMAX];
  int n_args;
  int ** natoms;						//array of number of atoms of each type in each species
  int speciesii;						//species type index
  int typeii, argii;
  int atomii;							//atom type index
  float L_x, L_y, L_z;						//temporary storage for box dimensions
  float xlo, xhi, ylo, yhi, zlo, zhi;
  string xyzfilename;
  bool trajtemplate_given=0;			//does user specify a trajectory template?
  string trajtemplate, logfilename;
  int xlo_pos,ylo_pos,zlo_pos,xhi_pos,yhi_pos,zhi_pos;
  bool found_bounds=0;

  n_args = tokenize(fileline, args);
  xyzfilename = args[0];
  if(n_args<2)
  {
	Control::throw_error("No log file specified.", true);
  }
  else {logfilename = args[1];}

  ifstream ifile(xyzfilename.c_str());
  if (!ifile)
  {
	Control::throw_error("XYZ Trajectory file not found.", true);
  }
  ifstream ifile1(logfilename.c_str());
  if (!ifile1)
  {
	Control::throw_error("Log file not found.", true);
  }


  if(n_args>2)
  {
    trajtemplate_given=1;
    trajtemplate = args[2];
  }

  /*read in species names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  n_species = n_args/2;
  atoms_per_species = new int [n_species];
  species_name = new string [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    species_name[speciesii] = args[speciesii*2];
    atoms_per_species[speciesii]=0;
  }
  n_molecules = new int[n_species];	//allocate memory for molecule-number array
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    n_molecules[speciesii] = atoi(args[speciesii*2+1].c_str());
  }

  natoms = new int*[n_species];			//allocate memory for atom-number-per-molecule array

  /*read in atom type names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_atomtypes = tokenize(line, args);
  atomtype_name = new string [n_atomtypes];
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    atomtype_name[typeii] = args[typeii];
  }

  /*read in atomtype counts for each species*/
  for(speciesii=0; speciesii<n_species; speciesii++)
  {
    natoms[speciesii] = new int [n_atomtypes];
    line = "";
    for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
//    getline(*file_in,line);
    line = Control::read_line();
    n_args = tokenize(line, args);
    if(n_args!=n_atomtypes){cout << "Error: number of atom counts ("<<n_args<<") does not match number of atom types("<<n_atomtypes<<").\n";exit(1);}
    for(atomii=0; atomii<n_atomtypes; atomii++)
    {
      natoms[speciesii][atomii] = atoi(args[atomii].c_str());
      atoms_per_species += natoms[speciesii][atomii];
    }
  }


  /*read log file to get box dimensions*/

  ifstream logfile(logfilename.c_str());
  if (!logfile.is_open()){cout << "Error opening log file " << logfilename << ".\n"; exit(1);}
  while(!logfile.eof())
  {
    line = "";
    getline(logfile,line);
    n_args = tokenize(line, args);
    if(n_args==0){continue;}
    else if(args[0]=="Step")
    {
      if(!in_string_array(args,n_args,"Xlo")||!in_string_array(args,n_args,"Xhi")||!in_string_array(args,n_args,"Ylo")||!in_string_array(args,n_args,"Yhi")||!in_string_array(args,n_args,"Zlo")||!in_string_array(args,n_args,"Zhi"))
      {
      }
      else
      {
	xlo_pos=find_in_string_array(args,n_args,"Xlo");
	xhi_pos=find_in_string_array(args,n_args,"Xhi");
	ylo_pos=find_in_string_array(args,n_args,"Ylo");
	yhi_pos=find_in_string_array(args,n_args,"Yhi");
	zlo_pos=find_in_string_array(args,n_args,"Zlo");
	zhi_pos=find_in_string_array(args,n_args,"Zhi");
	line = "";
	getline(logfile,line);
	n_args = tokenize(line, args);
	xlo = atof(args[xlo_pos].c_str());
	xhi = atof(args[xhi_pos].c_str());
	ylo = atof(args[ylo_pos].c_str());
	yhi = atof(args[yhi_pos].c_str());
	zlo = atof(args[zlo_pos].c_str());
	zhi = atof(args[zhi_pos].c_str());
	found_bounds=1;
	L_x = xhi-xlo;
	L_y = yhi-ylo;
	L_z = zhi-zlo;

	for(int timeii=0;timeii<n_timesteps;timeii++)
	{
	  box_size[timeii].set(L_x,L_y,L_z);
	  box_boundary[timeii][0].set(xlo, ylo, zlo);
	  box_boundary[timeii][1].set(xhi, yhi, zhi);
	}
	break;
      }
    }
  }

  if(!found_bounds)
  {
    cout << "Error: Thermo data not found in log file/\n";
    exit(0);
  }

  count_atoms(natoms);

  create_molecules(natoms);

  if(trajtemplate_given){read_xyz(xyzfilename,trajtemplate);}
  else {read_xyz(xyzfilename);}

  wrapped = 1;
  /*calculate unwrapped trajectories*/
  cout << "\nUnwrapping trajectories.";
  unwrap();

}


/*----------------------------------------------------------------------------------*/

/*Method to read in a LAMMPS xyz file.  The file must be formatted in an ordered manner if this method is to work correctly.  In particular the method assumes that the xyz file is grouped by species, with the species appearing in the order in which they were entered by the user.  This approach could be replaced by a method to read in an enhanced xyz file that stores molecular data.*/

void System::read_xyz(string xyzfilename)
{
  int file_atoms;			//variable to store the number of atoms listed in xyz file
  string trash;				//create garbage string
  int timestepii;			//index over timesteps
  int speciesii;			//index over species
  int moleculeii;			//index over molecules
  int atomii;				//index over atoms within molecule
  int type;				//type to pass to molecule
  Coordinate coordinate;		//coordinate to pass to molecule
  float x;				//temporary storage for x coordinate
  float y;				//temporary storage for y coordinate
  float z;				//temporary storage for z coordinate
  int * n_typeii;			//array of indices to track how many of each type have been passed to particular molecule
  int typeii;				//index over elements of above aray
  int timetally=0;

  ifstream filexyz(xyzfilename.c_str());
  ifstream * fileobject = &filexyz;

  n_typeii = new int [n_atomtypes];
  cout << "\nReading a " << n_timesteps <<" timestep trajectory of " << n_atoms << " atoms.\n";

  for(timestepii=0; timestepii<n_timesteps; timestepii++)
  {
    *fileobject >> file_atoms;		//read in number of atoms from file
    if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
    {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
    }

    *fileobject >> trash;				//read out trash line: "atoms"
    getline(*fileobject,trash);

    for(speciesii=0; speciesii<n_species; speciesii++)
    {
      for(moleculeii=0; moleculeii<n_molecules[speciesii]; moleculeii++)
      {
        for(typeii=0;typeii<n_atomtypes;typeii++)
        {n_typeii[typeii]=0;}  //initiate type count array to zero at start of each molecule

        for(atomii=0; atomii < ((molecules[speciesii][moleculeii]).atomcount());atomii++)
        {
	      *fileobject >> type >> x >> y >> z;
          if(type > n_atomtypes)
          {
            cout << "Atom type in trajectory file out of range!\n";
            exit(1);
          }
         coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
         (molecules[speciesii][moleculeii]).set_coordinate(type-1,n_typeii[type-1],coordinate,timestepii);	//send coordinates to atom
         n_typeii[type-1]++;		//increment count of atoms of this type
       }
     }
   }
    print_progress(++timetally, n_timesteps);
  }

    *fileobject >> trash;
    if(trash !="")
    {
        cout << "\n\nWARNING: trajectory file contains more frames than were used."<<endl;cout.flush();
    }




  (*fileobject).close();
  delete [] n_typeii;
}


/*----------------------------------------------------------------------------------*/

/*Method to read in a LAMMPS xyz file.  Atoms within every molecule of each species must be ordered the same, but order of molecules is given by extra file*/

void System::read_xyz(string xyzfilename, string structure_filename)
{
  int file_atoms;			//variable to store the number of atoms listed in xyz fileng
  int timestepii;			//index over timesteps
  int speciesii;			//index over species
  int moleculeii;			//index over molecules
  int atomii;				//index over atoms within molecule
  Coordinate coordinate;		//coordinate to pass to molecule
  float x;				//temporary storage for x coordinate
  float y;				//temporary storage for y coordinate
  float z;				//temporary storage for z coordinate
  int * n_typeii;			//array of indices to track how many of each type have been passed to particular molecule
  int typeii;				//index over elements of above aray
  int timetally=0;
  string line;
  string args [ARGMAX];
  int n_args;
  int n_moleculeblocks=0;
  int * moleculeblock_type;
  int * moleculeblock_size;
  int moleculeblockii, argii;
  int ** atomorder;			//store order of atoms within molecules of each species
  int * moleculecount;

  moleculecount = new int [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    moleculecount[speciesii]=0;
  }

  ifstream filexyz(xyzfilename.c_str());
  ifstream structurefile(structure_filename.c_str());

  atomorder = new int * [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    atomorder[speciesii] = new int [molecules[speciesii][0].atomcount()];
  }


  n_typeii = new int [n_atomtypes];
  cout << "\nReading a " << n_timesteps <<" timestep trajectory of " << n_atoms << " atoms.\n";



  /*read in file containing information on trajectory file structure*/
  while(!structurefile.eof())
  {
    getline(structurefile, line);
//    line=Control::replace_constants(line);
    n_moleculeblocks++;
  }
  structurefile.seekg(0,ios::beg);	//go back to beginning of file
  structurefile.clear();
  moleculeblock_type = new int [n_moleculeblocks];
  moleculeblock_size = new int [n_moleculeblocks];

  moleculeblockii=0;
  while(!structurefile.eof())
  {
    getline(structurefile, line);
//    line=Control::replace_constants(line);
    n_args = tokenize(line, args);
    if(n_args==0){continue;}
    else if(n_args!=2)
    {
      cout << "Error: incorrect number of arguments in structure file.  Each line should consist of 2 arguments: a species and a number of lines.\n";
      exit(1);
    }
    else
    {
      moleculeblock_type[moleculeblockii] = show_species_index(args[0]);
      if(moleculeblock_type[moleculeblockii]==-1){cout<<"Error: species "<<args[0]<<" not found.\n";exit(1);}
      moleculeblock_size[moleculeblockii] = atoi(args[1].c_str());
      moleculeblockii++;
    }
  }

  /*determine atom order template for each species*/
  getline(filexyz, line);
  n_args = tokenize(line, args);
  file_atoms = atoi(args[0].c_str());
  if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
  {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
  }

  getline(filexyz, line); 		//read out garbage header line
  for(moleculeblockii=0;moleculeblockii<n_moleculeblocks;moleculeblockii++)
  {
    for(moleculeii=0;moleculeii<moleculeblock_size[moleculeblockii];moleculeii++)
    {
      line = "";
      for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
      for(atomii=0;atomii<molecules[moleculeblock_type[moleculeblockii]][0].atomcount();atomii++)
      {
        getline(filexyz, line);
        if(moleculeii==0)
        {
          n_args = tokenize(line, args);
          atomorder[moleculeblock_type[moleculeblockii]][atomii]=show_atomtype_index(args[0]);
        }
      }
    }
  }


  /*do data readout based on above template*/
  filexyz.seekg(0,ios::beg);		//go back to beginning of file and prepare for data readout
  filexyz.clear();
  cout << filexyz.eof()<<"\n";

  for(timestepii=0; timestepii<n_timesteps; timestepii++)
  {
    for(speciesii=0;speciesii<n_species;speciesii++)
    {
      moleculecount[speciesii]=0;
    }
    /*parse section header*/

    line = "";
    for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
    getline(filexyz, line);
    n_args = tokenize(line, args);
    file_atoms = atoi(args[0].c_str());
    if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
    {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
    }
    getline(filexyz, line); 								//read out trash line: "atoms"

    for(moleculeblockii=0;moleculeblockii<n_moleculeblocks;moleculeblockii++)
    {
      for(moleculeii=0;moleculeii<moleculeblock_size[moleculeblockii];moleculeii++)
      {
        for(typeii=0;typeii<n_atomtypes;typeii++) {n_typeii[typeii]=0;}  //initiate type count array to zero at start of each molecule
        line = "";
        for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
        for(atomii=0;atomii<molecules[moleculeblock_type[moleculeblockii]][0].atomcount();atomii++)
        {
          getline(filexyz, line);
          n_args = tokenize(line, args);
          x = atof(args[1].c_str());
          y = atof(args[2].c_str());
          z = atof(args[3].c_str());
          coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
          (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_coordinate(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
          n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]]++;		//increment count of atoms of this type
        }
        moleculecount[moleculeblock_type[moleculeblockii]]++;
      }
    }
    print_progress(++timetally, n_timesteps);
  }

  filexyz.close();
  delete [] n_typeii;
}

#define MOLECULE_SIZE 5
#define SPECIES_SIZE 4
#define TYPE_SIZE 6



/*------------------------Methods to read custom LAMMPS trajectories --------------------- ------*/


void System::custom_prep(vector<string> file_in, string fileline)
{
    /** Prepare to read in custom LAMMPS trajectory
    * @param file_in - user-specified input file
    * @param fileline - line containing location of trajectory file(s)
    * @author David S. Simmons
    **/
  string line;
  string args [ARGMAX];
  int n_args;
  int ** natoms;		//array of number of atoms of each type in each species
  int speciesii;		//species type index
  int typeii, argii;
  int atomii;							//atom type index
  string customfilename;
  bool trajtemplate_given=0;			//does user specify a trajectory template?
  string trajtemplate;

  n_args = tokenize(fileline, args);
  customfilename = args[0];
  if(n_args>1)
  {
    trajtemplate_given=1;
    trajtemplate = args[1];
  }

  ifstream ifile(customfilename.c_str());
  if (!ifile)
  {
	Control::throw_error("Custom Trajectory file not found.", true);
  }

  /*read in species names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++)
  {
    args[argii]="";
  }
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_args = tokenize(line, args);
  n_species = n_args/2;
  atoms_per_species = new int [n_species];
  species_name = new string [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    species_name[speciesii] = args[speciesii*2];
    atoms_per_species[speciesii]=0;
  }
  n_molecules = new int[n_species];	//allocate memory for molecule-number array
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    n_molecules[speciesii] = atoi(args[speciesii*2+1].c_str());
  }

  natoms = new int*[n_species];			//allocate memory for atom-number-per-molecule array

  /*read in atom type names*/
  line = "";
  for(argii=0;argii<ARGMAX;argii++)
  {
    args[argii]="";
  }
//  getline(*file_in,line);
//  line=Control::replace_constants(line);
  line = Control::read_line();
  n_atomtypes = tokenize(line, args);
  atomtype_name = new string [n_atomtypes];
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    atomtype_name[typeii] = args[typeii];
  }

  /*read in atomtype counts for each species*/
  for(speciesii=0; speciesii<n_species; speciesii++)
  {
    natoms[speciesii] = new int [n_atomtypes];
    line = "";
    for(argii=0;argii<ARGMAX;argii++)
    {
      args[argii]="";
    }
//    getline(*file_in,line);
    line = Control::read_line();
    n_args = tokenize(line, args);
    if(n_args!=n_atomtypes)
    {
      cout << "Error: number of atom counts ("<<n_args<<") does not match number of atom types("<<n_atomtypes<<").\n";
      exit(1);
    }
    for(atomii=0; atomii<n_atomtypes; atomii++)
    {
      natoms[speciesii][atomii] = atoi(args[atomii].c_str());
      atoms_per_species += natoms[speciesii][atomii];
    }
  }

  count_atoms(natoms);

  create_molecules(natoms);

  if(trajtemplate_given)
  {
    read_custom(customfilename,trajtemplate);
  }
  else
  {
    read_custom(customfilename);
  }


}

/*----------------------------------------------------------------------------------*/
/* Method to do a better job at comparing float values, since sometimes AMDAT will say box
dimensions are changing when they really aren't.
Taken from here: (which includes a good description of the problem and offers several solutions) 
http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
Originally called "AlmostEqualRelative" on that site.
*/
bool System::floatCompare(float A, float B)
{
	float maxRelDiff = FLT_EPSILON;
    // Calculate the difference.
    float diff = fabs(A - B);
    A = fabs(A);
    B = fabs(B);
    // Find the largest
    float largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}


/*----------------------------------------------------------------------------------*/

/*Method to read in a LAMMPS custom file.  The file must be formatted in an ordered manner if this method is to work correctly.  In particular the method assumes that the xyz file is grouped by species, with the species appearing in the order in which they were entered by the user.  This approach could be replaced by a method to read in an enhanced xyz file that stores molecular data.*/

void System::read_custom(string xyzfilename)
{
  /** Read in custom LAMMPS trajectory file
  * @param xyzfilename location of trajectory file
  * @author David S. Simmons
  **/

  int file_atoms;			//variable to store the number of atoms listed in xyz file
  string trash;				//create garbage string
  int timestepii;			//index over timesteps
  int speciesii;			//index over species
  int moleculeii;			//index over molecules
  int atomii;				//index over atoms within molecule
  int type;				//type to pass to molecule
  Coordinate coordinate;		//coordinate to pass to molecule
  float x;				//temporary storage for x coordinate
  float y;				//temporary storage for y coordinate
  float z;				//temporary storage for z coordinate
  int * n_typeii;			//array of indices to track how many of each type have been passed to particular molecule
  int typeii;				//index over elements of above aray
  int timetally=0;
  float xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz;
  int n_columns;		//number of columns of atom data in custom dump file
  bool r_provided, rs_provided, ru_provided, rsu_provided;		//booleans specifying whether a complete set of wrapped, scaled wraped, unwrapped, and scaled unwrapped coordinates are provided by the trajectory file
  bool i_provided;		//specifies whether image indices are provided for unwrapping
  bool v_provided;	//bool specifying whether complete velocity is provided by trajectory file
  bool mass_provided;	//bool specifying whether mass is provided by trajectory file
  bool read_r, read_rs, read_ru, read_rsu, read_i;		//bools specifying which type of coordinates to read in
  int x_position, y_position, z_position, xs_position, ys_position, zs_position, xu_position, yu_position, zu_position, xsu_position, ysu_position, zsu_position, ix_position, iy_position, iz_position, type_position, vx_position, vy_position, vz_position, mass_position;
  bool calc_wrapped;

  string line;
  string args [ARGMAX];
  int n_args;

  ifstream filexyz(xyzfilename.c_str());
  ifstream * fileobject = &filexyz;

  n_typeii = new int [n_atomtypes];
  cout << "\nReading a " << n_timesteps <<" timestep trajectory of " << n_atoms << " atoms.\n";

  for(timestepii=0; timestepii<n_timesteps; timestepii++)
  {
    /*read in header information from custom trajectory file*/
     line = "";
     getline(*fileobject,line);		//read in "ITEM: TIMESTEP" line
     getline(*fileobject,line);		//read in timestep line
     getline(*fileobject,line);		//read in "ITEM: NUMBER OF ATOMS" line
     line = "";
     getline(*fileobject,line);		//read in number of atoms
     n_args = tokenize(line, args);
     file_atoms=atoi(args[0].c_str());
    if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
    {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
    }

    /*read in box bounds from trajectory file*/
    getline(*fileobject,line);		//read in "ITEM: BOX BOUNDS..." line
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    xlo = atof(args[0].c_str());
    xhi = atof(args[1].c_str());
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    ylo = atof(args[0].c_str());
    yhi = atof(args[1].c_str());
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    zlo = atof(args[0].c_str());
    zhi = atof(args[1].c_str());
    line = "";

//    if(timestepii==0)
//    {
      /*set box size*/
      Lx = xhi - xlo;
      Ly = yhi - ylo;
      Lz = zhi - zlo;
      box_size[timestepii].set(Lx,Ly,Lz);
      box_boundary[timestepii][0].set(xlo, ylo, zlo);
      box_boundary[timestepii][1].set(xhi, yhi, zhi);
//    }
//    else
//    {
//      if(np) /*if non-constant volume, set box size at each frame*/
//      {
//	Lx = xhi - xlo;
//	Ly = yhi - ylo;
//	Lz = zhi - zlo;
//	box_size[timestepii].set(Lx,Ly,Lz);
//	box_boundary[timestepii][0].set(xlo, ylo, zlo);
//	box_boundary[timestepii][1].set(xhi, yhi, zhi);
//    }
//      else /*if constant volume, return error if file box bounds do not equal previous bounds*/
      if(!np)
      {
//	if(box_boundary[0][0].show_x()!=xlo||box_boundary[0][1].show_x()!=xhi||box_boundary[0][0].show_y()!=ylo||box_boundary[0][1].show_y()!=yhi||box_boundary[0][0].show_z()!=zlo||box_boundary[0][1].show_z()!=zhi)
        if (!(floatCompare(box_boundary[0][0].show_x(), xlo)&&floatCompare(box_boundary[0][1].show_x(), xhi)&&floatCompare(box_boundary[0][0].show_y(), ylo)&&floatCompare(box_boundary[0][1].show_y(), yhi)&&floatCompare(box_boundary[0][0].show_z(), zlo)&&floatCompare(box_boundary[0][1].show_z(), zhi)))
        {
            cout << "Error: The box boundaries provided in the custom file are not constant. For varying-volume trajectory, please select system_np system type.\n";
            exit(0);
        }
      }
 //   }


    /*read in and parse line specifying data types in custom dump file*/
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    n_columns = n_args - 2;	//determine number of columns of atom data

    /*If this is the initial timestep, determine which types of trajectory data are to be read in and which columns they are in*/
    if(timestepii==0)
    {
      /*check which coordinate types are provided*/
      r_provided = in_string_array(args,n_args,"x")&&in_string_array(args,n_args,"y")&&in_string_array(args,n_args,"z");		//wrapped
      rs_provided = in_string_array(args,n_args,"xs")&&in_string_array(args,n_args,"ys")&&in_string_array(args,n_args,"zs");		//wrapped, scaled
      ru_provided = in_string_array(args,n_args,"xu")&&in_string_array(args,n_args,"yu")&&in_string_array(args,n_args,"zu");		//unwrapped
      rsu_provided = in_string_array(args,n_args,"xsu")&&in_string_array(args,n_args,"ysu")&&in_string_array(args,n_args,"zsu");	//unwrapped, scaled
      i_provided = in_string_array(args,n_args,"ix")&&in_string_array(args,n_args,"iy")&&in_string_array(args,n_args,"iz");		//image index

      /*check if velocities are provided; if so, note their columns*/
      v_provided = in_string_array(args,n_args,"vx")&&in_string_array(args,n_args,"vy")&&in_string_array(args,n_args,"vz");
      if(v_provided)
      {
	vx_position = find_in_string_array(args,n_args,"vx")-2;
	vy_position = find_in_string_array(args,n_args,"vz")-2;
	vz_position = find_in_string_array(args,n_args,"vz")-2;
      }

      /*Check if mass is provided; if so, note its column*/
      mass_provided = in_string_array(args,n_args,"mass");
      if(mass_provided) mass_position = find_in_string_array(args,n_args,"mass")-2;

      /*figure out which coordinate types to read in*/
      if(!r_provided&&!rs_provided&&!ru_provided&&!rsu_provided)	//return error if no complete set of coordinates is provided
      {
	cout << "Error: no complete set of coordinates provided in trajectory file.\n";
	exit(0);
      }
      else								//otherwise, determine which coordinate types to read in and their location
      {
	calc_wrapped=true;
	read_r=read_rs=read_ru=read_rsu=read_i=false;
	/*Only read in at most one type of wrapped coordinates - unscaled or scaled; always choose unscaled if it is provided*/
	if(r_provided)
	{
	  read_r=true;
	  calc_wrapped=false;
	  x_position = find_in_string_array(args,n_args,"x")-2;
	  y_position = find_in_string_array(args,n_args,"y")-2;
	  z_position = find_in_string_array(args,n_args,"z")-2;
	}
	else if(rs_provided)
	{
	  read_rs=true;
	  calc_wrapped=false;
	  xs_position = find_in_string_array(args,n_args,"xs")-2;
	  ys_position = find_in_string_array(args,n_args,"ys")-2;
	  zs_position = find_in_string_array(args,n_args,"zs")-2;
	}

	/*Only read in at most one type of unwrapped coordinates - unscaled or scaled; always choose unscaled if it is provided*/
	if(ru_provided)
	{
	  read_ru=true;
	  xu_position = find_in_string_array(args,n_args,"xu")-2;
	  yu_position = find_in_string_array(args,n_args,"yu")-2;
	  zu_position = find_in_string_array(args,n_args,"zu")-2;
	}
	else if(rsu_provided)
	{
	  read_rsu=true;
	  xsu_position = find_in_string_array(args,n_args,"xsu")-2;
	  ysu_position = find_in_string_array(args,n_args,"ysu")-2;
	  zsu_position = find_in_string_array(args,n_args,"zsu")-2;
	}
	else if(i_provided&&(r_provided||rs_provided))
	{
	  read_i=true;
	  ix_position = find_in_string_array(args,n_args,"ix")-2;
	  iy_position = find_in_string_array(args,n_args,"iy")-2;
	  iz_position = find_in_string_array(args,n_args,"iz")-2;
	}

	if(np)
	{
	  if(!ru_provided&&!rsu_provided&&!i_provided)
	  {
	    cout << "Warning:  unwrapping of coordinates without image indices for variable-volume trajectories is incorrect. All calculations depending on unwrapped coordinates are therefore incorrect!\n";
	  }
	}
      }

      /*Find and store position of type column; if not present, return error*/
      if(in_string_array(args,n_args,"type"))
      {
	type_position=find_in_string_array(args,n_args,"type")-2;
      }
      else
      {
	cout << "Error: atom type data not provided.\n";
	exit(0);
      }
    }

    /*loop over atoms at this time*/
    for(speciesii=0; speciesii<n_species; speciesii++)
    {
      for(moleculeii=0; moleculeii<n_molecules[speciesii]; moleculeii++)
      {
        for(typeii=0;typeii<n_atomtypes;typeii++)
        {n_typeii[typeii]=0;}  //initiate type count array to zero at start of each molecule

        for(atomii=0; atomii < ((molecules[speciesii][moleculeii]).atomcount());atomii++)
        {
	  line = "";
	  getline(*fileobject,line);
	  n_args = tokenize(line, args);


	  /*read type and check whether it is valid*/
	  type = show_atomtype_index(args[type_position]);

	  if(type > n_atomtypes)
          {
            cout << "Atom type " << type << " in trajectory file out of range!\n";
            exit(1);
          }

	  if(read_r) /*read wrapped, unscaled coordinates if appropriate*/
	  {
	    x = atof(args[x_position].c_str());
	    y = atof(args[y_position].c_str());
	    z = atof(args[z_position].c_str());
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_coordinate(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_rs) /*read wrapped, scaled coordinates if appropriate*/
	  {
	    x = box_size[timestepii].show_x()*atof(args[xs_position].c_str())+box_boundary[timestepii][0].show_x();
	    y = box_size[timestepii].show_y()*atof(args[ys_position].c_str())+box_boundary[timestepii][0].show_y();
	    z = box_size[timestepii].show_z()*atof(args[zs_position].c_str())+box_boundary[timestepii][0].show_z();
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_coordinate(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }

	  if(read_i)
	  {
	    x=x+atof(args[ix_position].c_str())*box_size[timestepii].show_x();
	    y=y+atof(args[iy_position].c_str())*box_size[timestepii].show_y();
	    z=z+atof(args[iz_position].c_str())*box_size[timestepii].show_z();
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_unwrapped(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_ru)	/*read unwrapped, scaled coordinates if appropriate*/
	  {
	    x = atof(args[xu_position].c_str());
	    y = atof(args[yu_position].c_str());
	    z = atof(args[zu_position].c_str());
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_unwrapped(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_rsu)	/*read unwrapped, scaled coordinates if appropriate*/
	  {
	    x = box_size[timestepii].show_x()*atof(args[xsu_position].c_str())+box_boundary[timestepii][0].show_x();;
	    y = box_size[timestepii].show_y()*atof(args[ysu_position].c_str())+box_boundary[timestepii][0].show_y();;
	    z = box_size[timestepii].show_z()*atof(args[zsu_position].c_str())+box_boundary[timestepii][0].show_z();;
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_unwrapped(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }

	  if(calc_wrapped) /*calculate wrapped from unwrapped coordinates, if necessary*/
	  {
	    coordinate -= box_size[timestepii]*((coordinate-box_boundary[timestepii][0])/box_size[timestepii]).coord_floor();
	    (molecules[speciesii][moleculeii]).set_coordinate(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }

          if(v_provided)	/*read velocities if they are provided*/
	  {
	    x = atof(args[vx_position].c_str());
	    y = atof(args[vy_position].c_str());
	    z = atof(args[vz_position].c_str());
	    coordinate.set(x,y,z);		//store velocities temporarily in coordinate object
	    (molecules[speciesii][moleculeii]).set_velocity(type,n_typeii[type],coordinate,timestepii);	//send coordinates to atom
	  }
         if(mass_provided)
	 {
	   molecules[speciesii][moleculeii].show_atom_trajectory(type,n_typeii[type])->set_mass(atof(args[mass_position].c_str()));
	 }

         n_typeii[type]++;		//increment count of atoms of this type
       }
     }
   }
    print_progress(++timetally, n_timesteps);
  }

  (*fileobject).close();
  delete [] n_typeii;

  /*generate wrapped or unwrapped coordinates as necessary*/
  if(!ru_provided&&!rsu_provided&&!i_provided)
  {
    /*calculate unwrapped trajectories*/
    cout << "\nUnwrapping trajectories.";
    unwrap();
  }
  else
  {
    unwrapped = 1;
  }
//  if(!r_provided&&!rs_provided)
//  {
//    /*calculate wrapped trajectories*/
//    cout << "\nWrapping trajectories.";
//    wrap();
//  }
//  else
//  {
    wrapped = 1;
//  }
}




/*--------------------------------------------------------------------------------------------*/


/*Method to read in a LAMMPS custom file.  Atoms within every molecule of each species must be ordered the same, but order of molecules is given by extra template file*/

void System::read_custom(string xyzfilename, string structure_filename)
{
  /** Read in custom LAMMPS trajectory file with a template file specifying order of molecules in trajectory
  * @param xyzfilename location of trajectory file
  * @param structure_filename location of template file
  * @author David S. Simmons
  **/
  int file_atoms;			//variable to store the number of atoms listed in xyz fileng
  int timestepii;			//index over timesteps
  int speciesii;			//index over species
  int moleculeii;			//index over molecules
  int atomii;				//index over atoms within molecule
  Coordinate coordinate;		//coordinate to pass to molecule
  float x;				//temporary storage for x coordinate
  float y;				//temporary storage for y coordinate
  float z;				//temporary storage for z coordinate
  int * n_typeii;			//array of indices to track how many of each type have been passed to particular molecule
  int typeii;				//index over elements of above aray
  int timetally=0;
  string line;
  string args [ARGMAX];
  int n_args;
  int n_moleculeblocks=0;
  int * moleculeblock_type;
  int * moleculeblock_size;
  int moleculeblockii, argii;
  int ** atomorder;			//store order of atoms within molecules of each species
  int * moleculecount;
  int type;
  float xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz;
  int n_columns;		//number of columns of atom data in custom dump file
  bool r_provided, rs_provided, ru_provided, rsu_provided;		//booleans specifying whether a complete set of wrapped, scaled wraped, unwrapped, and scaled unwrapped coordinates are provided by the trajectory file
  bool i_provided;		//specifies whether image indices are provided for unwrapping
  bool v_provided;	//bool specifying whether complete velocity is provided by trajectory file
  bool mass_provided;	//bool specifying whether mass is provided by trajectory file
  bool read_r, read_rs, read_ru, read_rsu, read_i;		//bools specifying which type of coordinates to read in
  int x_position, y_position, z_position, xs_position, ys_position, zs_position, xu_position, yu_position, zu_position, xsu_position, ysu_position, zsu_position, type_position, vx_position, vy_position, vz_position, ix_position, iy_position, iz_position, mass_position;
  bool calc_wrapped;

  moleculecount = new int [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    moleculecount[speciesii]=0;
  }

  ifstream filexyz(xyzfilename.c_str());
  ifstream * fileobject = &filexyz;
  ifstream structurefile(structure_filename.c_str());

  atomorder = new int * [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    atomorder[speciesii] = new int [molecules[speciesii][0].atomcount()];
  }

  n_typeii = new int [n_atomtypes];
  cout << "\nReading a " << n_timesteps <<" timestep trajectory of " << n_atoms << " atoms.\n";



  /*read in template file containing information on trajectory file structure*/
  while(!structurefile.eof())
  {
    getline(structurefile, line);
//    line=Control::replace_constants(line);
    n_moleculeblocks++;
  }
  structurefile.seekg(0,ios::beg);	//go back to beginning of file
  structurefile.clear();
  moleculeblock_type = new int [n_moleculeblocks];
  moleculeblock_size = new int [n_moleculeblocks];

  moleculeblockii=0;
  while(!structurefile.eof())
  {
    getline(structurefile, line);
//    line=Control::replace_constants(line);
    n_args = tokenize(line, args);
    if(n_args==0){continue;}
    else if(n_args!=2)
    {
      cout << "Error: incorrect number of arguments in structure file.  Each line should consist of 2 arguments: a species and a number of lines.\n";
      exit(1);
    }
    else
    {
      moleculeblock_type[moleculeblockii] = show_species_index(args[0]);
      if(moleculeblock_type[moleculeblockii]==-1){cout<<"Error: species "<<args[0]<<" not found.\n";exit(1);}
      moleculeblock_size[moleculeblockii] = atoi(args[1].c_str());
      moleculeblockii++;
    }
  }



   /*determine atom order template for each species*/
   /*header lines*/
    getline(*fileobject,line);		//read in "ITEM: TIMESTEP" line
    getline(*fileobject,line);		//read in timestep line
    getline(*fileobject,line);		//read in "ITEM: NUMBER OF ATOMS" line
    line = "";
    getline(*fileobject,line);		//read in number of atoms
    n_args = tokenize(line, args);
    file_atoms=atoi(args[0].c_str());
    if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
    {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
    }
    getline(*fileobject,line);		//read in "ITEM: BOX BOUNDS..." line
    getline(*fileobject,line);
    getline(*fileobject,line);
    getline(*fileobject,line);
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    /*Find and store position of type column; if not present, return error*/
    if(in_string_array(args,n_args,"type"))
    {
      type_position=find_in_string_array(args,n_args,"type")-2;
    }
    else
    {
      cout << "Error: atom type data not provided.\n";
      exit(0);
    }
  for(moleculeblockii=0;moleculeblockii<n_moleculeblocks;moleculeblockii++)
  {
    for(moleculeii=0;moleculeii<moleculeblock_size[moleculeblockii];moleculeii++)
    {
      line = "";
      for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
      for(atomii=0;atomii<molecules[moleculeblock_type[moleculeblockii]][0].atomcount();atomii++)
      {
        getline(filexyz, line);
        if(moleculeii==0)
        {
          n_args = tokenize(line, args);
	  type=show_atomtype_index(args[type_position]);
	  if(type > n_atomtypes)
          {
            cout << "Atom type " << type << " in trajectory file out of range!\n";
            exit(1);
          }
          atomorder[moleculeblock_type[moleculeblockii]][atomii]=type;
        }
      }
    }
  }


  /*do data readout based on above template*/
  filexyz.seekg(0,ios::beg);		//go back to beginning of file and prepare for data readout
  filexyz.clear();
//  cout << filexyz.eof()<<"\n";
  for(timestepii=0; timestepii<n_timesteps; timestepii++)
  {
    /*zero out count of molecules at this timestep*/
    for(speciesii=0;speciesii<n_species;speciesii++)
    {
      moleculecount[speciesii]=0;
    }

   /*read in header information from custom trajectory file*/
     line = "";
     getline(*fileobject,line);		//read in "ITEM: TIMESTEP" line
     getline(*fileobject,line);		//read in timestep line
     getline(*fileobject,line);		//read in "ITEM: NUMBER OF ATOMS" line
     line = "";
     getline(*fileobject,line);		//read in number of atoms
     n_args = tokenize(line, args);
     file_atoms=atoi(args[0].c_str());
    if(file_atoms!=n_atoms)		//check if the number of atoms listed in file is consistent with the molecule and atom counts given by the user
    {
      cout << "The number of atoms listed in the xyz file is inconsistent with user input: ";	//if not, give error...
      cout << file_atoms << " != " << n_atoms << "\n";
      exit(0);											//and terminate program.
    }

    /*read in box bounds from trajectory file*/
    getline(*fileobject,line);		//read in "ITEM: BOX BOUNDS..." line
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    xlo = atof(args[0].c_str());
    xhi = atof(args[1].c_str());
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    ylo = atof(args[0].c_str());
    yhi = atof(args[1].c_str());
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    zlo = atof(args[0].c_str());
    zhi = atof(args[1].c_str());
    line = "";

//    if(timestepii==0)
//    {
      /*set box size*/
      Lx = xhi - xlo;
      Ly = yhi - ylo;
      Lz = zhi - zlo;
      box_size[timestepii].set(Lx,Ly,Lz);
      box_boundary[timestepii][0].set(xlo, ylo, zlo);
      box_boundary[timestepii][1].set(xhi, yhi, zhi);
//    }
//    else
//    {
//      if(np) /*if non-constant volume, set box size at each frame*/
//      {
//	Lx = xhi - xlo;
//	Ly = yhi - ylo;
//	Lz = zhi - zlo;
//	box_size[timestepii].set(Lx,Ly,Lz);
//	box_boundary[timestepii][0].set(xlo, ylo, zlo);
//	box_boundary[timestepii][1].set(xhi, yhi, zhi);
//    }
//      else /*if constant volume, return error if file box bounds do not equal previous bounds*/
      if(!np)
      {
//	if(box_boundary[0][0].show_x()!=xlo||box_boundary[0][1].show_x()!=xhi||box_boundary[0][0].show_y()!=ylo||box_boundary[0][1].show_y()!=yhi||box_boundary[0][0].show_z()!=zlo||box_boundary[0][1].show_z()!=zhi)
        if (!(floatCompare(box_boundary[0][0].show_x(), xlo)&&floatCompare(box_boundary[0][1].show_x(), xhi)&&floatCompare(box_boundary[0][0].show_y(), ylo)&&floatCompare(box_boundary[0][1].show_y(), yhi)&&floatCompare(box_boundary[0][0].show_z(), zlo)&&floatCompare(box_boundary[0][1].show_z(), zhi)))
        {
            cout << "Error: The box boundaries provided in the custom file are not constant. For varying-volume trajectory, please select system_np system type.\n";
            exit(0);
        }
      }
 //   }


    /*read in and parse line specifying data types in custom dump file*/
    line = "";
    getline(*fileobject,line);
    n_args = tokenize(line, args);
    n_columns = n_args - 2;	//determine number of columns of atom data

    /*If this is the initial timestep, determine which types of trajectory data are to be read in and which columns they are in*/
    if(timestepii==0)
    {
      /*check which coordinate types are provided*/
      r_provided = in_string_array(args,n_args,"x")&&in_string_array(args,n_args,"y")&&in_string_array(args,n_args,"z");		//wrapped
      rs_provided = in_string_array(args,n_args,"xs")&&in_string_array(args,n_args,"ys")&&in_string_array(args,n_args,"zs");		//wrapped, scaled
      ru_provided = in_string_array(args,n_args,"xu")&&in_string_array(args,n_args,"yu")&&in_string_array(args,n_args,"zu");		//unwrapped
      rsu_provided = in_string_array(args,n_args,"xsu")&&in_string_array(args,n_args,"ysu")&&in_string_array(args,n_args,"zsu");	//unwrapped, scaled
      i_provided = in_string_array(args,n_args,"ix")&&in_string_array(args,n_args,"iy")&&in_string_array(args,n_args,"iz");		//image index

      /*check if velocities are provided; if so, note their columns*/
      v_provided = in_string_array(args,n_args,"vx")&&in_string_array(args,n_args,"vy")&&in_string_array(args,n_args,"vz");
      if(v_provided)
      {
	vx_position = find_in_string_array(args,n_args,"vx")-2;
	vy_position = find_in_string_array(args,n_args,"vz")-2;
	vz_position = find_in_string_array(args,n_args,"vz")-2;
      }

      /*Check if mass is provided; if so, note its column*/
      mass_provided = in_string_array(args,n_args,"mass");
      if(mass_provided) mass_position = find_in_string_array(args,n_args,"mass")-2;

      /*figure out which coordinate types to read in*/
      if(!r_provided&&!rs_provided&&!ru_provided&&!rsu_provided)	//return error if no complete set of coordinates is provided
      {
	cout << "Error: no complete set of coordinates provided in trajectory file.\n";
	exit(0);
      }
      else								//otherwise, determine which coordinate types to read in and their location
      {
	calc_wrapped=true;
	read_r=read_rs=read_ru=read_rsu=read_i=false;
	/*Only read in at most one type of wrapped coordinates - unscaled or scaled; always choose unscaled if it is provided*/
	if(r_provided)
	{
	  read_r=true;
	  calc_wrapped=false;
	  x_position = find_in_string_array(args,n_args,"x")-2;
	  y_position = find_in_string_array(args,n_args,"y")-2;
	  z_position = find_in_string_array(args,n_args,"z")-2;
	}
	else if(rs_provided)
	{
	  read_rs=true;
	  calc_wrapped=false;
	  xs_position = find_in_string_array(args,n_args,"xs")-2;
	  ys_position = find_in_string_array(args,n_args,"ys")-2;
	  zs_position = find_in_string_array(args,n_args,"zs")-2;
	}


	/*Only read in at most one type of unwrapped coordinates - unscaled or scaled; always choose unscaled if it is provided*/
	if(ru_provided)
	{
	  read_ru=true;
	  xu_position = find_in_string_array(args,n_args,"xu")-2;
	  yu_position = find_in_string_array(args,n_args,"yu")-2;
	  zu_position = find_in_string_array(args,n_args,"zu")-2;
	}
	else if(rsu_provided)
	{
	  read_rsu=true;
	  xsu_position = find_in_string_array(args,n_args,"xsu")-2;
	  ysu_position = find_in_string_array(args,n_args,"ysu")-2;
	  zsu_position = find_in_string_array(args,n_args,"zsu")-2;
	}
	else if(i_provided&&(r_provided||rs_provided))
	{
	  read_i=true;
	  ix_position = find_in_string_array(args,n_args,"ix")-2;
	  iy_position = find_in_string_array(args,n_args,"iy")-2;
	  iz_position = find_in_string_array(args,n_args,"iz")-2;
	}


	if(np)
	{
	  if(!ru_provided&&!rsu_provided&&!i_provided)
	  {
	    cout << "Warning:  unwrapping of coordinates without image indices for variable-volume trajectories is incorrect. All calculations depending on unwrapped coordinates are therefore incorrect!\n";
	  }
	}
      }

      /*Find and store position of type column; if not present, return error*/
      if(in_string_array(args,n_args,"type"))
      {
	type_position=find_in_string_array(args,n_args,"type")-2;
      }
      else
      {
	cout << "Error: atom type data not provided.\n";
	exit(0);
      }
    }

  //[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]
    /*loop over atoms at this time*/
    for(moleculeblockii=0;moleculeblockii<n_moleculeblocks;moleculeblockii++)
    {
      for(moleculeii=0;moleculeii<moleculeblock_size[moleculeblockii];moleculeii++)
      {
        for(typeii=0;typeii<n_atomtypes;typeii++) {n_typeii[typeii]=0;}  //initiate type count array to zero at start of each molecule
	line = "";
        for(argii=0;argii<ARGMAX;argii++){args[argii]="";}
        for(atomii=0;atomii<molecules[moleculeblock_type[moleculeblockii]][0].atomcount();atomii++)
        {
	  line = "";
	  getline(*fileobject,line);
	  n_args = tokenize(line, args);


	  /*read type and check whether it is valid*/
	  type = show_atomtype_index(args[type_position]);
	  if(type > n_atomtypes)
          {
            cout << "Atom type " << type << " in trajectory file out of range!\n";
            exit(1);
          }

	  if(read_r) /*read wrapped, unscaled coordinates if appropriate*/
	  {
	    x = atof(args[x_position].c_str());
	    y = atof(args[y_position].c_str());
	    z = atof(args[z_position].c_str());
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	   (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_coordinate(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_rs) /*read wrapped, scaled coordinates if appropriate*/
	  {
	    x = box_size[timestepii].show_x()*atof(args[xs_position].c_str())+box_boundary[timestepii][0].show_x();
	    y = box_size[timestepii].show_y()*atof(args[ys_position].c_str())+box_boundary[timestepii][0].show_y();
	    z = box_size[timestepii].show_z()*atof(args[zs_position].c_str())+box_boundary[timestepii][0].show_z();
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	   (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_coordinate(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }

	  if(read_i)
	  {
	    x=x+atof(args[ix_position].c_str())*box_size[timestepii].show_x();
	    y=y+atof(args[iy_position].c_str())*box_size[timestepii].show_y();
	    z=z+atof(args[iz_position].c_str())*box_size[timestepii].show_z();
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_unwrapped(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_ru)	/*read unwrapped, scaled coordinates if appropriate*/
	  {
	    x = atof(args[xu_position].c_str());
	    y = atof(args[yu_position].c_str());
	    z = atof(args[zu_position].c_str());
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_unwrapped(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }
	  else if(read_rsu)	/*read unwrapped, scaled coordinates if appropriate*/
	  {
	    x = box_size[timestepii].show_x()*atof(args[xsu_position].c_str())+box_boundary[timestepii][0].show_x();;
	    y = box_size[timestepii].show_y()*atof(args[ysu_position].c_str())+box_boundary[timestepii][0].show_y();;
	    z = box_size[timestepii].show_z()*atof(args[zsu_position].c_str())+box_boundary[timestepii][0].show_z();;
	    coordinate.set(x,y,z);		//store coordinates temporarily in coordinate object
	    (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_unwrapped(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }

	  if(calc_wrapped) /*calculate wrapped from unwrapped coordinates, if necessary*/
	  {
	    coordinate -= box_size[timestepii]*((coordinate-box_boundary[timestepii][0])/box_size[timestepii]).coord_floor();
	    (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_coordinate(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }

          if(v_provided)	/*read velocities if they are provided*/
	  {
	    x = atof(args[vx_position].c_str());
	    y = atof(args[vy_position].c_str());
	    z = atof(args[vz_position].c_str());
	    coordinate.set(x,y,z);		//store velocities temporarily in coordinate object
	    (molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]]).set_velocity(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]],coordinate,timestepii);	//send coordinates to atom
	  }

         if(mass_provided)
	 {
	   molecules[moleculeblock_type[moleculeblockii]][moleculecount[moleculeblock_type[moleculeblockii]]].show_atom_trajectory(atomorder[moleculeblock_type[moleculeblockii]][atomii],n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]])->set_mass(atof(args[mass_position].c_str()));
	 }

         n_typeii[atomorder[moleculeblock_type[moleculeblockii]][atomii]]++;		//increment count of atoms of this type
       }
       moleculecount[moleculeblock_type[moleculeblockii]]++;
     }
   }
    print_progress(++timetally, n_timesteps);
  }

  filexyz.close();
  delete [] n_typeii;

  /*generate wrapped or unwrapped coordinates as necessary*/
  if(!ru_provided&&!rsu_provided&&!i_provided)
  {
    /*calculate unwrapped trajectories*/
    cout << "\nUnwrapping trajectories.";
    unwrap();
  }
  else
  {
    unwrapped = 1;
  }
//  if(!r_provided&&!rs_provided)
//  {
//    /*calculate wrapped trajectories*/
//    cout << "\nWrapping trajectories.";
//    wrap();
//  }
//  else
//  {
    wrapped = 1;
//  }
}



/*--------------------------------------------------------------------------------------------*/


#define MOLECULE_SIZE 5
#define SPECIES_SIZE 4
#define TYPE_SIZE 6


#ifndef TACC
void System::xtc_prep(vector<string> file_in, string fileline)
{
  /** Prepare to read in xtc trajectory file
  * @param file_in user spectified input file
  * @param fileline line from input file specifying trajectory file location(s)
  * @author David S. Simmons
  **/
  string line;
  string args [ARGMAX];
  int n_args;
  string gro_file, xtc_file;
  int ** atomidentifier;

  n_args = tokenize(fileline, args);

  if(n_args!=2){cout<<"Error: incorrect number of files given for xtc trajectory.  File line must contain two files: an xtc file and a gro file.\n";exit(1);}

  xtc_file = args[0];
  gro_file = args[1];

  ifstream ifile(xtc_file.c_str());
  if (!ifile)
  {
	Control::throw_error("XTC Trajectory file not found.", true);
  }
  ifstream ifile1(gro_file.c_str());
  if (!ifile1)
  {
	Control::throw_error("GRO file not found.", true);
  }

  atomidentifier = read_gro(gro_file);
  read_xtc_format(xtc_file,atomidentifier);

  /*calculate unwrapped trajectories*/
  cout << "\nUnwrapping trajectories.";
  unwrap();

}
#endif

/*--------------------------------------------------------------------------------------------*/


int** System::read_gro(string filename)
{
  /** Read in .gro file to get structure and format data before reading .xtc trajectory file
  * @param fileline line from input file specifying trajectory file location(s)
  * @author David S. Simmons
  **/
  string line;
  string args [ARGMAX];
  int n_args;
  int ** natoms;						//array of number of atoms of each type in each species
  int ii, moleculeii,typeii;
  int molecule_size=MOLECULE_SIZE;
  int species_size=SPECIES_SIZE;
  int type_size=TYPE_SIZE;

  int current_molecule=0;
  string current_species, current_type;
  int species_match, type_match;
  int molecule_num;
  int * this_molecule_count;

  /*set maximum number of species and atom_types to be read in*/
  int max_species=256;
  int max_types=256;
  float x,y,z;

  int ** atomidentifier;

  /*allocate memory based on these maximums*/
  species_name = new string [max_species];
  n_molecules = new int [max_species];
  natoms = new int*[max_species];
  atomtype_name = new string [max_types];
  this_molecule_count = new int [max_types];

  /*open .gro file*/
  ifstream filegro(filename.c_str());

  /*read title line - not used*/
  getline(filegro, line);
  n_args = tokenize(line, args);

  /*read number of atoms*/
  getline(filegro, line);
  n_args = tokenize(line, args);
  n_atoms = atoi(args[0].c_str());

  atomidentifier = new int* [n_atoms];
  for(ii=0;ii<n_atoms;ii++)
  {
    atomidentifier[ii] = new int [4];
  }

  /*loop over coordinate lines to determine the species and atom types*/
  n_species=0;
  n_atomtypes=0;
  for(ii=0;ii<n_atoms;ii++)
  {
    /*read molecule index*/
    getline(filegro,line);
    moleculeii=atoi((line.substr(0,molecule_size)).c_str());	/*molecule index*/
    current_species=line.substr(molecule_size,species_size);
    n_args = tokenize(current_species, args);
    current_species=args[0];
    current_type=line.substr(molecule_size+species_size,type_size);
    n_args = tokenize(current_type, args);
    current_type = args[0];

    /*Update lists of species and molecules, as appropriate*/
    if(moleculeii>current_molecule)  //if this is a new molecule, assign to proper species
    {
      current_molecule=moleculeii;  //update current molecule
      species_match = show_species_index(current_species);  //check if the most recently read species matches any present species

      if(species_match==-1)
      {
        species_name[n_species]=current_species;	//add species to list
        n_molecules[n_species]=1;			//initiate count of molecules
        natoms[n_species]=new int[max_types];		//assign memory for atom type counts for this species
        for(typeii=0;typeii<max_types;typeii++)
        {
          natoms[n_species][typeii]=0;
        }
        n_species++;					//increment species count
      }
      else
      {
        n_molecules[species_match]++;			//update count of molecules
      }

      for(typeii=0;typeii<max_types;typeii++)
      {
        this_molecule_count[typeii]=0;
      }
      molecule_num=n_molecules[show_species_index(current_species)]-1;
    }

    /*Update lists of atoms, as appropriate*/
    if(species_match==-1) //if the current molecule is a new species, assign atom to proper type
    {
      type_match = show_atomtype_index(current_type);
      if(type_match==-1)
      {
        atomtype_name[n_atomtypes]=current_type;	//add atom type to list
        natoms[show_species_index(current_species)][n_atomtypes]++;		//increment count of this atom type in this species

        n_atomtypes++;					//increment atom type count
      }
      else
      {
        natoms[show_species_index(current_species)][type_match]++;
      }
    }

    //add this atom to list of atom addresses, and then increment the list of atoms in this molecule
    atomidentifier[ii][0] = show_species_index(current_species);
    atomidentifier[ii][1] = molecule_num;
    atomidentifier[ii][2] = show_atomtype_index(current_type);
    atomidentifier[ii][3] = this_molecule_count[atomidentifier[ii][2]];
    this_molecule_count[atomidentifier[ii][2]]++;

  }

  /*get box size*/
  getline(filegro,line);
  n_args = tokenize(line, args);
  x = atof(args[0].c_str());
  y = atof(args[1].c_str());
  z = atof(args[2].c_str());
  for(int timeii=0;timeii<n_timesteps;timeii++) box_size[timeii].set(x,y,z);


  /*create data structure based on the information obtained from this .gro file*/
  create_molecules(natoms);

  cout <<"\n" << n_species << " found in .gro file:";

  for(ii=0;ii<n_species;ii++)
  {
    cout<<"\n" << species_name[ii];
  }
  return atomidentifier;
}


/*----------------------------------------------------------------------------------*/

#ifndef TACC
void System::read_xtc_format(string filename, int** atomidentifier)
{
  int atomcount;
  XDRFILE * trajectory;
  int step;
  float time;
  matrix boxdim;
  rvec * coordinates;
  float *prec;
  float precision;
  int errorcode;
  int timeii,atomii;
  const char mode = 'r';
  Coordinate coordinate;

  trajectory = xdrfile_open(filename.c_str(),&mode);

  errorcode = read_xtc_natoms(const_cast<char *>(filename.c_str()),&atomcount);

  if(atomcount!=n_atoms)
  {
    cout<<"Error: Number of atoms given by .gro file (" << n_atoms<<") and .xtc file (" << atomcount << ") are not equal.\n";
    exit(1);
  }

  coordinates = new rvec [n_atoms];
  prec=&precision;

  //step = 0;
  //time = 0;


  for(timeii=0;timeii<n_timesteps;timeii++)
  {
    errorcode = read_xtc(trajectory,atomcount,&step,&time,boxdim,coordinates,prec);

    if(atomcount!=n_atoms)
    {
      cout<<"Error: Number of atoms given by .gro file (" << n_atoms<<") and .xtc file (" << atomcount << ") are not equal.\n";
      exit(1);
    }

    for(atomii=0;atomii<n_atoms;atomii++)
    {
      coordinate.set(coordinates[atomii][0],coordinates[atomii][1],coordinates[atomii][2]);		//store coordinates temporarily in coordinate object
      (molecules[atomidentifier[atomii][0]][atomidentifier[atomii][1]]).set_coordinate(atomidentifier[atomii][2],atomidentifier[atomii][3],coordinate,timeii);
    }


    //errorcode = xdrfile_close(trajectory);
  }
  delete[] coordinates;
}
#endif

/*----------------------------------------------------------------------------------*/



/*Determine times corresponding to timesteps*/
void System::create_timelist()
{
  int timeii=0;
  int expii = 1;
  int blockii;
  float block_starttime=0;

  timelist = new float [n_timesteps];
  timelist[timeii] = 0;

  for(blockii=0;blockii<n_exponentials;blockii++)
  {
	for(expii=1;expii<=n_exponential_steps;expii++)
	{
		timeii++;
		if(pow(exp_base,expii-1+first_exponent) <= expii)
		{
			timelist[timeii] = block_starttime+float(expii)*time_unit;
		}
		else
		{
			timelist[timeii] = block_starttime+float(floor(pow(exp_base,expii-1+first_exponent)))*time_unit;
		}
	}
	block_starttime = timelist[timeii];
  }
}



/*----------------------------------------------------------------------------------*/
Coordinate System::min_box_dimensions()const
{
  Coordinate min(0,0,0);
  min.smallest(box_size,n_timesteps);
  return min;
}




/*----------------------------------------------------------------------------------*/


int System::block(int blockii,int* timelist)const
{
	int blocksize,timeii,starttime;

	blocksize = n_exponential_steps+1;

	timelist = new int [blocksize];

	starttime = blockstart(blockii);

	for(timeii=0;timeii<=n_exponential_steps;timeii++)
	{
		timelist[timeii]=starttime+timeii;
	}

	return blocksize;
}

/*----------------------------------------------------------------------------------*/


void System::big_block(int blockii,int* timelist)const
{
	int blocksize,timeii,starttime;

	blocksize = n_exponential_steps+n_exponentials-blockii;

	starttime = blockstart(blockii);

	for(timeii=0;timeii<=n_exponential_steps;timeii++)
	{
		timelist[timeii]=starttime+timeii;
	//	cout << "\n" << timelist[timeii];
	}

	for(timeii=n_exponential_steps+1;timeii<blocksize;timeii++)
	{
		timelist[timeii]=starttime+n_exponential_steps*(timeii-n_exponential_steps+1);
	//	cout << "\n" << timelist[timeii];
	}
	//cout<<"\n";

	//exit (0);
}

/*----------------------------------------------------------------------------------*/


/*Method to count up atoms and arrive at total atom number*/
void System::count_atoms(int ** natoms)
{
  int speciesii;						//species type index
  int atomii;							//atom type index
  n_atoms = 0;

  for(speciesii=1; speciesii<=n_species; speciesii++)
  {
    for(atomii=1; atomii<=n_atomtypes; atomii++)
    {
      n_atoms += n_molecules[speciesii-1]*natoms[speciesii-1][atomii-1];
    }
  }
}


/*----------------------------------------------------------------------------------*/


/*Method to create molecule array and pass atomic structure down to molecules*/
void System::create_molecules(int ** natoms)
{
  int speciesii;				//species type index
  int moleculeii;				//index over molecules of species

  //create non-square array of molecule objects, segregated by type
  molecules = new Molecule * [n_species];
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    molecules[speciesii] = new Molecule[n_molecules[speciesii]];
    for(moleculeii=0; moleculeii<n_molecules[speciesii]; moleculeii++)	//create array of molecule objects
    {
      (molecules[speciesii][moleculeii]).set(n_atomtypes,natoms[speciesii],n_timesteps);	//tell molecule object how to create atom array
    }
  }
}


/*----------------------------------------------------------------------------------*/


//System::~System()
//{
//  clear_memory();
//}





/*----------------------------------------------------------------------------------*/

void System::clear_memory()
{
  int speciesii;
  for(speciesii=1;speciesii<n_species;speciesii++)
  {
    delete [] molecules[speciesii];
  }
  delete [] molecules;
  delete [] timelist;
  delete [] n_molecules;
}




/*----------------------------------------------------------------------------------*/


/*Method to command molecules to unwrap atom trajectories*/
void System::unwrap()
{
  int speciesii;
  int moleculeii;
  int progressii = 0;

  if(!unwrapped)
  {
    cout << "\nUnwrapping coordinates...\n";

    for(speciesii=0; speciesii<n_species; speciesii++)
    {
      for(moleculeii=0; moleculeii<n_molecules[speciesii]; moleculeii++)
      {
        progressii += (molecules[speciesii][moleculeii]).unwrap_atoms(box_size[0]);	//unwrap molecule and increment progress index
        print_progress(progressii,n_atoms);
      }
    }
    unwrapped = 1;		//mark system as unwrapped
  }
}



/*----------------------------------------------------------------------------------*/


/*Method to command molecules to wrap atom trajectories*/
void System::wrap()
{
  int speciesii;
  int moleculeii;
  int progressii = 0;

  if(!wrapped)
  {
    cout << "\nWrapping coordinates...\n";

    for(speciesii=0; speciesii<n_species; speciesii++)
    {
      for(moleculeii=0; moleculeii<n_molecules[speciesii]; moleculeii++)
      {
        progressii += (molecules[speciesii][moleculeii]).wrap_atoms(box_size,box_boundary);	//wrap molecule and increment progress index
        print_progress(progressii,n_atoms);
      }
    }
    wrapped = 1;		//mark system as wrapped
  }
}


/*----------------------------------------------------------------------------------*/



/*calculate array of times corresponding to displacement timesteps*/
float * System::displacement_times() const
{
	int timeii;

	float * timegap;

	timegap = new float [n_timegaps];

	for(timeii=0; timeii<n_exponential_steps;timeii++)
	{
		timegap[timeii] = timelist[timeii+int(frt)]-timelist[int(frt)];
	}

	for(timeii=1; timeii<n_exponentials; timeii++)
	{
		timegap[timeii+n_exponential_steps-1] = timelist[n_exponential_steps*timeii+int(frt)]-timelist[int(frt)];
	}
	if(!frt) {timegap[n_timegaps-1] = timelist[n_timesteps-1]-timelist[0];}

	return timegap;

}

/*-----------------------------------------------------------------------------------*/


/*calculate time corresponding to given displacement timestep*/
float System::displacement_times(int timeii) const
{
	float timegap;

	if(timeii>=0&&timeii<n_exponential_steps)
	{
	  timegap = timelist[timeii+int(frt)]-timelist[int(frt)];
	}
	else if(timeii>=n_exponential_steps&&timegap<n_timegaps-1+int(frt))
	{
	  timegap = timelist[n_exponential_steps*timeii+int(frt)]-timelist[int(frt)];
	}
	else if(timegap=n_timegaps-1+int(frt)&&!frt)
	{
	  timegap = timelist[n_timesteps-1]-timelist[0];
	}
	else
	{
	  cout << "\nError: timegap time information requested for nonexistent timegap.\n"; 
	  exit(0);
	}

	return timegap;

}


/*-----------------------------------------------------------------------------------*/

void System::show_times(int ntimes, int * timeindices, float * times)const
{
  int timeii;

  for(timeii=0;timeii<ntimes;timeii++)
  {
    times[timeii]=timelist[timeindices[timeii]];
  }
}




/*-----------------------------------------------------------------------------------*/


int System::show_species_index(string name) const
{
  int speciesii;
  int species_index = -1;
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    if(species_name[speciesii]==name)
    {
      species_index=speciesii;
      break;
    }
  }
  return species_index;
}



/*-----------------------------------------------------------------------------------*/


int System::show_atomtype_index(string name) const
{
  int typeii;
  int atomtype_index = -1;
  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    if(atomtype_name[typeii]==name)
    {
      atomtype_index=typeii;
      break;
    }
  }

  return atomtype_index;
}


/*-----------------------------------------------------------*/
/*---------------------Time loops----------------------------*/
/*-----------------------------------------------------------*/


/*Method to loop over displacement times for a single atom defined by species, molecule_index, atom_type, and atim_index.  At each instance of the loop it passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_loop(Analysis* analysis, Trajectory * traj, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int timegapii;						//index over displacement timestep
	int block_timegapii;
	int displacement_count;
	for(timegapii=0;timegapii<n_exponential_steps;timegapii++)  //loop over exponential time step spacings within each block
	{
		displacement_count=0;
		for(blockii=0;blockii<n_exponentials;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			analysis->displacementkernel(timegapii,thisii,nextii,traj);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
//			cout << thisii << "\t" << nextii << "\n";
		
		}
	}
//	cout << "exponential\n";
	for(timegapii=n_exponential_steps; timegapii<n_timegaps-1+int(frt); timegapii++)  //loop over linear time step spacings between blocks
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				analysis->displacementkernel(timegapii,thisii,nextii,traj);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
//				cout << thisii << "\t" << nextii << "\n";
			}
			if(displacement_count == displacement_limit) break;
		}
	}
//	exit(1);

	/*Run last displacement, which is from first to last configuration only.*/
	if(!frt) {analysis->displacementkernel(n_timegaps-1,0,n_timesteps-1,traj);}
}



/*----------------------------------------------------------------------------------*/


/*Method to loop over times to produce a single displacement time for a single atom defined by species, molecule_index, atom_type, and atim_index.  It passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_loop(Analysis* analysis, Trajectory * traj, int timegapii, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int block_timegapii;
	int displacement_count;

	if(timegapii<n_exponential_steps)
	{
		displacement_count=0;
		for(blockii=0;blockii<n_exponentials;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			analysis->displacementkernel(timegapii,thisii,nextii,traj);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
		}
		
	}
	else if(timegapii>=n_exponential_steps&&timegapii<n_timegaps-1+int(frt))
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				analysis->displacementkernel(timegapii,thisii,nextii,traj);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
			}
		
			if(displacement_count == displacement_limit) break;
		}
	}
	else if(timegapii==n_timegaps-1+int(frt))
	{
	  analysis->displacementkernel(n_timegaps-1,0,n_timesteps-1,traj);
	}
	else
	{
		cout << "Error: requested timegap out of range.";
		exit(1);
	}
}


/*----------------------------------------------------------------------------------*/





/*Method to loop over displacement times for a particle list.  At each instance of the loop it passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_loop_list(Analysis* analysis, Particle_List* particle_list, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int timegapii;						//index over displacement timestep
	int block_timegapii;
	int displacement_count;

	for(timegapii=0;timegapii<n_exponential_steps;timegapii++)  //loop over exponential time step spacings within each block
	{
		displacement_count=0;
		for(blockii=0;blockii<n_exponentials;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
			analysis->list_displacementkernel(timegapii,thisii,nextii);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
//			cout << thisii << "\t" << nextii << "\n";
		
		}
	}
//	cout << "exponential\n";
	for(timegapii=n_exponential_steps; timegapii<n_timegaps-1+int(frt); timegapii++)  //loop over linear time step spacings between blocks
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
				analysis->list_displacementkernel(timegapii,thisii,nextii);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
//				cout << thisii << "\t" << nextii << "\n";
			}
		
			if(displacement_count == displacement_limit) break;
		}
	}
//	exit(1);

	/*Run last displacement, which is from first to last configuration only.*/
	//if(!frt) {analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1,particle_list);}
	if(!frt) {analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1);}
}


/*----------------------------------------------------------------------------------*/



/*Method to loop over displacement times.  At each instance of the loop it passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_list(Analysis* analysis, bool fullblock)const
{
//	int blockii;
//	int expii;

//	int timegapii;						//index over displacement timestep
//	int block_timegapii;
	//int displacement_count;
	cout << "CALLED displacement_list function (This is for testing only)!" << endl;
	if(analysis->isThreadSafe() && omp_get_max_threads() > 1)
		cout << "Analysis is thread safe, parallelizing." << endl;
	else if(analysis->isThreadSafe() && omp_get_max_threads() == 1)
		cout << "Analysis is thread safe but only one thread is permitted, running serially." << endl;
	else
		cout << "Analysis is not thread safe, running serially." << endl;

	{
		{
			int thisii;
			int nextii;
			#pragma omp parallel for schedule(dynamic) if(analysis->isThreadSafe()) // TODO: Test if we can use the old loop
			for(int timegapii=0;timegapii<n_exponential_steps;timegapii++)  //loop over exponential time step spacings within each block
			{
				int displacement_count=0;
				bool abort = false;
				for(int blockii=0;blockii<n_exponentials;blockii++)
				{
					if (!abort)
					{
						thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
						nextii = thisii+timegapii;		//calculate dispaced index
						analysis->list_displacementkernel(timegapii,thisii,nextii);
						#pragma omp atomic
						displacement_count++;
						abort = (displacement_count == displacement_limit && displacement_limit != 0);
						#pragma omp flush(abort)
					}
					//if(displacement_count == displacement_limit) break;
		//			cout << thisii << "\t" << nextii << "\n";
				}
			}
			cout << "Part 1 done!" << endl;
		}
//	cout << "exponential\n";


		
		{

			#pragma omp parallel for schedule(dynamic)  if(analysis->isThreadSafe()) // This makes this loop execute in parallel, splitting by time values.
			for(int timegapii=n_exponential_steps; timegapii<n_timegaps-1+int(frt); timegapii++)  //loop over linear time step spacings between blocks
			{
				int displacement_count=0;
				int block_timegapii = timegapii - n_exponential_steps + 1;
				bool abort = false;
				
				for(int expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
				{
					if (!abort)
					{
						for(int blockii=0; blockii<n_exponentials-block_timegapii; blockii++) // If this loop could be parallelized then you'd see even higher speed increases, but I've had a lot of data corruption when trying
						{
							if (!abort)
							{
								int thisii = n_exponential_steps*blockii+expii+int(frt);
								int nextii = thisii + n_exponential_steps*block_timegapii;
								analysis->list_displacementkernel(timegapii,thisii,nextii);
								#pragma omp atomic
								displacement_count++;
								abort = (displacement_count == displacement_limit && displacement_limit != 0);
								#pragma omp flush(abort)
//								analysis->list_displacementkernel(timegapii,thisii,nextii);
//								displacement_count++;
//						//		if (displacement_count == displacement_limit) break;
//								abort = (displacement_count == displacement_limit);
//				//				cout << thisii << "\t" << nextii << "\n";
							}
						}
					//	if (displacement_count == displacement_limit) break;
					}
				}
			}
			cout << "Part 2 done!" << endl;
		}
	}
//	exit(1);

	/*Run last displacement, which is from first to last configuration only.*/
	if(!frt) {analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1);}
}


/*----------------------------------------------------------------------------------*/

/*Method to loop over times to produce a single displacement time for a particle list.  It passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_loop_list(Analysis* analysis, Particle_List* particle_list, int timegapii, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int block_timegapii;
	int displacement_count;

	if(timegapii<n_exponential_steps)
	{
		displacement_count=0;

		for(blockii=0;blockii<n_exponentials;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
			analysis->list_displacementkernel(timegapii,thisii,nextii);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
		}
		
	}
	else if(timegapii>=n_exponential_steps&&timegapii<n_timegaps-1+int(frt))
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
				analysis->list_displacementkernel(timegapii,thisii,nextii);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
			}

		
			if(displacement_count == displacement_limit) break;
		}
	}
	else if(timegapii==n_timegaps-1+int(frt))
	{
		//analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1,particle_list);
		analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1);
	}
	else
	{
		cout << "Error: requested timegap out of range.";
		exit(1);
	}
}





/*----------------------------------------------------------------------------------*/



/*Method to loop over times to produce a single displacement time for a trajectory list.  It passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_list(Analysis* analysis, int timegapii, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int block_timegapii;
	int displacement_count;

	if(timegapii<n_exponential_steps)
	{
		displacement_count=0;
		for(blockii=0;blockii<n_exponentials;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
			analysis->list_displacementkernel(timegapii,thisii,nextii);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
		
		}
	}
	else if(timegapii>=n_exponential_steps&&timegapii<n_timegaps-1+int(frt))
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
				analysis->list_displacementkernel(timegapii,thisii,nextii);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
			}
		
			if(displacement_count == displacement_limit) break;
		}
	}
	else if(timegapii==n_timegaps-1+int(frt))
	{
		//analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1,particle_list);
		analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1);
	}
	else
	{
		cout << "Error: requested timegap out of range.";
		exit(1);
	}
}





/*----------------------------------------------------------------------------------*/


/*Method to loop over times to produce a single displacement time for a trajectory list, using only a specified range of blocks.  It passes the displacement time and two corresponding time indices to a displacementkernel method of an Analysis class.  If fullblock=1 (default), then for timegaps from one exponential block to another all times are used.  If fullblock = 0, crossblock timegaps only use the first timestep of each block.*/
void System::displacement_list(Analysis* analysis, int timegapii, int firstblock, int lastblock, bool fullblock)const
{
	int blockii;
	int expii;
	int thisii;
	int nextii;
	int block_timegapii;
	int displacement_count;

	if(firstblock<0){cout << "Error: first block cannot be less than zero";exit(1);}
	if(lastblock>=n_exponentials){cout << "Error: last block greater than number of blocks";exit(1);}

	if(timegapii<n_exponential_steps)
	{
		displacement_count=0;
		for(blockii=firstblock;blockii<=lastblock;blockii++)
		{
			thisii = n_exponential_steps*blockii+int(frt);	//calculate starting index of this block
			nextii = thisii+timegapii;		//calculate dispaced index
			//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
			analysis->list_displacementkernel(timegapii,thisii,nextii);
			displacement_count++;
			if(displacement_count == displacement_limit) break;
		
		}
	}
	else if(timegapii>=n_exponential_steps&&timegapii<n_timegaps-1+int(frt))
	{
		displacement_count=0;
		block_timegapii = timegapii - n_exponential_steps + 1;
		//if(lastblock>n_exponentials-block_timegapii){lastblock-n_exponentials-block_timegapii;}
		for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
		{
			for(blockii=firstblock; blockii<lastblock; blockii++)
			{
				thisii = n_exponential_steps*blockii+expii+int(frt);
				nextii = thisii + n_exponential_steps*block_timegapii;
				//analysis->list_displacementkernel(timegapii,thisii,nextii,particle_list);
				analysis->list_displacementkernel(timegapii,thisii,nextii);
				displacement_count++;
				if(displacement_count == displacement_limit) break;
			}
		
			if(displacement_count == displacement_limit) break;
		}
	}
	else if(timegapii==n_timegaps-1+int(frt))
	{
		//analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1,particle_list);
		analysis->list_displacementkernel(n_timegaps-1,0,n_timesteps-1);
	}
	else
	{
		cout << "Error: requested timegap out of range.";
		exit(1);
	}
}



/*method to calculate and return the number of time separation data points going into each time displacement*/
int * System::timegap_weighting(bool fullblock) const
{
  int * weighting;
  weighting = new int [n_timegaps];
  int timegapii;
  int block_timegapii;

  /*calculate weighting for timesteps within exponential blocks*/
  for(timegapii=0;timegapii<n_exponential_steps;timegapii++)
  {
    weighting[timegapii]=n_exponentials;
    if(weighting[timegapii] > displacement_limit&&displacement_limit!=0){weighting[timegapii]=displacement_limit;}  //check for limit on number of time data points per calculation and adjust weighting if appropriate
//    cout << timegapii <<"\t"<<weighting[timegapii] << "\n";
  }

  /*calculate weighting for timesteps across exponential blocks*/
  for(timegapii=n_exponential_steps;timegapii<n_timegaps-1+int(frt);timegapii++)
  {
    block_timegapii = timegapii - n_exponential_steps+1;
    weighting[timegapii]=(n_exponentials-block_timegapii)+int(fullblock)*(n_exponentials-block_timegapii)*(n_exponential_steps-1)+1;
    if(weighting[timegapii] > displacement_limit&&displacement_limit!=0){weighting[timegapii]=displacement_limit;}	//check for limit on number of time data points per calculation and adjust weighting if appropriate
//    cout << timegapii <<"\t"<< weighting[timegapii] << "\n";
  }

  /*set weighting for largest timestep*/
  if(!frt) {weighting[n_timegaps-1]=1;}
//  cout << n_timegaps-1 <<"\t"<< weighting[n_timegaps-1] << "\n";

  return weighting;
}



/*method to calculate and return the number of time separation data points going into given time displacement*/
int System::timegap_weighting(int timegap, bool fullblock) const
{
  int weighting;
  int timegapii;
  int block_timegapii;

  /*calculate weighting for timesteps within exponential blocks*/
  if(timegap>=0 && timegap < n_exponential_steps)
  {
    weighting=n_exponentials;
    if(weighting > displacement_limit&&displacement_limit!=0){weighting=displacement_limit;}  //check for limit on number of time data points per calculation and adjust weighting if appropriate
  }
  else if(timegap>=n_exponential_steps&&timegap<n_timegaps-1+int(frt))
  {
    block_timegapii = timegapii - n_exponential_steps+1;
    weighting=(n_exponentials-block_timegapii)*((int(fullblock)*(n_exponential_steps-1))+1);
    if(weighting > displacement_limit&&displacement_limit!=0){weighting=displacement_limit;}	//check for limit on number of time data points per calculation and adjust weighting if appropriate
  }
  else if(timegap=n_timegaps-1+int(frt)&&!frt)
  {
    weighting=1;
  }
  else
  {
    cout << "\nError: timegap weighting requested for nonexistent timegap.\n"; 
    exit(0);
  }
  

  return weighting;
}



/*------------------------------------------------------------*/
/*---------------------Atoms loops----------------------------*/
/*------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------*/
/*loop over centers of mass of all molecules in system*/
void System::loop_all_moleculecom(Analysis* analysis)const
{
  int speciesii;
  int moleculeii;
  Trajectory * traj;
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
    {
      print_progress(moleculeii+1,n_molecules[speciesii]);
      traj = &molecules[speciesii][moleculeii];
      analysis->atomkernel(traj);
    }
  }
}


/*loop over centers of mass of all molecules of a species*/
void System::loop_species_moleculecom(Analysis* analysis, int species_index)const
{
  int moleculeii;
  Trajectory * traj;
  for(moleculeii=0;moleculeii<n_molecules[species_index];moleculeii++)
  {
    print_progress(moleculeii+1,n_molecules[species_index]);
    traj = &molecules[species_index][moleculeii];
    analysis->atomkernel(traj);
  }
}


/*loop over all molecules of a species for an atom at a given position*/
void System::loop_atom_species(Analysis* analysis,int species_index, int atomtype, int atomindex)const
{
  int moleculeii;
  Trajectory * traj;
  cout<<"\n";cout.flush();
  for(moleculeii=0;moleculeii<n_molecules[species_index];moleculeii++)
  {
    print_progress(moleculeii+1,n_molecules[species_index]);
    traj = molecules[species_index][moleculeii].show_atom_trajectory(atomtype,atomindex);
    analysis->atomkernel(traj);
  }
}



/*-----------------------------------------------------------------------------------*/


//loop over atoms of given type within one molecule
void System::loop_type_molecule(Analysis* analysis, int species_index, int molecule_index, int atomtype)const
{
  int atomii;
  int mol_typecount = (molecules[species_index][molecule_index]).typecount(atomtype);
  Trajectory * traj;

  for(atomii=0;atomii<mol_typecount;atomii++)
  {
    traj = molecules[species_index][molecule_index].show_atom_trajectory(atomtype,atomii);
    analysis->atomkernel(traj);
  }
}



/*-----------------------------------------------------------------------------------*/


//loop over all atoms within molecule
void System::loop_molecule(Analysis* analysis, int species_index, int molecule_index)const
{
  int typeii;

  for(typeii=0;typeii<n_atomtypes;typeii++)
  {
    loop_type_molecule(analysis, species_index, molecule_index, typeii);
  }
}


/*-----------------------------------------------------------------------------------*/


//loop over all atoms within all molecules of a given species
void System::loop_species(Analysis* analysis, int species_index)const
{
  int moleculeii;
  cout<<"\n";cout.flush();
  //for(moleculeii=0;moleculeii<n_molecules[species_index];moleculeii++)
  for(moleculeii=0;moleculeii<n_molecules[species_index];moleculeii++)
  {
    loop_molecule(analysis, species_index, moleculeii);
    print_progress(moleculeii+1,n_molecules[species_index]);

  }

}




/*-----------------------------------------------------------------------------------*/



//loop all atoms of a given type within all molecules of a given species
void System::loop_type_species(Analysis* analysis,int species_index, int atomtype)const
{
  int moleculeii;
  cout<<"\n";cout.flush();
  for(moleculeii=0;moleculeii<n_molecules[species_index];moleculeii++)
  {
    loop_type_molecule(analysis, species_index, moleculeii, atomtype);
    print_progress(moleculeii+1,n_molecules[species_index]);
  }
}



/*-----------------------------------------------------------------------------------*/


//loop over all atoms of a given type within the entire system
void System::loop_type_system(Analysis* analysis, int atomtype)const
{
  int speciesii;

  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    cout << "\nSpecies " << speciesii+1 << "\n";
    loop_type_species(analysis, speciesii, atomtype);
  }
}



/*-----------------------------------------------------------------------------------*/


//loop over all atoms in the system
void System::loop_system(Analysis* analysis)const
{
  int speciesii;
  for(speciesii=0;speciesii<n_species;speciesii++)
  {
    cout << "\nSpecies " << (speciesii+1) << "\n";
    loop_species(analysis,speciesii);
  }
}



/*-----------------------------------------------------------------------------------*/


/*Method to impose periodic boundary conditions, in case lammps has let any particles slip out (it does not impose periodic boundaries every timestep)*/
void System::boxify()
{
 float xmin, ymin, zmin, xmax, ymax, zmax;
 int speciesii, moleculeii, typeii, atomii, timeii;
 Coordinate coordinate;
 Coordinate unitx(1,0,0);
 Coordinate unity(0,1,0);
 Coordinate unitz(0,0,1);
 Coordinate xshifter, yshifter, zshifter;
 int xshift, yshift, zshift;

 if (!boxified)
 {
   cout << "\nBoxifying system.";
   for(timeii=0;timeii<n_timesteps;timeii++)
   {
     xmin = box_boundary[timeii][0].show_x();
     ymin = box_boundary[timeii][0].show_y();
     zmin = box_boundary[timeii][0].show_z();
     xmax = box_boundary[timeii][1].show_x();
     ymax = box_boundary[timeii][1].show_y();
     zmax = box_boundary[timeii][1].show_z();

     xshifter=unitx*(xmax-xmin);
     yshifter=unity*(ymax-ymin);
     zshifter=unitz*(zmax-zmin);
     for(speciesii=0;speciesii<n_species;speciesii++)
     {
       for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)
       {
	 for(typeii=0;typeii<n_atomtypes;typeii++)
	 {
	   for(atomii=0;atomii<molecules[speciesii][moleculeii].typecount(typeii);atomii++)
	   {
	      xshift=yshift=zshift=0;
              coordinate = molecules[speciesii][moleculeii].show_atom_trajectory(typeii,atomii)->show_coordinate(timeii);
              if(coordinate.show_x()<xmin||coordinate.show_x()>xmax||coordinate.show_y()<ymin||coordinate.show_y()>ymax||coordinate.show_z()<zmin||coordinate.show_z()>zmax)
              {
                if(coordinate.show_x()>xmax)
                {
                  xshift = -1;
                }
                else if (coordinate.show_x()<xmin)
                {
                  xshift = 1;
                }

                if(coordinate.show_y()>ymax)
                {
                  yshift = -1;
                }
                else if (coordinate.show_y()<ymin)
                {
                  yshift = 1;
                }


                if(coordinate.show_z()>zmax)
                {
                  zshift = -1;
                }
                else if (coordinate.show_z()<zmin)
                {
                  zshift = 1;
                }
                coordinate+=(xshifter*xshift+yshifter*yshift+zshifter*zshift);
                molecules[speciesii][moleculeii].show_atom_trajectory(typeii,atomii)->set(coordinate,timeii);
              }
            }
          }
        }
      }
    }
  }
}





/*-----------------------------------------------------------------------------------*/



/*Method to pass on a single coordinate from a single atom trajectory.*/
Coordinate System::show_unwrapped(int species_index, int molecule_index, int atom_type, int atom_index, int timestep) const
{
  return molecules[species_index][molecule_index].show_unwrapped(atom_type,atom_index,timestep);
}



/*-----------------------------------------------------------------------------------*/



/*Method to write trajectory to Francis Starr's binary multifile format; each exponential block has its own file*/
void System::write_starr()const
{
  ofstream output;
  int blockii, expstepii;
  int timeii=0;
  int block_starttime;;
  int speciesii, moleculeii,typeii, atomii, timestep;
  string filename;
  int maxtime;
  stringstream tmax;
  stringstream starttime;
  double density;
  Coordinate coordinate;
  float x, y, z;

  if(!frt){cout<<"Warning: Time scheme not fully compatible with selected write method; discarding zeroth time configuration.\n";}

  maxtime = int(float(timelist[n_exponential_steps]) / float(time_unit));
  tmax << maxtime;

  for(blockii=0;blockii<n_exponentials;blockii++)	//loop over exponential blocks; each block gets its own fiel
  {
    block_starttime = float(timelist[n_exponential_steps]*blockii)/time_unit;	//calculate time just before start of block

    /*do a bunch of conversions to generate correctly formatted filename*/
    /*Format: starttimestep-timesteps_spanned.pos*/
    starttime.str("");
    starttime << (block_starttime+1);
    filename.clear();
    filename = starttime.str();
    filename += "-";
    filename += tmax.str();
    filename += ".pos";
    FILE* output;
    output = fopen(filename.c_str(),"w");

    //output.open(filename.c_str(),ios::out|ios::binary);

    for(expstepii=0;expstepii<n_exponential_steps;expstepii++)  //loop over time indices within block
    {
      timeii++;	//increment overall time index
      timestep = int(float(timelist[timeii])/time_unit);	//determine present timestep
      /*write start timestep, number of atoms, and density to binary file*/
      //output.write(reinterpret_cast<char*>(&timestep),sizeof(int));
      //output.write(reinterpret_cast<char*>(&n_atoms), sizeof(int));
      //output.write(reinterpret_cast<char*>(&density), sizeof(double));
      density = rho[timeii];
      fwrite(&timestep, sizeof(int),1,output);
      fwrite(&n_atoms, sizeof(int),1,output);
      fwrite(&density, sizeof(double),1,output);

      for(speciesii=0;speciesii<n_species;speciesii++)
      {
        for(moleculeii=0;moleculeii<n_molecules[speciesii];moleculeii++)	//loop over system molecules
        {
          for(typeii=0;typeii<n_atomtypes;typeii++)		//loop over atom types
          {
            for(atomii=0;atomii<molecules[speciesii][moleculeii].typecount(typeii);atomii++)	//loop over atoms of type
            {
              /*get coordinates for present atom at present time*/
              coordinate = molecules[speciesii][moleculeii].show_atom_trajectory(typeii,atomii)->show_coordinate(timeii);
              x = float(coordinate.show_x());
              y = float(coordinate.show_y());
              z = float(coordinate.show_z());
              /*write coordinates to binary file*/
              //output.write(reinterpret_cast<char*>(&x),sizeof(float));
              //output.write(reinterpret_cast<char*>(&y),sizeof(float));
              //output.write(reinterpret_cast<char*>(&z),sizeof(float));
              //cout << "\n" <<x<<"\t"<<y<<"\t"<<z<<"\n";
              fwrite(&x, sizeof(float),1,output);
              fwrite(&y, sizeof(float),1,output);
              fwrite(&z, sizeof(float),1,output);
              //exit(0);
            }
          }
        }
      }
    }
    //output.close();
    fclose(output);
  }

}


/*-----------------------------------------------------------------------------------*/


void System::write_single_particle(int trajii, string filename)const
{
	trajectorylist[trajii]->write(filename);
}
