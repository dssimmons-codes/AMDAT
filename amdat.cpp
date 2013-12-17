/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*A software package for analysis of molecular dynamics trajectories*/
/*Written by David S. Simmons*/
/*Major contributors: Michael Marvin, Mark Mackura*/
/*Other contributors: Daniel Hunsicker, Ryan Lang*/


/*'Main'*/

#include "control.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "version.h"

#include <omp.h>

using namespace std;

#include "progress.h"
#define LISTSIZE 100

bool fileoutput;

int main(int argc, char *argv[])
{
  string * constants;
  constants = new string [LISTSIZE];
  string * constant_names;
  constant_names = new string [LISTSIZE];
  char * inputfile;
  bool input, procs_given;
  string t_file = "NULL";
  input=0;
  int n_constants;
  n_constants=0;
  for(int argii=0;argii<argc;argii++) // loop over tokens
  {
    if(argv[argii][0] == '-') // flags
    {
      if(argv[argii][1] == 'v')
      {
        cout << "\nAmorphous Molecular Dynamics Analysis Toolkit (AMDAT) v."<<VERSION;
        cout << "\nRelease date: "<<DATE<<endl<<endl;
        exit(0);
      }
      else if(argv[argii][1] == 'c') // constant flag
      {
        constant_names[n_constants] = argv[argii+1];
        constants[n_constants] = argv[argii+2];
        n_constants++;
      }
      else if(argv[argii][1] == 't') // temporary input file flag
      {
          t_file = argv[argii+1];
      }
      else if(argv[argii][1] == 'n') // number of processor cores to use flag
      {
          procs_given=true;
          omp_set_num_threads(atoi(argv[argii+1]));
      }
      else if(argv[argii][1] == 'h') // help flag
      {
        cout<<"\nAmorphous Molecular Dynamics Analysis Toolkit (AMDAT) v."<<VERSION;
        cout<<"\nRelease date: "<<DATE;
        cout<<"\n\nAMDAT is a toolkit of analysis methods appropriate for molecular dynamics simulations of amorphous matter.\nMore detailed descriptions of syntax and analysis methods can be found in the included documentation file.\n\n";
        cout<< "Examples:\n\n";
        cout<<"./AMDAT -n 8 -i input.amdat\t\t\t\t\t# Runs AMDAT with 8 processing cores available for multithreaded analysis methods.\n";
        cout<<"./AMDAT -c CONSTANT 33 -c TEMP 1.00 -c QVECTOR_PATH /data/qvector -i input.amdat\t# Replaces constants in amdat.input with those defined, and then runs program on resulting commands.\n";
        cout<<"./AMDAT -i input.amdat\t\t\t\t\t\t\t\t\t# Runs AMDAT based on commands from input.amdat\n";
        cout<<"./AMDAT -i input.amdat -w 10\t\t\t\t\t\t\t\t# Waits for 10 seconds, then runs AMDAT based on commands from input.amdat\n";
        cout<<"./AMDAT -i input.amdat > amdat.log 2>&1\t\t\t\t\t\t\t# Runs AMDAT based on commands from input.amdat with console output redirected to amdat.log\n";
        cout<<"./AMDAT -v \t\t\t\t\t\t\t\t\t\t# Displays version information and exits program.\n\n";
        cout<<"Flags:\n\n";
        cout<<"-i \t\t Input file. REQUIRED to run progam"<<endl;
        cout<<"-v \t\t Display version information and exits"<<endl;
        cout<<"-c \t\t Define constant"<<endl;
        cout<<"-n \t\t Set the number of available processing cores"<<endl;
        cout<<"-w \t\t Makes the system wait before reading in"<<endl;
        cout<<"-h, --help \t Display this message"<<endl<<endl;

        exit(0);
      }
      else if(argv[argii][1] == 'i')
      {
        inputfile=argv[argii+1];
        input=1;
      }
      else if(argv[argii][1] == 'w')
      {
        int time;
        time=atoi(argv[argii+1]);
        sleep(time);
      }
    }
    if(argv[argii][0] == '-') // flags
    {
      if(argv[argii][1] == '-') // flags
      {
        if(argv[argii][2] == 'h') // help flag
        {
        cout<<"\nAmorphous Molecular Dynamics Analysis Toolkit (AMDAT) v."<<VERSION;
        cout<<"\nRelease date: "<<DATE;
        cout<<"\n\nAMDAT is a toolkit of analysis methods appropriate for molecular dynamics simulations of amorphous matter.\nMore detailed descriptions of syntax and analysis methods can be found in the included documentation file.\n\n";
        cout<< "Examples:\n\n";
        cout<<"./AMDAT -n 8 -i input.amdat\t\t\t\t\t# Runs AMDAT with 8 processing cores available for multithreaded analysis methods.\n";
        cout<<"./AMDAT -c CONSTANT 33 -c TEMP 1.00 -c QVECTOR_PATH /data/qvector -i input.amdat\t# Replaces constants in amdat.input with those defined, and then runs program on resulting commands.\n";
        cout<<"./AMDAT -i input.amdat\t\t\t\t\t\t\t\t\t# Runs AMDAT based on commands from input.amdat\n";
        cout<<"./AMDAT -i input.amdat -w 10\t\t\t\t\t\t\t\t# Waits for 10 seconds, then runs AMDAT based on commands from input.amdat\n";
        cout<<"./AMDAT -i input.amdat > amdat.log 2>&1\t\t\t\t\t\t\t# Runs AMDAT based on commands from input.amdat with console output redirected to amdat.log\n";
        cout<<"./AMDAT -v \t\t\t\t\t\t\t\t\t\t# Displays version information and exits program.\n\n";
        cout<<"Flags:\n\n";
        cout<<"-i \t\t Input file. REQUIRED to run progam"<<endl;
        cout<<"-v \t\t Display version information and exits"<<endl;
        cout<<"-c \t\t Define constant"<<endl;
        cout<<"-n \t\t Set the number of available processing cores"<<endl;
        cout<<"-w \t\t Makes the system wait before reading in"<<endl;
        cout<<"-h, --help \t Display this message"<<endl<<endl;

        exit(0);
        }
      }
    }
  }

  if(!input)
  {
    cout << "No input file specified.  Please include input file name in command line.\n";
    exit(1);
  }
  if(!procs_given)
    omp_set_num_threads(1);

  cout<<"Amorphous Molecular Dynamics Analysis Toolkit (AMDAT) v."<<VERSION<<endl;
  Control control(inputfile,constants,constant_names,n_constants,t_file);
  return 0;
}


void print_progress(int done,int total)
{
	float percent;
	if(!fileoutput)
	{
	  percent = int(100*float(done)/float(total));
	  cout.width(3);
	  cout << "\b\b\b\b" << percent << "%";
	  cout.flush();
	}
}
