/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class containing information about system composition*/
/*Written by Daniel A. Hunsicker*/


#include <stdlib.h>
#include "composition.h"
#include "version.h"
#include <omp.h>

using namespace std;

Composition::Composition()
{
    system=0;
    n_molecules=0;
    n_atomtypes=0;
    n_times=0;
    total_atoms=0;
    volume=0;
    average_density=0;
    time_average_comp = new float [n_atomtypes];
    time_scheme=-1;


    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        time_average_comp[typeii]=0;
    }

    current_comp = new float* [n_times];
    current_density = new float [n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
        current_comp[timeii] = new float [n_atomtypes];
	current_density[timeii]=0;
        for (int typeii=0; typeii<n_atomtypes;typeii++)
        {
        current_comp[timeii][typeii]=0;
        }
    }

}

Composition::Composition(System * sys, int n_xbins, int n_ybins, int n_zbins,float lx, float ly, float lz, int timescheme)
{
  int xbins,ybins,zbins;
  float x,y,z;
  time_scheme=timescheme;
  xbins = n_xbins;
  ybins = n_ybins;
  zbins = n_zbins;
  x=lx;
  y=ly;
  z=lz;

    system = sys;
    volume = (x/float(xbins))*(y/float(ybins))*(z/float(zbins));
    n_atomtypes = system->show_n_atomtypes();
    n_times = determine_n_times();
    total_atoms=0;
    average_density=0;

    time_average_comp = new float [n_atomtypes];
    current_density = new float [n_times];

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        time_average_comp[typeii]=0;
    }

    current_comp = new float* [n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
        current_comp[timeii] = new float [n_atomtypes];
	current_density[timeii]=0;

        for (int typeii=0; typeii<n_atomtypes;typeii++)
        {
        current_comp[timeii][typeii]=0;
        }
    }
}



Composition::Composition(const Composition & copy)
:Analysis_Onetime(copy)
{
    system = copy.system;
    n_atomtypes = copy.n_atomtypes;;
    n_times = copy.n_times;
    n_molecules=copy.n_molecules;
    n_times=copy.n_times;
    total_atoms = copy.total_atoms;
    volume = copy.volume;
    average_density = copy.average_density;
    time_scheme=copy.time_scheme;

    time_average_comp = new float [n_atomtypes];

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        time_average_comp[typeii]=copy.time_average_comp[typeii];
    }

    current_comp = new float* [n_times];
    current_density = new float [n_times];

    for (int timeii=0; timeii<n_times; timeii++)
    {
        current_comp[timeii] = new float [n_atomtypes];
	current_density[timeii] = copy.current_density[timeii];

        for (int typeii=0; typeii<n_atomtypes;typeii++)
        {
        current_comp[timeii][typeii]=copy.current_comp[timeii][typeii];
        }
    }
}


Composition Composition::operator=(const Composition & copy)
{
  if(this!=&copy)
  {

    delete [] time_average_comp;
    for(int timeii=0; timeii<n_times;timeii++)
    {
        delete [] current_comp[timeii];
    }
    delete [] current_comp;
    delete [] current_density;

    system = copy.system;
    n_atomtypes = copy.n_atomtypes;;
    n_times = copy.n_times;
    n_molecules=copy.n_molecules;
    n_times=copy.n_times;
    total_atoms = copy.total_atoms;
    volume = copy.volume;
    average_density = copy.average_density;
    time_scheme=copy.time_scheme;

    time_average_comp = new float [n_atomtypes];

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        time_average_comp[typeii]=copy.time_average_comp[typeii];
    }

    current_comp = new float* [n_times];
    current_density = new float [n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
        current_comp[timeii] = new float [n_atomtypes];

        for (int typeii=0; typeii<n_atomtypes;typeii++)
        {
        current_comp[timeii][typeii]=copy.current_comp[timeii][typeii];
	current_density[timeii] = copy.current_density[timeii];
        }
    }
  }
  return *this;
}


Composition::~Composition()
{
  delete [] time_average_comp;
  for (int timeii=0; timeii<n_times;timeii++)
    {
        delete [] current_comp[timeii];
    }
    delete [] current_comp;
    delete [] current_density;
}


void Composition::timekernel(int timeii)
{
  current_total_atoms=trajectory_list->show_n_trajectories(timeii);
  current_time = timeii;
  current_density[current_time] = current_total_atoms/volume;

  trajectory_list->listloop(this,current_time);

  total_atoms += current_total_atoms;
  average_density+=current_density[current_time];

  for (int typeii=0;typeii<n_atomtypes;typeii++)
  {
    time_average_comp[typeii]+=current_comp[current_time][typeii];

    current_comp[current_time][typeii]/=float(current_total_atoms);

  }
}


void Composition::postprocess_list()
{
    average_density/=n_times;
    for (int typeii=0;typeii<n_atomtypes;typeii++)
    {
	time_average_comp[typeii]/= total_atoms;
    }
}

void Composition::listkernel(Trajectory* traj)
{
    int traj_type;

    traj_type = traj->show_type()-1;
    #pragma omp atomic
    current_comp[current_time][traj_type]++;
}


void Composition::write(string filename)
{

    cout << "\nWriting composition to file " << filename << "." << endl;

    ofstream output(filename.c_str());

    output << "Composition data created by AMDAT v." << VERSION << endl << endl;

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
    output << "composition of species "<< typeii+1 <<":\t" << time_average_comp[typeii] << endl;

    }
    output << "number density:\t"<< average_density << endl;
    output << "particles:\t"<< total_atoms << endl;
    output <<endl;

//
//
//    output <<"time";
//    for (int typeii=0; typeii<n_atomtypes;typeii++)
//    {
//        output << "\tspecies "<< typeii+1;
//    }
//    output<<"\tdensity";
//    output<<endl;
//
//    for (int timeii=0;timeii<n_times;timeii++)
//    {
//        output << system->show_time(timeii) <<"\t";
//        for (int typeii=0; typeii<n_atomtypes;typeii++)
//        {
//            output << current_comp[timeii][typeii]<< "\t";
//        }
//        output << current_density[timeii];
//        output<<endl;
//    }
}

