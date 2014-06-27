/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate multibody mean gyration radius*/
/*Written by David S. Simmons*/

#include <sstream>
#include <iostream>
#include <string>
#include <math.h>

#include "gyration_radius.h"
#include "version.h"
#include "system.h"

using namespace std;



Gyration_Radius::Gyration_Radius()
{
  system=0;
  multibody_list=0;
  gyration_radius = 0;
  weighting=0;
}

Gyration_Radius::Gyration_Radius(const Gyration_Radius & copy)
{
  system=copy.system;
  multibody_list=copy.multibody_list;
  gyration_radius = copy.gyration_radius;
  weighting=copy.weighting;
}


Gyration_Radius Gyration_Radius::operator=(const Gyration_Radius & copy)
{
  if(this!=&copy)
  {
  system=copy.system;
  multibody_list=copy.multibody_list;
  gyration_radius = copy.gyration_radius;
  weighting=copy.weighting;
  }
  return *this;
}


Gyration_Radius::Gyration_Radius(System * sys)
{
  system=sys;
  multibody_list=0;
  gyration_radius = 0;
  weighting=0;
}



void Gyration_Radius::analyze(Multibody_List * mblist)
{
  int timeii;
  multibody_list = mblist;
  
  for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    weighting+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0,timeii,0);
  }
  postprocess();
}




void Gyration_Radius::listkernel(Multibody * multibody, int timegap, int currenttime, int nexttime)
{
  gyration_radius+=multibody->square_gyration_radius(currenttime);
}


void Gyration_Radius::postprocess()
{
  gyration_radius/=float(weighting);
  gyration_radius = pow(gyration_radius,0.5);
}



void Gyration_Radius::write(string filename)
{
  cout << "\nWriting gyration radius to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Gyration radius calculated by AMDAT v." << VERSION << "\n";
  
  output << gyration_radius;
}
