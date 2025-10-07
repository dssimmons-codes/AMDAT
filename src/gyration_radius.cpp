// Copyright (C) 2010-2025 David S. Simmons

/*Molecular Dynamics Analysis Toolkit (AMDAT)*/
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
  max_n=0;
  rg_by_n = new float [1];
  rg_by_n[0]=0;
  weighting_by_n = new int [1];
  weighting_by_n[0]=0;
}

Gyration_Radius::Gyration_Radius(const Gyration_Radius & copy)
{
  system=copy.system;
  multibody_list=copy.multibody_list;
  gyration_radius = copy.gyration_radius;
  weighting=copy.weighting;
  max_n=copy.max_n;
  rg_by_n = new float [max_n];
  weighting_by_n = new int [max_n];
  for(int ii=0;ii<max_n;ii++)
  {
    rg_by_n[ii]=copy.rg_by_n[ii];
    weighting_by_n[ii]=copy.rg_by_n[ii];
  }
}


Gyration_Radius Gyration_Radius::operator=(const Gyration_Radius & copy)
{
  if(this!=&copy)
  {
  system=copy.system;
  multibody_list=copy.multibody_list;
  gyration_radius = copy.gyration_radius;
  weighting=copy.weighting;
  delete [] rg_by_n;
  delete [] weighting_by_n;
  max_n=copy.max_n;
  rg_by_n = new float [max_n];
  weighting_by_n = new int [max_n];
  for(int ii=0;ii<max_n;ii++)
  {
    rg_by_n[ii]=copy.rg_by_n[ii];
    weighting_by_n[ii]=copy.rg_by_n[ii];
  }
  
  }
  return *this;
}


Gyration_Radius::Gyration_Radius(System * sys)
{
  system=sys;
  multibody_list=0;
  gyration_radius = 0;
  weighting=0;
  max_n=0;
  rg_by_n = new float [1];
  rg_by_n[0]=0;
  weighting_by_n = new int [1];
  rg_by_n[0]=0;
  
}



void Gyration_Radius::analyze(Multibody_List * mblist)
{
  int timeii;
  multibody_list = mblist;
  
  delete [] rg_by_n;
  delete [] weighting_by_n;
 
  max_n = mblist->maxsize()+1;
  
  rg_by_n = new float [max_n];
  weighting_by_n = new int [max_n];
  
  for(int ii = 0; ii<max_n; ii++)
  {
    rg_by_n[ii]=0;
    weighting_by_n[ii] = 0;
  }
  
  for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    weighting+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0,timeii,0);
  }
  postprocess();
 
}




void Gyration_Radius::listkernel(Multibody * multibody, int timegap, int currenttime, int nexttime)
{
  int n_bodies;
  float rgsq;
  
  n_bodies = multibody->show_n_bodies();
  rgsq = multibody->square_gyration_radius(currenttime);
  
  weighting_by_n[n_bodies]++;
  rg_by_n[n_bodies]+=rgsq;
  
  gyration_radius+=rgsq;
  
}


void Gyration_Radius::postprocess()
{
  gyration_radius/=float(weighting);
  gyration_radius = pow(gyration_radius,0.5);
  
  for(int ii = 0; ii<max_n; ii++)
  {
    rg_by_n[ii]/=float(weighting_by_n[ii]);
    rg_by_n[ii]=pow(rg_by_n[ii],0.5);
  }
}



void Gyration_Radius::write(string filename)
{
  
  int ii;
  
  cout << "\nWriting gyration radius to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Gyration radius calculated by AMDAT v." << VERSION << "\n";
  
  output << "bodies";
  
  for(ii = 0; ii<max_n; ii++)
  {
    output << "\t" << ii;
  }
  
  output << "\nfraction";
  
  for(ii = 0; ii<max_n; ii++)
  {
    output << "\t" << float(weighting_by_n[ii])/float(weighting);
  }
  
  output<<"\nRg";
  
  for(ii = 0; ii<max_n; ii++)
  {
    output << "\t" << rg_by_n[ii];
  }
  
 
  
  output << "\n\noverall mean\t" << gyration_radius;
}

void Gyration_Radius::write(ofstream& output)const
{
  cout << "\nWriting gyration radius to file.";

  output << "Gyration radius calculated by AMDAT v." << VERSION << "\n";
  
  output << gyration_radius;
}
