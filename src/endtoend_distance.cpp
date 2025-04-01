/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to calculate multibody mean end-to-end distance*/
/*Written by David S. Simmons*/

#include <sstream>
#include <iostream>
#include <string>
#include <math.h>

#include "endtoend_distance.h"
#include "version.h"
#include "system.h"

using namespace std;



Endtoend_Distance::Endtoend_Distance()
{
  system=0;
  multibody_list=0;
  endtoend_distance = 0;
  weighting=0;
  max_n=0;
  ree_by_n = new float [1];
  ree_by_n[0]=0;
  weighting_by_n = new int [1];
  weighting_by_n[0]=0;
}

Endtoend_Distance::Endtoend_Distance(const Endtoend_Distance & copy)
{
  system=copy.system;
  multibody_list=copy.multibody_list;
  endtoend_distance = copy.endtoend_distance;
  weighting=copy.weighting;
  max_n=copy.max_n;
  ree_by_n = new float [max_n];
  weighting_by_n = new int [max_n];
  for(int ii=0;ii<max_n;ii++)
  {
    ree_by_n[ii]=copy.ree_by_n[ii];
    weighting_by_n[ii]=copy.ree_by_n[ii];
  }
}


Endtoend_Distance Endtoend_Distance::operator=(const Endtoend_Distance & copy)
{
  if(this!=&copy)
  {
  system=copy.system;
  multibody_list=copy.multibody_list;
  endtoend_distance = copy.endtoend_distance;
  weighting=copy.weighting;
  delete [] ree_by_n;
  delete [] weighting_by_n;
  max_n=copy.max_n;
  ree_by_n = new float [max_n];
  weighting_by_n = new int [max_n];
  for(int ii=0;ii<max_n;ii++)
  {
    ree_by_n[ii]=copy.ree_by_n[ii];
    weighting_by_n[ii]=copy.ree_by_n[ii];
  }

  }
  return *this;
}


Endtoend_Distance::Endtoend_Distance(System * sys)
{
  system=sys;
  multibody_list=0;
  endtoend_distance = 0;
  weighting=0;
  max_n=0;
  ree_by_n = new float [1];
  ree_by_n[0]=0;
  weighting_by_n = new int [1];
  ree_by_n[0]=0;

}



void Endtoend_Distance::analyze(Multibody_List * mblist)
{
  int timeii;
  multibody_list = mblist;

  delete [] ree_by_n;
  delete [] weighting_by_n;

  max_n = mblist->maxsize()+1;

  ree_by_n = new float [max_n];
  weighting_by_n = new int [max_n];

  for(int ii = 0; ii<max_n; ii++)
  {
    ree_by_n[ii]=0;
    weighting_by_n[ii] = 0;
  }

  for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    weighting+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0,timeii,0);
  }
  postprocess();

}

void Endtoend_Distance::listkernel(Multibody * multibody, int timegap, int currenttime, int nexttime)
{
  int n_bodies;
  float reesq;

  n_bodies = multibody->show_n_bodies();
  reesq = multibody->square_endtoend_distance(currenttime);

  weighting_by_n[n_bodies]++;
  ree_by_n[n_bodies]+=reesq;

  endtoend_distance+=reesq;

}

void Endtoend_Distance::postprocess()
{
  endtoend_distance/=float(weighting);
  endtoend_distance = pow(endtoend_distance,0.5);

  for(int ii = 0; ii<max_n; ii++)
  {
    ree_by_n[ii]/=float(weighting_by_n[ii]);
    ree_by_n[ii]=pow(ree_by_n[ii],0.5);
  }
}

void Endtoend_Distance::write(string filename)
{

  int ii;

  cout << "\nWriting end-to-end distance to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "end-to-end distance calculated by AMDAT v." << VERSION << "\n";

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
    output << "\t" << ree_by_n[ii];
  }

  output << "\n\noverall mean\t" << endtoend_distance;
}

void Endtoend_Distance::write(ofstream& output)const
{
  cout << "\nWriting end-to-end distance to file.";

  output << "end-to-end distance calculated by AMDAT v." << VERSION << "\n";

  output << endtoend_distance;
}
