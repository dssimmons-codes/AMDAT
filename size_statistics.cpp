/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class Size_Statistics: computes size distributions and statistics for multibodies*/
/*Written by David S. Simmons*/

#include "size_statistics.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "system.h"

using namespace std;




Size_Statistics::Size_Statistics()
{
  weighting = 0;
  moments = new float [1];
  n_moments=0;
}



Size_Statistics::Size_Statistics(const Size_Statistics & copy):Multibody_Analysis(copy)
{
  int momentii;
  
  size_count=copy.size_count;
  weighting=copy.weighting;
  n_moments=copy.n_moments;
  moments=new float [n_moments];
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    moments[momentii]=copy.moments[momentii];
  }
}

Size_Statistics Size_Statistics::operator=(const Size_Statistics & copy)
{
  int momentii;
  
  if(this!=&copy)
  {
    Multibody_Analysis::operator=(copy);
    size_count=copy.size_count;
    weighting=copy.weighting;
    n_moments=copy.n_moments;
    moments=new float [n_moments];
  
    for(momentii=0;momentii<n_moments;momentii++)
    {
      moments[momentii]=copy.moments[momentii];
    }
  }
  
  return *this;
}

Size_Statistics::~Size_Statistics()
{
  delete [] moments;
}



Size_Statistics::Size_Statistics(System*syst,int momentcount)
{
  int momentii;
  system=syst;
  n_moments=momentcount;
  weighting=0;
  moments = new float [n_moments];
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    moments[momentii]=0;
  }
}


void Size_Statistics::set(System*syst,int momentcount)
{
  int momentii;
  system=syst;
  n_moments=momentcount;
  delete [] moments;
  weighting=0;
  moments = new float [n_moments];
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    moments[momentii]=0;
  }
}

void Size_Statistics::analyze(Multibody_List * mblist)
{
  int timeii;
  multibody_list=mblist;
  weighting=0;
  
  for(timeii=0;timeii<system->show_n_timesteps();timeii++)
  {
    weighting+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0,timeii,0);
  }
  postprocess();
}



void Size_Statistics::listkernel(Multibody * multibody, int timegap, int currenttime, int nexttime)
{
  int size;
  size = multibody->show_n_bodies();
  if(size>=size_count.size())
  {
    size_count.resize(size+1,0);
  }
  size_count[size]++;
}


void Size_Statistics::postprocess()
{
  int sizeii, momentii;
  for(momentii=0;momentii<n_moments;momentii++)
  {
    moments[momentii]=0;
  }
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    moments[0]+=size_count[sizeii];
    size_count[sizeii]/=float(weighting);
    for(momentii=1;momentii<n_moments;momentii++)
    {
      moments[momentii]+=pow(float(sizeii),float(momentii))*size_count[sizeii];
    }
  }
}


void Size_Statistics::write(string filename)const
{
  int sizeii, momentii;
  
  cout << "\nWriting multibody size statistics to file.";
  ofstream output(filename.c_str());
  output << "Size statistics data created by AMDAT v." << VERSION << "\n";
  
  output << "Size\t";
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    output << sizeii << "\t";
  }
  
  output<<"\nFrequency\t";
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    output << size_count[sizeii] << "\t";
  }
  
  output<<"\n\nMoment\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<momentii<<"\t";
  }
  
  output << "\nValue\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<moments[momentii]<<"\t";
  }
  
  output << "\n\nTotal_Multibodies\t" << weighting;
  output << "\nMean_Multibodies\t" << weighting/system->show_n_timesteps();
}




void Size_Statistics::write(ofstream& output)const
{
  int sizeii, momentii;
  
  cout << "\nWriting multibody size statistics to file.";

  output << "Size statistics data created by AMDAT v." << VERSION << "\n";
  
  output << "Size\t";
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    output << sizeii << "\t";
  }
  
  output<<"\nFrequency\t";
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    output << size_count[sizeii] << "\t";
  }
  
  output<<"\n\nMoment\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<momentii<<"\t";
  }
  
  output << "\nValue\t";
  
  for(momentii=0;momentii<n_moments;momentii++)
  {
    output<<moments[momentii]<<"\t";
  }
  
  output << "\n\nTotal_Multibodies\t" << weighting;
  output << "\nMean_Multibodies\t" << weighting/system->show_n_timesteps();
  
  
  
}