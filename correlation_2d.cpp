/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include "correlation_2d.h"
#include "version.h"

using namespace std;


Correlation_2D::Correlation_2D()
{
  system=0;

  n_times=0;
  firsttime=0;
  lasttime=-1;

  n_spacebins=0;
  first_bin_index=0;
  last_bin_index=-1;
  bin_size=0;

  wavevectors=0;
  timegap_weighting=0;
  n_atoms_represented=0;

  n_atoms = new float [n_times];
  timetable = new float [n_times];
  correlation = new float * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    correlation[timeii]=new float [n_spacebins];
  }


}


Correlation_2D::Correlation_2D(const Correlation_2D & copy)
{
  system=copy.system;

  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;

  n_spacebins=copy.n_spacebins;
  first_bin_index=copy.first_bin_index;
  last_bin_index=copy.last_bin_index;
  bin_size=copy.bin_size;

  wavevectors=copy.wavevectors;
  timegap_weighting=copy.timegap_weighting;
  n_atoms_represented=copy.n_atoms_represented;

  n_atoms = new float [n_times];
  timetable = new float [n_times];
  correlation = new float * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = copy.n_atoms[timeii];
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=copy.correlation[timeii][wavenumberii];
    }
  }
}


#ifdef NEVER

Correlation_2D::~Correlation_2D()
{
  delete [] n_atoms;
  delete [] timetable;
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  delete [] correlation;
}

#endif

Correlation_2D Correlation_2D::operator =(const Correlation_2D & copy)
{
  if(this!=&copy)
  {
  system=copy.system;

  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;

  n_spacebins=copy.n_spacebins;
  first_bin_index=copy.first_bin_index;
  last_bin_index=copy.last_bin_index;
  bin_size=copy.bin_size;

  wavevectors=copy.wavevectors;
  timegap_weighting=copy.timegap_weighting;
  n_atoms_represented=copy.n_atoms_represented;

  n_atoms = new float [n_times];
  timetable = new float [n_times];
  correlation = new float * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = copy.n_atoms[timeii];
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=copy.correlation[timeii][wavenumberii];
    }
  }
  }

  return *this;
}



/*----------------------------------------------------------------------------------*/


/*default postprocess normalization*/
void Correlation_2D::postprocess()
{
  int timeii, binii;

  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    for(binii=0;binii<n_spacebins;binii++)
    {
            correlation[timeii][binii]/=(float(n_atoms_represented)*float(timegap_weighting[timeii]));		//normalize by number of atoms (which is a constant)
    }
  }
}



/*----------------------------------------------------------------------------------*/


/*default postprocess normalization*/
void Correlation_2D::postprocess_list()
{
	int timeii, binii;

	for(timeii=firsttime;timeii<=lasttime;timeii++)
	{
		for(binii=0;binii<n_spacebins;binii++)
		{

                correlation[timeii][binii]/=float(n_atoms[timeii]);		//normalize by average number of atoms at each timegap
                //correlation[timeii][binii]/=(n_atoms[timeii]*float(timegap_weighting[timeii]));		//normalize by average number of atoms at each timegap

		}
	}
}



/*----------------------------------------------------------------------------------*/

/*Add two correlation objects*/
Correlation_2D Correlation_2D::operator+ (const Correlation_2D & increment)const
{
  Correlation_2D temp;
  int timeii, binii;

  temp.wavevectors = wavevectors;
  if(increment.wavevectors != wavevectors)
  {
    cout << "Warning: Wave vectors for correlations do not match!\n";
  }

  if(increment.first_bin_index==first_bin_index&&increment.last_bin_index==last_bin_index&&increment.firsttime==firsttime&&increment.lasttime==lasttime)
  {
    temp.first_bin_index = first_bin_index;
    temp.last_bin_index = last_bin_index;
    temp.firsttime = firsttime;
    temp.lasttime = lasttime;
    temp.n_times = n_times;
    temp.n_spacebins = n_spacebins;
  }
  else
  {
    cout << "Error: range of correlations do not match.\n";
    exit(1);
  }

  if(increment.bin_size == bin_size)
  {
    temp.bin_size = bin_size;
  }
  else
  {
    cout << "Error: bin sizes do not match.\n";
    exit(1);
  }

   if(increment.n_atoms_represented == n_atoms_represented)
   {
     temp.n_atoms_represented = n_atoms_represented;
   }
   else
   {
     cout << "Error: number of atoms represented by correlation functions do not match.  Perhaps you intended to use '&' instead of '+'?\n";
     exit(1);
   }

   for(timeii=0;timeii<n_times;timeii++)
   {
     if(increment.n_atoms[timeii] == n_atoms[timeii])
     {
       temp.n_atoms[timeii] = n_atoms[timeii];
     }
     else
     {
       cout << "Error: number of atoms represented by correlation functions do not match.  Perhaps you intended to use '&' instead of '+'?\n";
       exit(1);
     }
   }

  for(timeii=0;timeii<n_times;timeii++)
  {
    if(timetable[timeii]==increment.timetable[timeii])
    {
      temp.timetable[timeii] = timetable[timeii];
    }
    else
    {
      cout << "Error: Timetables of correlations to be added do not match.\n";
      exit(1);
    }

    for(binii=first_bin_index;binii<=last_bin_index;binii++)
    {
      temp.correlation[timeii][binii] = correlation[timeii][binii] + increment.correlation[timeii][binii];
    }
  }
  return temp;
}



/*----------------------------------------------------------------------------------*/



/*Subtract two correlation objects*/
Correlation_2D Correlation_2D::operator- (const Correlation_2D & decrement)const
{
  Correlation_2D temp;
  int timeii, binii;

  temp.wavevectors = wavevectors;
  if(decrement.wavevectors != wavevectors)
  {
    cout << "Warning: Wave vectors for correlations do not match!\n";
  }

  if(decrement.first_bin_index==first_bin_index&&decrement.last_bin_index==last_bin_index&&decrement.firsttime==firsttime&&decrement.lasttime==lasttime)
  {
    temp.first_bin_index = first_bin_index;
    temp.last_bin_index = last_bin_index;
    temp.firsttime = firsttime;
    temp.lasttime = lasttime;
    temp.n_times = n_times;
    temp.n_spacebins = n_spacebins;
  }
  else
  {
    cout << "Error: range of correlations do not match.\n";
    exit(1);
  }

  if(decrement.bin_size == bin_size)
  {
    temp.bin_size = bin_size;
  }
  else
  {
    cout << "Error: bin sizes do not match.\n";
    exit(1);
  }

   if(decrement.n_atoms_represented == n_atoms_represented)
   {
     temp.n_atoms_represented = n_atoms_represented;
   }
   else
   {
     cout << "Error: number of atoms represented by correlation function do not match.  Perhaps you intended to use '|' instead of '-'?\n";
     exit(1);
   }

   for(timeii=0;timeii<n_times;timeii++)
   {
     if(decrement.n_atoms[timeii] == n_atoms[timeii])
     {
       temp.n_atoms[timeii] = n_atoms[timeii];
     }
     else
     {
       cout << "Error: number of atoms represented by correlation functions do not match.  Perhaps you intended to use '&' instead of '+'?\n";
       exit(1);
     }
   }

  for(timeii=0;timeii<n_times;timeii++)
  {
    if(timetable[timeii]==decrement.timetable[timeii])
    {
      temp.timetable[timeii] = timetable[timeii];
    }
    else
    {
      cout << "Error: Timetables of correlations to be added do not match.\n";
      exit(1);
    }

    for(binii=first_bin_index;binii<=last_bin_index;binii++)
    {
      temp.correlation[timeii][binii] = correlation[timeii][binii] - decrement.correlation[timeii][binii];
    }
  }
  return temp;
}



/*----------------------------------------------------------------------------------*/



Correlation_2D Correlation_2D::operator& (const Correlation_2D & increment)const
{
  Correlation_2D temp;
  int timeii, binii;

  temp.wavevectors = wavevectors;
  if(increment.wavevectors != wavevectors)
  {
    cout << "Warning: Wave vectors for correlations do not match!\n";
  }

  if(increment.first_bin_index==first_bin_index&&increment.last_bin_index==last_bin_index&&increment.firsttime==firsttime&&increment.lasttime==lasttime)
  {
    temp.first_bin_index = first_bin_index;
    temp.last_bin_index = last_bin_index;
    temp.firsttime = firsttime;
    temp.lasttime = lasttime;
    temp.n_times = n_times;
    temp.n_spacebins = n_spacebins;
  }
  else
  {
    cout << "Error: range of correlations do not match.\n";
    exit(1);
  }

  if(increment.bin_size == bin_size)
  {
    temp.bin_size = bin_size;
  }
  else
  {
    cout << "Error: bin sizes do not match.\n";
    exit(1);
  }

  temp.n_atoms_represented = increment.n_atoms_represented + n_atoms_represented;

  for(timeii=0;timeii<n_times;timeii++)
  {
    temp.n_atoms[timeii] = increment.n_atoms[timeii] + n_atoms[timeii];
    if(timetable[timeii]==increment.timetable[timeii])
    {
      temp.timetable[timeii] = timetable[timeii];
    }
    else
    {
      cout << "Error: Timetables of correlations to be added do not match.\n";
      exit(1);
    }

    for(binii=first_bin_index;binii<=last_bin_index;binii++)
    {
      temp.correlation[timeii][binii] = (correlation[timeii][binii]*n_atoms[timeii] + increment.correlation[timeii][binii]*increment.n_atoms[timeii])/temp.n_atoms[timeii];
    }
  }

  return temp;
}




/*----------------------------------------------------------------------------------*/



Correlation_2D Correlation_2D::operator| (const Correlation_2D & increment)const
{
  Correlation_2D temp;
  int timeii, binii;

  temp.wavevectors = wavevectors;
  if(increment.wavevectors != wavevectors)
  {
    cout << "Warning: Wave vectors for correlations do not match!\n";
  }

  if(increment.first_bin_index==first_bin_index&&increment.last_bin_index==last_bin_index&&increment.firsttime==firsttime&&increment.lasttime==lasttime)
  {
    temp.first_bin_index = first_bin_index;
    temp.last_bin_index = last_bin_index;
    temp.firsttime = firsttime;
    temp.lasttime = lasttime;
    temp.n_times = n_times;
    temp.n_spacebins = n_spacebins;
  }
  else
  {
    cout << "Error: range of correlations do not match.\n";
    exit(1);
  }

  if(increment.bin_size == bin_size)
  {
    temp.bin_size = bin_size;
  }
  else
  {
    cout << "Error: bin sizes do not match.\n";
    exit(1);
  }

  temp.n_atoms_represented = n_atoms_represented- increment.n_atoms_represented;

  for(timeii=0;timeii<n_times;timeii++)
  {
    temp.n_atoms[timeii] = increment.n_atoms[timeii] + n_atoms[timeii];
    if(timetable[timeii]==increment.timetable[timeii])
    {
      temp.timetable[timeii] = timetable[timeii];
    }
    else
    {
      cout << "Error: Timetables of correlations to be added do not match.\n";
      exit(1);
    }

    for(binii=first_bin_index;binii<=last_bin_index;binii++)
    {
      temp.correlation[timeii][binii] = (correlation[timeii][binii]*n_atoms[timeii] - increment.correlation[timeii][binii]*increment.n_atoms[timeii])/temp.n_atoms[timeii];
    }
  }

  return temp;
}



/*----------------------------------------------------------------------------------*/


/*Write correlation object to file*/
void Correlation_2D::write(string filename)const
{
  int timeii, binii;

  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<"."<<endl;cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Correlation data created by AMDAT v." << VERSION << "\n";
  output << "\t";
  for(binii=first_bin_index;binii<=last_bin_index;binii++)
  {
    output << wavevectors->show_approx_wavenumber(binii) << "\t";  //error here
  }
  output << "\n";

  /*Write correlation function to file*/
  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    output << timetable[timeii] << "\t"; 	//first column is time data
    for(binii=0;binii<n_spacebins;binii++)
    {
      output << correlation[timeii][binii] << "\t";
//      output << correlation[timeii][binii] << "\t"; //print unnormalized data to file
    }
    output << "\n";
  }
  output.close();
}

void Correlation_2D::write(ofstream& output)const
{
  int timeii, binii;

  cout << "\nWriting to file."<<endl;

  /*Write first row - list of bin numbers*/
  output << "Correlation data created by AMDAT v." << VERSION << "\n";
  output << "\t";
  for(binii=first_bin_index;binii<=last_bin_index;binii++)
  {
    output << wavevectors->show_approx_wavenumber(binii) << "\t";  //error here
  }
  output << "\n";

  /*Write correlation function to file*/
  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    output << timetable[timeii] << "\t"; 	//first column is time data
    for(binii=0;binii<n_spacebins;binii++)
    {
      output << correlation[timeii][binii] << "\t";
//      output << correlation[timeii][binii] << "\t"; //print unnormalized data to file
    }
    output << "\n";
  }
}
