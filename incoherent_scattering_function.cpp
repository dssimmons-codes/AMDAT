/*Written by David S. Simmons*/

#include <math.h>
#include <stdlib.h>
#include "incoherent_scattering_function.h"
#include "system.h"
#include "version.h"

using namespace std;


Incoherent_Scattering_Function::Incoherent_Scattering_Function()
{

  system = 0;
  wavevectors = 0;
  fullblock = 0;
  
  bin_size=0;
  
  first_bin_index = 0;
  last_bin_index = 0;
  n_spacebins = 0;
  
  n_times = 0;
  firsttime = 0;
  lasttime = -1;
  
  timetable = 0;
  timegap_weighting=0;
  
  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



Incoherent_Scattering_Function::Incoherent_Scattering_Function(const Incoherent_Scattering_Function & copy)
{
  system=copy.system;
  wavevectors = copy.wavevectors;
  bin_size = copy.bin_size;
  first_bin_index = copy.first_bin_index;
  last_bin_index = copy.last_bin_index;
  n_spacebins = copy.n_spacebins;
  fullblock = copy.fullblock;
  
  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;
  timetable = copy.timetable;
  timegap_weighting = system->timegap_weighting(fullblock);
  
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
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

Incoherent_Scattering_Function::~Incoherent_Scattering_Function()
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  delete [] correlation;
  delete [] n_atoms;
}


#endif


Incoherent_Scattering_Function Incoherent_Scattering_Function::operator =(const Incoherent_Scattering_Function& copy)
{
  if(this!=&copy)
  {
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  delete [] correlation;
  delete [] n_atoms;
  system=copy.system;
  wavevectors = copy.wavevectors;
  bin_size = copy.bin_size;
  first_bin_index = copy.first_bin_index;
  last_bin_index = copy.last_bin_index;
  n_spacebins = copy.n_spacebins;
  fullblock = copy.fullblock;
  
  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;
  timetable = copy.timetable;
  timegap_weighting = system->timegap_weighting(fullblock);
  
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
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


Incoherent_Scattering_Function::Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, bool fblock)
{
  int timeii, wavenumberii;

  system = sys;

  wavevectors = wv;
  fullblock = fblock;
  bin_size = wavevectors->show_delta_wavenumber();
  timegap_weighting = system->timegap_weighting(fullblock);
  n_atoms_represented = 0;
  
  
  first_bin_index = 0;
  last_bin_index = wavevectors->show_n_wavenumbers()-1;
  n_spacebins = last_bin_index-first_bin_index+1;

  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



Incoherent_Scattering_Function::Incoherent_Scattering_Function(System * sys, const Wave_Vectors * wv, int inner, int outer, bool fblock)
{

  int timeii, wavenumberii;

  system = sys;

  wavevectors = wv;

  bin_size = wavevectors->show_delta_wavenumber();
  fullblock = fblock;
  
  timegap_weighting = system->timegap_weighting(fullblock);
  n_atoms_represented = 0;

  first_bin_index = inner;
  last_bin_index = outer;
  n_spacebins = last_bin_index-first_bin_index+1;

  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();

   /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



/*--------------------------------------------*/
/*---Methods to use system loops over atoms---*/
/*--------------------------------------------*/


void Incoherent_Scattering_Function::preprocess()
{
}

void Incoherent_Scattering_Function::atomkernel(Trajectory * traj)
{
  system->displacement_loop(this,traj,fullblock);
  n_atoms_represented++;
}




void Incoherent_Scattering_Function::displacementkernel(int timegapii, int thisii, int nextii, Trajectory * traj)
{
  int wavenumberii;
  int wavevectorii;
  Coordinate const * vectorlist;
  int vectorcount;
  Coordinate coordinate1;
  Coordinate coordinate2;


  /*increment fourier bins*/
  for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
  {
    vectorlist = wavevectors->vectorlist(wavenumberii+first_bin_index);
    vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index);
    for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
    {
      coordinate1 = traj->show_unwrapped(thisii);
      coordinate2 = traj->show_unwrapped(nextii);
      correlation[timegapii][wavenumberii] += cos(vectorlist[wavevectorii]&(coordinate2-coordinate1)) / (float(vectorcount));
    }
  }
}




/*-----------------------------------------------------*/
/*---Methods to use trajectory list loops over atoms---*/
/*-----------------------------------------------------*/

void Incoherent_Scattering_Function::analyze(Trajectory_List * t_list)
{
	trajectory_list=t_list;

	system->displacement_list(this,fullblock);
	postprocess_list();
}

void Incoherent_Scattering_Function::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	timegap=timegapii;
//	//cout<<trajectory_list->show_n_trajectories(currenttime)<<"\t"<<timegap_weighting[timegap]<<"\t";
//	//n_atoms[timegap]+=float(trajectory_list->show_n_trajectories(currenttime))/float(timegap_weighting[timegap]);
//	//n_atoms[timegap]+=float(trajectory_list->show_n_trajectories(currenttime));
//	//cout << n_atoms[timegap]<<"\n";
//	trajectory_list->listloop(this,currenttime);
	trajectory_list->listloop(this,timegapii,thisii,nextii);

}

void Incoherent_Scattering_Function::listkernel(Trajectory* current_trajectory)
{
	int wavenumberii;
	int wavevectorii;
	Coordinate const * vectorlist;
	int vectorcount;
	Coordinate coordinate1;
	Coordinate coordinate2;
	
	double tempcorrelation=0;

	n_atoms[timegap]++;
	/*increment fourier bins*/
	for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
	{
		tempcorrelation=0;
		vectorlist = wavevectors->vectorlist(wavenumberii+first_bin_index);
		vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index);
		for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
		{
			coordinate1 = current_trajectory->show_unwrapped(currenttime);
			coordinate2 = current_trajectory->show_unwrapped(nexttime);
			tempcorrelation += double(cos(vectorlist[wavevectorii]&(coordinate2-coordinate1))) / double(vectorcount);
		}
		correlation[timegap][wavenumberii] += float(tempcorrelation);
	}	
}

void Incoherent_Scattering_Function::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
	int wavenumberii;
	int wavevectorii;
	Coordinate const * vectorlist;
	int vectorcount;
	Coordinate coordinate1;
	Coordinate coordinate2;
	
	double tempcorrelation=0;

	n_atoms[timegapii]++;
	/*increment fourier bins*/
	for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
	{
		tempcorrelation=0;
		vectorlist = wavevectors->vectorlist(wavenumberii+first_bin_index); //first_bin_index is a global variable! Watch out when threading!
		vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index); //first_bin_index is a global variable! Watch out when threading!
		for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
		{
			coordinate1 = current_trajectory->show_unwrapped(thisii);
			coordinate2 = current_trajectory->show_unwrapped(nextii);
			tempcorrelation += double(cos(vectorlist[wavevectorii]&(coordinate2-coordinate1))) / double(vectorcount);
		}
		correlation[timegapii][wavenumberii] += float(tempcorrelation);
	}	
}

/*-----------------------------------------------------*/
/*---------Methods to use trajectory list bins---------*/
/*-----------------------------------------------------*/


void Incoherent_Scattering_Function::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  
  trajectory_list=t_list;
  list_displacementkernel(timegapii, thisii, nextii);

}

void Incoherent_Scattering_Function::postprocess_bins()
{
  postprocess_list();
}
