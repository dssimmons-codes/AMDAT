/*Value_List class - stores a value for each trajectory in the system*/
/*Amorphous Molecular Dynamics Analysis Toolkit (aMDAT)*/
/*Written by David S. Simmons*/

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>

#include "boolean_list.h"
#include "system.h"
#include "trajectory_list.h"

#ifndef VALUE_LIST
#define VALUE_LIST

namespace std
{

template <class valType>
class Value_List
{
  
protected:
  System * syst;

  vector<vector<valType>> values;	//array of values stored: [internal_time][trajectory_ID]
  //valType ** values;			//array of values stored: [internal_time][trajectory_ID]
  mutable Boolean_List * included;		//array of boolean lists specifying which trajectories are in value list at each time: [internal_time]

  //data on time scheme
  int n_times;				//number of unique times stored
  int * time_conversion;	//time conversion table converting system time index to internal storage time index
  bool * defined_times;		//stores a 1 at times for which the value_list was originall calculated

  bool threshold(int, int, bool, valType)const;
  bool threshold(int, int, valType, valType)const;
  
  void update_size();

  public:

  Value_List();
  Value_List(System*);
  Value_List(const Value_List<valType>&);
  Value_List operator = (const Value_List &);

  void set(System*);
  void set(int,int,valType);

  valType show_value(int timeii, int trajii){return values[convert_time(timeii)][trajii];};

  void operator () (int timeii, int trajii){return values[convert_time(timeii)][trajii];};

  //bool is_included(int trajii){return included(int trajii);};

  int n_included()const;
  float mean()const;
  float power_mean(float power)const;
  float variance()const;
  valType max()const;	
  valType min()const;
  void distribution(string filename, float max, int n_bins)const;	//new - code this

  Value_List operator * (float)const;						
  void operator *= (float);	
  //Value_List operator * (const Value_List&)const;
  //void operator *= (const Value_List&);				
  Value_List power(float);	
  float static_crosscorrelation(const Value_List&)const;
  void dynamic_crosscorrelation_function(string filename, const Value_List&)const;
  void dynamic_autocorrelation_function(string filename)const;
  

  int convert_time(int timeii)const{return time_conversion[timeii];};	//convert requested time (Where the index is the time index from the system object) to internal time index

  string write_pdb(int, string, int, string*)const;
  string write_pdb(int, string, int)const;
  string write_pdb(int, string, int, string*,valType)const;
  string write_pdb(int, string, int,valType)const;
  void write_data_file(string filename)const;
  
  void construct_t_list(bool, valType,Trajectory_List*);
  void construct_t_list(valType, valType,Trajectory_List*);
  void percentile_t_list(bool, float,Trajectory_List*);
  void percentile_t_list(float, float,Trajectory_List*);

  virtual void set_time_conv(){};
  
 // Value_List operator *(Trajectory_List)const;
};


template <class valType>
Value_List<valType>::Value_List()
{
  syst = 0;
  n_times = 0;
  
  included = new Boolean_List[n_times];
  time_conversion = new int [0];
  defined_times = new bool [0];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
  }
}

template <class valType>
Value_List<valType>::Value_List(System * sys)
{
  syst=sys;
  n_times = syst->show_n_timesteps();
  
  included = new Boolean_List[n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
  }
  
  time_conversion = new int [syst->show_n_timesteps()];
  defined_times = new bool [syst->show_n_timesteps()];
  for(int timeii=0;syst->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=timeii;
    defined_times[timeii]=1;
  }
  
}

template <class valType>
Value_List<valType>::Value_List(const Value_List<valType>& copy)
{
  int timeii;
  
  syst=copy.syst;
  n_times = copy.n_times;

  included = new Boolean_List [n_times];
  time_conversion = new int [syst->show_n_timesteps()];
  defined_times = new bool [syst->show_n_timesteps()];
  values=copy.values;
  for(timeii=0;timeii<n_times;timeii++)
  {
	  included[timeii] = copy.included[timeii];
  }
  
  for(timeii=0;timeii<syst->show_n_timesteps();timeii++)
  {
    time_conversion[timeii] = copy.time_conversion[timeii];
    defined_times[timeii] = copy.defined_times[timeii];
  }

}

template <class valType>
Value_List<valType> Value_List<valType>::operator = (const Value_List<valType> & copy)
{
  int timeii;

  if(this!=&copy)
  {
    syst=copy.syst;
    n_times = copy.n_times;

    included = new Boolean_List [n_times];
    time_conversion = new int [syst->show_n_timesteps()];
    defined_times = new bool [syst->show_n_timesteps()];
    values=copy.values;
    for(timeii=0;timeii<n_times;timeii++)
    {
      included[timeii] = copy.included[timeii];
    }
  
    for(timeii=0;timeii<syst->show_n_timesteps();timeii++)
    {
      time_conversion[timeii] = copy.time_conversion[timeii];
      defined_times[timeii] = copy.defined_times[timeii];
    }

  }

  return *this;

}



template <class valType>
void Value_List<valType>::set(System * sys)
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
	  delete [] values[timeii];
  }
  delete [] values;
  delete [] included;
  delete [] time_conversion;
  delete [] defined_times;

  syst=sys;
  n_times = syst->show_n_timesteps();
  
  included = new Boolean_List[n_times];
  values = new valType * [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
  }
  
  time_conversion = new int [syst->show_n_timesteps()];
  defined_times = new bool [syst->show_n_timesteps()];
  for(int timeii=0;syst->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=timeii;
    defined_times[timeii]=1;
  }
}



template <class valType>
void Value_List<valType>::set(int timeii, int trajii, valType val)
{
  included[timeii](trajii,1);
  values[timeii][trajii] = val;
}




template <class valType>
float Value_List<valType>::mean()const
{
  float avg=0;
  int weighting=0;
  //sum over values of trajectories that are included

  for (int timeii=0; timeii<n_times; timeii++)
  {
	for (int trajii=0;trajii<values[timeii].size();trajii++)
	{
	 avg += float(values[timeii][trajii])*float(included[timeii](trajii));
	 weighting+=int(included[timeii](trajii));
	}
  }

  //normalize
  avg/=float(weighting);

  return avg;
}



template <class valType>
float Value_List<valType>::power_mean(float power)const
{
  float avg=0;
  int weighting=0;
  //sum over values of trajectories that are included

  for (int timeii=0; timeii<n_times; timeii++)
  {
	for (int trajii=0;trajii<values[timeii].size();trajii++)
	{
	 avg += pow(float(values[timeii][trajii]),power)*float(included[timeii](trajii));
	 weighting+=int(included[timeii](trajii));
	}
  }

  //normalize
  avg/=float(weighting);

  return avg;
}


template <class valType>
float Value_List<valType>::variance()const
{
  float avg=mean();
  int weighting=0;
  float var=0;
  //sum over values of trajectories that are included

  for (int timeii=0; timeii<n_times; timeii++)
  {
	for (int trajii=0;trajii<values[timeii].size();trajii++)
	{
	 var += pow(float(values[timeii][trajii])*float(included[timeii](trajii))-avg,2.0);
	 weighting+=int(included[timeii](trajii));
	}
  }

  //normalize
  var/=float(weighting);

  return var;
}


template <class valType>
valType Value_List<valType>::max()const
{
  valType maximum;

  maximum = values[included[0].first_included()];
  for (int timeii=0; timeii<n_times; timeii++)
  {
	for (int trajii=included[timeii].first_included()+1;trajii<values[timeii].size();trajii++)
	{
		if(included[timeii](trajii)){maximum=(values[timeii][trajii]>maximum?values[timeii][trajii]:maximum);}
	}
  }
  return maximum;
}

template <class valType>
valType Value_List<valType>::min()const
{
  valType minimum;

  minimum = values[included[0].first_included()];
  for (int timeii=0; timeii<n_times; timeii++)
  {
	for (int trajii=included[timeii].first_included()+1;trajii<values[timeii].size();trajii++)
	{
		if(included[timeii](trajii)){minimum=(values[timeii][trajii]>minimum?values[timeii][trajii]:minimum);}
	}
  }
  return minimum;
}

template <class valType>
void Value_List<valType>::distribution(string filename, float max, int n_bins)const
{
	float * dist;
	int timeii, trajii, binii;
	float binsize = max/float(n_bins);
	int bin;
	int weighting = 0;

	dist = new float [n_bins+1];	//last bin is overflow
	for(binii=0;binii<=n_bins;binii++);
	{
		dist[binii]=0;		//zero out bin values
	}

	for(timeii=0;timeii<n_times;timeii++)
	{
		for(trajii=0;trajii<values[timeii].size();trajii++)
		{
			bin = int(float(values[timeii][trajii])/binsize)-1;
			if(bin>n_bins)bin=n_bins;
			dist[binii]++;
			weighting+=int(included[timeii](trajii));
		}	
	}

	for(binii=0;binii<=n_bins;binii++)
	{
		dist[binii]/=float(weighting);
	}
}

template <class valType>
Value_List<valType> Value_List<valType>::operator * (float multiplier)const
{
	Value_List<valType> templist(*this);
	
	for(int timeii=0;timeii<n_times;timeii++)
	{
		for(int trajii=0;trajii<values[timeii].size();trajii++)
		{
			if(included[timeii](trajii))
			{
				templist.values[timeii][trajii]*=multiplier;
			}
		}
	}

	return templist;
}

#ifdef NEVER

template <class valType>
Value_List<valType> Value_List<valType>::operator * (const Value_List & multiplier)const
{
//  int overlap_times = -1;
  if(syst!=multiplier.syst)
  {
    cout << "Error: Value_Lists to be multiplied have different systems.\n";
    exit(0);
  }
  
  Value_List<valType> templist(*this);
//  for(timeii=0;timeii<syst->show_n_timesteps();timeii++)
//  {
//    templist.defined_times[timeii]=defined_times[timeii]&&multiplier.defined_times[timeii];
//    overlap_times+=int(defined_times[timeii]&&multiplier.defined_times[timeii])
//  }
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
    for(int trajii=0;trajii<values[timeii].size();trajii++)
    {
      if(included[timeii](trajii))
      {
	templist.values[timeii][trajii]*=valType(multiplier.values[timeii][trajii]);
      }
    }
  }
  
  return templist;
}

template <class valType>
void Value_List<valType>::operator *= (const Value_List & multiplier)
{

	if(syst!=multiplier.syst)
	{
	  cout << "Error: Value_Lists to be multiplied have different systems.\n";
	  exit(0);
	}
  
	for(int timeii=0;timeii<n_times;timeii++)
	{
		for(int trajii=0;trajii<values[timeii].size();trajii++)
		{
			if(included[timeii](trajii))
			{
				values[timeii][trajii]*=valType(multiplier.values[timeii][trajii]);
			}
		}
	}
}

#endif

template <class valType>
void Value_List<valType>::operator *= (float multiplier)
{
	for(int timeii=0;timeii<n_times;timeii++)
	{
		for(int trajii=0;trajii<values[timeii].size();trajii++)
		{
			if(included[timeii](trajii))
			{
				values[timeii][trajii]*=multiplier;
			}
		}
	}
}


template <class valType>
Value_List<valType> Value_List<valType>::power(float power)
{
	Value_List<valType> templist(*this);

	for(int timeii=0;timeii<n_times;timeii++)
	{
		for(int trajii=0;trajii<values[timeii].size();trajii++)
		{
			if(included[timeii](trajii))
			{
				templist.values[timeii][trajii]=valType(pow(float(templist.values[timeii][trajii]),power));
			}
		}
	}
	return templist;
}


template <class valType>
float Value_List<valType>::static_crosscorrelation(const Value_List& target)const
{
  float correlation=0;
  int n_atoms=0;
  
  if(syst!=target.syst)
  {
    cout << "Error: Value_Lists to be multiplied have different systems.\n";
    exit(0);
  }
  
  int n_trajectories = syst->show_n_trajectories();
  ((Value_List<valType>)(*this)).update_size();
  ((Value_List<valType>)(target)).update_size();
  
  for(int timeii=0;timeii<syst->show_n_timesteps();timeii++)
  {
    for(int trajii=0;trajii<n_trajectories;trajii++)
    {
      if(included[convert_time(timeii)](trajii)&&target.included[target.convert_time(timeii)](trajii))
      {
	correlation+=float(values[convert_time(timeii)][trajii])*float(target.values[target.convert_time(timeii)][trajii]);
	n_atoms++;
      }
    }
  }
  
  correlation/=(float(n_atoms)*mean()*target.mean());
  
  return correlation;
}

template <class valType>
void Value_List<valType>::dynamic_crosscorrelation_function(string filename, const Value_List& target)const
{
  int blockii, expii, thisii, nextii;
  int timegapii;						//index over displacement timestep
  int block_timegapii;
  int displacement_count;
  bool fullblock=0;
  float * timetable; 
  ofstream output;
  
  if(syst!=target.syst)
  {
    cout << "Error: Value_Lists to be multiplied have different systems.\n";
    exit(0);
  }
  
  
  int displacement_limit = syst->show_displacement_limit();
  int n_timegaps = syst->show_n_timegaps();
  int n_exponential_steps = syst->show_n_exponential_steps();
  bool frt = syst->show_frt();
  int n_exponentials = syst->show_n_exponentials();
  
  int n_trajectories = syst->show_n_trajectories();
  ((Value_List<valType>)(*this)).update_size();
  ((Value_List<valType>)(target)).update_size();
  
  float * correlation;
  int * n_atoms;
  correlation = new float [n_timegaps];
  n_atoms = new int [n_timegaps];
  for(timegapii=0;timegapii<n_timegaps;timegapii++)
  {
    correlation[timegapii]=0;
    n_atoms[timegapii]=0;
  }
  
  for(timegapii=0;timegapii<n_exponential_steps;timegapii++)  //loop over exponential time step spacings within each block
  {
    displacement_count=0;
    for(blockii=0;blockii<n_exponentials;blockii++)
    {
      thisii = convert_time(n_exponential_steps*blockii+int(frt));	//calculate starting index of this block
      nextii = target.convert_time(n_exponential_steps*blockii+int(frt)+timegapii);		//calculate dispaced index
      
      for(int trajii=0;trajii<n_trajectories;trajii++)
      {
	if(included[thisii](trajii)&&target.included[nextii](trajii))
        {
	  correlation[timegapii]+=float(values[thisii][trajii])*float(target.values[nextii][trajii]);
	  n_atoms[timegapii]++;
	}
      }
      displacement_count++;
      if(displacement_count == displacement_limit) break;
    }
  }
  for(timegapii=n_exponential_steps; timegapii<n_timegaps-1+int(frt); timegapii++)  //loop over linear time step spacings between blocks
  {
    displacement_count=0;
    block_timegapii = timegapii - n_exponential_steps + 1;
    for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
    {
      for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
      {
	thisii = convert_time(n_exponential_steps*blockii+int(frt));	//calculate starting index of this block
	nextii = target.convert_time(n_exponential_steps*blockii+int(frt)+timegapii);		//calculate dispaced index
	for(int trajii=0;trajii<n_trajectories;trajii++)
	{
	  if(included[thisii](trajii)&&target.included[nextii](trajii))
	  {
	    correlation[timegapii]+=float(values[thisii][trajii])*float(target.values[nextii][trajii]);
	    n_atoms[timegapii]++;
	  }
	}
	displacement_count++;
	if(displacement_count == displacement_limit) break;
      }
      if(displacement_count == displacement_limit) break;
    }
  }
  /*Run last displacement, which is from first to last configuration only.*/
  if(!frt)
  {
    thisii = convert_time(0);	//calculate starting index of this block
    nextii = target.convert_time(syst->show_n_timesteps()-1);		//calculate dispaced index
    for(int trajii=0;trajii<n_trajectories;trajii++)
    {
      if(included[thisii](trajii)&&target.included[nextii](trajii))
      {
        correlation[n_timegaps-1]+=float(values[thisii][trajii])*float(target.values[nextii][trajii]);
        n_atoms[timegapii]++;
      }
    }
  }
  
  timetable = syst->displacement_times();
  
  output.open(filename.c_str(),ios::out | ios::app);
  for(timegapii=0;timegapii<n_timegaps;timegapii++)
  {
    correlation[timegapii]/=(float(n_atoms[timegapii])*mean()*target.mean());
    output << timetable[timegapii]<< "\t" << correlation[timegapii]<<"\n";
  }
  output.close();
}




template <class valType>
void Value_List<valType>::dynamic_autocorrelation_function(string filename)const
{
int blockii, expii, thisii, nextii;
  int timegapii;						//index over displacement timestep
  int block_timegapii;
  int displacement_count;
  bool fullblock=0;
  float * timetable; 
  ofstream output;
  
  int displacement_limit = syst->show_displacement_limit();
  int n_timegaps = syst->show_n_timegaps();
  int n_exponential_steps = syst->show_n_exponential_steps();
  bool frt = syst->show_frt();
  int n_exponentials = syst->show_n_exponentials();
  
  ((Value_List<valType>)(*this)).update_size();
  int n_trajectories = syst->show_n_trajectories();
  
  float * correlation;
  int * n_atoms;
  correlation = new float [n_timegaps];
  n_atoms = new int [n_timegaps];
  for(timegapii=0;timegapii<n_timegaps;timegapii++)
  {
    correlation[timegapii]=0;
    n_atoms[timegapii]=0;
  }
  
  for(timegapii=0;timegapii<n_exponential_steps;timegapii++)  //loop over exponential time step spacings within each block
  {
    displacement_count=0;
    for(blockii=0;blockii<n_exponentials;blockii++)
    {
      thisii = convert_time(n_exponential_steps*blockii+int(frt));	//calculate starting index of this block
      nextii = convert_time(n_exponential_steps*blockii+int(frt)+timegapii);		//calculate dispaced index
      for(int trajii=0;trajii<n_trajectories;trajii++)
      {
	if(included[thisii](trajii)&&included[nextii](trajii))
	{
	  correlation[timegapii]+=float(values[thisii][trajii])*float(values[nextii][trajii]);
	  n_atoms[timegapii]++;
	}
      }
      displacement_count++;
      if(displacement_count == displacement_limit) break;
    }
  }
  for(timegapii=n_exponential_steps; timegapii<n_timegaps-1+int(frt); timegapii++)  //loop over linear time step spacings between blocks
  {
    displacement_count=0;
    block_timegapii = timegapii - n_exponential_steps + 1;
    for(expii=0;expii<((int(fullblock)*(n_exponential_steps-1))+1);expii++)
    {
      for(blockii=0; blockii<n_exponentials-block_timegapii; blockii++)
      {
	thisii = convert_time(n_exponential_steps*blockii+int(frt));	//calculate starting index of this block
	nextii = convert_time(n_exponential_steps*blockii+int(frt)+timegapii);		//calculate dispaced index
	for(int trajii=0;trajii<n_trajectories;trajii++)
	{
	  if(included[thisii](trajii)&&included[nextii](trajii))
	  {
	    correlation[timegapii]+=float(values[thisii][trajii])*float(values[nextii][trajii]);
	    n_atoms[timegapii]++;
	  }
	}    
	displacement_count++;
	if(displacement_count == displacement_limit) break;
      }
      if(displacement_count == displacement_limit) break;
    }
  }
  /*Run last displacement, which is from first to last configuration only.*/
  if(!frt)
  {
    thisii = convert_time(0);	//calculate starting index of this block
    nextii = convert_time(syst->show_n_timesteps()-1);		//calculate dispaced index
    for(int trajii=0;trajii<n_trajectories;trajii++)
    {
      if(included[thisii](trajii)&&included[nextii](trajii))
      {
	correlation[n_timegaps-1]+=float(values[thisii][trajii])*float(values[nextii][trajii]);
	n_atoms[timegapii]++;
      }
    }
  }
  
  timetable = syst->displacement_times();
  
  output.open(filename.c_str(),ios::out | ios::app);
  
  for(timegapii=0;timegapii<n_timegaps;timegapii++)
  {
    correlation[timegapii]/=(float(n_atoms[timegapii])*mean()*mean());
    output << timetable[timegapii]<< "\t" << correlation[timegapii]<<"\n";
  }
  
  output.close();
}

#ifdef NEVEr

template <class valType>
Value_List<valType> Value_List<valType>::operator *(Trajectory_List t_list)const
{
  int system_timeii;
  int internal_time1, internal_time2;
  
  Value_List temp;
  
  update_size();
  
  temp.syst = syst;
  
  temp.time_conversion = new int [syst->show_n_timesteps()];
  temp.defined_times = new bool [syst->show_n_timesteps()];
  
  temp.n_times = n_times;
  
  /*loop to count the number of independent internal times the new Value_list must have*/
  for(system_timeii=0;system_timeii<syst->show_n_timesteps();system_timeii++)
  {
    internal_time1 = convert_time(system_timeii);
    internal_time2 = t_list.convert_time(system_timeii);   //problem with access

    if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
    {
      new_internal_time++;
    }

    last_internal_time1=internal_time1;
    last_internal_time2=internal_time2;
  }
  n_internal_times = new_internal_time;

  
  temp.included = new Boolean_List[n_internal_times];
  temp.values = new valType * [n_internal_times];
  for(int timeii=0;timeii<n_internal_times;timeii++)
  {
    temp.included[timeii].set(syst);
  }
  
  
  new_internal_time = -1;
  last_internal_time1 = -1;
  last_internal_time2 = -1;

  /*calculate new values based on intersection with trajectory*/
  for(system_timeii=0;system_timeii<sys->show_n_timesteps();system_timeii++)
  {
    internal_time1 = convert_time(system_timeii);
    internal_time2 = comparison.convert_time(system_timeii);

    if((internal_time1!=last_internal_time1)||(internal_time2!=last_internal_time2))
    {
      new_internal_time++;
      temp.included[new_internal_time] = included[internal_time1]&&t_list.included[internal_time2];
      for(int trajii=0;trajii<n_trajectories;trajii++)
      {
	temp.values[timeii][trajii]=values[timeii][trajii]*valType(temp.t_list[internal_time2](trajii));
      }
    }
    temp.time_conversion[system_timeii]=new_internal_time;	//set value in time conversion table to current new_internal_time for this system_timeii

    last_internal_time1=internal_time1;
    last_internal_time2=internal_time2;
  }

  
  
  for(int timeii=0;syst->show_n_trajectories();timeii++)
  {
    temp.defined_times[timeii]=defined_times[timeii];
  }
  
  
  return temp;
}

#endif

template <class valType>
string Value_List<valType>::write_pdb(int valtime, string file_name_stem, int time, string* atomnames)const
{
    string filename;
    string time_str;
    stringstream time_out;
    string atom_name;
    float x_coord;
    float y_coord;
    float z_coord;
   int type;
   int speciesID;
    time_out << time;
    time_str = time_out.str();

    filename = file_name_stem.append(".");
    filename = filename.append(time_str);
    filename = filename.append(".pdb");
    
    int timeii=convert_time(valtime);

    FILE * pdbfile;
    pdbfile = fopen (filename.c_str(),"w");

    if(pdbfile==NULL)
    {
        cout << "file  not opened"<< endl; cout.flush();
        exit(1);
    }

    for(int trajii=0;trajii<values[timeii].size();trajii++)
    {
        if (included[timeii](trajii))
        {


            speciesID = trajii+1;
            x_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_x();
            y_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_y();
            z_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_z();

            type = syst->show_trajectory(trajii)->show_type()-1;

//            cout<<"x_coord = "<<x_coord<<endl;cout.flush();
//            cout<<"y_coord = "<<y_coord<<endl;cout.flush();
//            cout<<"z_coord = "<<z_coord<<endl;cout.flush();
//            cout<<"speciesID = "<<speciesID<<endl;cout.flush();
//            cout<<"atomname = "<<atomnames[type]<<endl;cout.flush();
//            cout<<"value = "<<values[trajii]<<endl;cout.flush();

            fprintf(pdbfile, "HETATM%5i %4s MOL A   1    %8.4f%8.4f%8.4f  1.00%6f          %2s  \n", speciesID, atomnames[type].c_str(), x_coord, y_coord, z_coord, values[timeii][trajii], atomnames[type].c_str());


//fprintf(pdbfile, "HETATM%5i %4s MOL A   1    \n", speciesID, atomnames[type].c_str());
        }
    }
    fclose(pdbfile);


    return filename;

}


/*write in which atomnames are not defined by user*/
template <class valType>
string Value_List<valType>::write_pdb(int valtime, string file_name_stem, int time)const
{
    string filename;
    string time_str;
    stringstream time_out;
    string atom_name;
    float x_coord;
    float y_coord;
    float z_coord;
    int type;
    int speciesID;
    string atomnames [18];

    atomnames[0] = "H";
    atomnames[1] = "He";
    atomnames[2] = "Li";
    atomnames[3] = "Be";
    atomnames[4] = "B";
    atomnames[5] = "C";
    atomnames[6] = "N";
    atomnames[7] = "O";
    atomnames[8] = "F";
    atomnames[9] = "Ne";
    atomnames[10] = "Na";
    atomnames[11] = "Mg";
    atomnames[12] = "Al";
    atomnames[13] = "Si";
    atomnames[14] = "P";
    atomnames[15] = "S";
    atomnames[16] = "Cl";
    atomnames[17] = "Ar";

    int timeii=convert_time(valtime);
    
    time_out << time;
    time_str = time_out.str();

    filename = file_name_stem.append(".");
    filename = filename.append(time_str);
    filename = filename.append(".pdb");


    FILE * pdbfile;
    pdbfile = fopen (filename.c_str(),"w");

    cout<<"\n"<<values[timeii].size()<<"\t"<<timeii;
    for(int trajii=0;trajii<values[timeii].size();trajii++)
    {
   //   cout<<"\t"<<included[timeii](trajii);;
        if (included[timeii](trajii))
        {


            speciesID = trajii+1;
            x_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_x();
            y_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_y();
            z_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_z();

            type = syst->show_trajectory(trajii)->show_type()-1;

	    if(type>17){cout << "Error: More atom types in use than names available";exit(1);}

            fprintf(pdbfile, "HETATM%5i %4s MOL A   1    %8.4f%8.4f%8.4f  1.00%6f          %2s  \n", speciesID, atomnames[type].c_str(), x_coord, y_coord, z_coord, values[timeii][trajii], atomnames[type].c_str());


//fprintf(pdbfile, "HETATM%5i %4s MOL A   1    \n", speciesID, atomnames[type].c_str());
        }
    }
    fclose(pdbfile);


    return filename;

}



template <class valType>
string Value_List<valType>::write_pdb(int valtime, string file_name_stem, int time, string* atomnames,valType maxval)const
{
    string filename;
    string time_str;
    stringstream time_out;
    string atom_name;
    float x_coord;
    float y_coord;
    float z_coord;
    int type;
    int speciesID;
    time_out << time;
    time_str = time_out.str();
    valType val;

    filename = file_name_stem.append(".");
    filename = filename.append(time_str);
    filename = filename.append(".pdb");

    int timeii=convert_time(valtime);
    
    FILE * pdbfile;
    pdbfile = fopen (filename.c_str(),"w");


    for(int trajii=0;trajii<values[timeii].size();trajii++)
    {
        if (included[timeii](trajii))
        {

	    val = values[timeii][trajii];
	    if(val > maxval){val=maxval;}

            speciesID = trajii+1;
            x_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_x();
            y_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_y();
            z_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_z();

            type = syst->show_trajectory(trajii)->show_type()-1;

//            cout<<"x_coord = "<<x_coord<<endl;cout.flush();
//            cout<<"y_coord = "<<y_coord<<endl;cout.flush();
//            cout<<"z_coord = "<<z_coord<<endl;cout.flush();
//            cout<<"speciesID = "<<speciesID<<endl;cout.flush();
//            cout<<"atomname = "<<atomname[type]<<endl;cout.flush();
//            cout<<"value = "<<values[trajii]<<endl;cout.flush();

            fprintf(pdbfile, "HETATM%5i %4s MOL A   1    %8.4f%8.4f%8.4f  1.00%6f          %2s  \n", speciesID, atomnames[type].c_str(), x_coord, y_coord, z_coord, val, atomnames[type].c_str());


//fprintf(pdbfile, "HETATM%5i %4s MOL A   1    \n", speciesID, atomnames[type].c_str());
        }
    }
    fclose(pdbfile);


    return filename;

}


/*write in which atomnames are not defined by user*/
template <class valType>
string Value_List<valType>::write_pdb(int valtime, string file_name_stem, int time,valType maxval)const
{
    string filename;
    string time_str;
    stringstream time_out;
    string atom_name;
    float x_coord;
    float y_coord;
    float z_coord;
    int type;
    int speciesID;
    string atomnames [18];
    valType val;


    atomnames[0] = "H";
    atomnames[1] = "He";
    atomnames[2] = "Li";
    atomnames[3] = "Be";
    atomnames[4] = "B";
    atomnames[5] = "C";
    atomnames[6] = "N";
    atomnames[7] = "O";
    atomnames[8] = "F";
    atomnames[9] = "Ne";
    atomnames[10] = "Na";
    atomnames[11] = "Mg";
    atomnames[12] = "Al";
    atomnames[13] = "Si";
    atomnames[14] = "P";
    atomnames[15] = "S";
    atomnames[16] = "Cl";
    atomnames[17] = "Ar";
    
    int timeii=convert_time(valtime);

    time_out << time;
    time_str = time_out.str();

    filename = file_name_stem.append(".");
    filename = filename.append(time_str);
    filename = filename.append(".pdb");


    FILE * pdbfile;
    pdbfile = fopen (filename.c_str(),"w");


    for(int trajii=0;trajii<values[timeii].size();trajii++)
    {
        if (included[timeii](trajii))
        {

	    val = values[timeii][trajii];
	    if(val > maxval){val=maxval;}

            speciesID = trajii+1;
            x_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_x();
            y_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_y();
            z_coord = syst->show_trajectory(trajii)->show_coordinate(time).show_z();

            type = syst->show_trajectory(trajii)->show_type()-1;

	    if(type>17){cout << "Error: More atom types in use than names available";exit(1);}

            fprintf(pdbfile, "HETATM%5i %4s MOL A   1    %8.4f%8.4f%8.4f  1.00%6f          %2s  \n", speciesID, atomnames[type].c_str(), x_coord, y_coord, z_coord, val, atomnames[type].c_str());


//fprintf(pdbfile, "HETATM%5i %4s MOL A   1    \n", speciesID, atomnames[type].c_str());
        }
    }
    fclose(pdbfile);


    return filename;

}


template <class valType>
void Value_List<valType>::write_data_file(string filename)const
{
  
}



template <class valType>
bool Value_List<valType>::threshold(int timeii, int trajii,bool greater, valType threshold)const
{
    bool threshbool=0;
    if (greater)
    {
        if (included[timeii](trajii))
        {
            if (values[timeii][trajii] >= threshold)
            {
                threshbool = 1;
            }
        }
    }
    else
    {
        if (included[timeii](trajii))
        {
            if (values[timeii][trajii] < threshold)
            {
                threshbool = 1;
            }
        }
    }

    return threshbool;

}


template <class valType>
bool Value_List<valType>::threshold(int timeii, int trajii,valType low, valType high)const
{

        bool threshbool=0;
        if (included[timeii](trajii))
        {
            if (values[timeii][trajii]>=low)
            {
                if (values[timeii][trajii]<=high)
            {
                threshbool = 1;
            }
            }
            else
            {

            }
        }

        return threshbool;
}



template <class valType>
void Value_List<valType>::construct_t_list(bool greater, valType thresh, Trajectory_List* t_list)
{
    Boolean_List* thresholded;
    thresholded = new Boolean_List [n_times];
    bool threshbool=0;


    for (int timeii=0; timeii<n_times;timeii++)
    {
        thresholded[timeii].set(syst);
    }

    update_size();
    int n_trajectories = syst->show_n_trajectories();

        for (int timeii=0; timeii<n_times;timeii++)
        {
            for (int trajii=0;trajii<n_trajectories;trajii++)
            {
                threshbool = threshold(timeii,trajii,greater,thresh);
                thresholded[timeii](trajii,threshbool);

            }
       }



    t_list->set(syst, n_times, n_trajectories, thresholded, time_conversion);


    delete [] thresholded;

}




template <class valType>
void Value_List<valType>::construct_t_list(valType low, valType high,Trajectory_List* t_list)
{
    Boolean_List* thresholded;
    thresholded = new Boolean_List [n_times];
    bool threshbool=0;


    for (int timeii=0; timeii<n_times;timeii++)
    {
        thresholded[timeii].set(syst);
    }

    update_size();
    int n_trajectories = syst->show_n_trajectories();

    for (int timeii=0; timeii<n_times;timeii++)
    {
        for (int trajii=0;trajii<n_trajectories;trajii++)
        {
            threshbool = threshold(timeii,trajii,low,high);
            thresholded[timeii](trajii,threshbool);
        }
    }

    t_list->set(syst, n_times, n_trajectories, thresholded, time_conversion);

    delete [] thresholded;

}



template <class valType>
void Value_List<valType>::percentile_t_list(bool greater, float percentile, Trajectory_List* t_list)
{
    Boolean_List* thresholded;
    thresholded = new Boolean_List [n_times];
    bool threshbool=0;
    int timeii, trajii;
    
    vector<valType> temp_vals;
    
    valType * threshold;
    threshold = new valType [n_times];

    update_size();
    
    for(timeii=0; timeii<n_times; timeii++)
    {
      temp_vals.clear();
      for(trajii=0;trajii<values[timeii].size();trajii++)
      {
	if((included[timeii])(trajii))
	{
	  temp_vals.push_back(values[timeii][trajii]);
	}
      }
      sort(temp_vals.begin(),temp_vals.end());
      
      if(greater)
      {
	threshold[timeii]=temp_vals[int(ceil(temp_vals.size()*percentile/100.0))];
      }
      else
      {
	threshold[timeii]=temp_vals[int(floor(temp_vals.size()*percentile/100.0))];
      }
    }

    for (int timeii=0; timeii<n_times;timeii++)
    {
        thresholded[timeii].set(syst);
    }

    int n_trajectories = syst->show_n_trajectories();

        for (int timeii=0; timeii<n_times;timeii++)
        {
            for (int trajii=0;trajii<n_trajectories;trajii++)
            {
                threshbool = threshold(timeii,trajii,greater,threshold[timeii]);
                thresholded[timeii](trajii,threshbool);

            }
       }



    t_list->set(syst, n_times, n_trajectories, thresholded, time_conversion);


    delete [] thresholded;

}


template <class valType>
void Value_List<valType>::percentile_t_list(float low_percentile, float high_percentile,Trajectory_List* t_list)
{
    Boolean_List* thresholded;
    thresholded = new Boolean_List [n_times];
    bool threshbool=0;
    int timeii,trajii;
    
    vector<valType> temp_vals;
    
    valType * low;
    valType * high;
    low = new valType [n_times];
    high = new valType [n_times];
    
    update_size();
    
    if(low_percentile<0||high_percentile<0)
    {
      cout << "\nError: percentile thresholds must be non-negative.\n";
      exit(0);
    }
    else if(low_percentile>100||high_percentile>100)
    {
      cout << "\nError: percentile thresholds must be equal to or less than 100.\n";
      exit(0);
    }
    else if(low_percentile>high_percentile)
     {
      cout << "\nError: lower threshold must be less than higher threshold.\n";
      exit(0);
    }

    
    
    for(timeii=0; timeii<n_times; timeii++)
    {
      temp_vals.clear();
      for(trajii=0;trajii<values[timeii].size();trajii++)
      {
	if((included[timeii])(trajii))
	{
	  temp_vals.push_back(values[timeii][trajii]);
	}
      }
      sort(temp_vals.begin(),temp_vals.end());
      
      low[timeii]=temp_vals[int(ceil(temp_vals.size()*low_percentile/100.0))];
      high[timeii]=temp_vals[int(floor(temp_vals.size()*high_percentile/100.0))];

    }
    
    
    for (int timeii=0; timeii<n_times;timeii++)
    {
        thresholded[timeii].set(syst);
    }

    
    int n_trajectories = syst->show_n_trajectories();

    for (int timeii=0; timeii<n_times;timeii++)
    {
        for (int trajii=0;trajii<n_trajectories;trajii++)
        {
            threshbool = threshold(timeii,trajii,low,high);
            thresholded[timeii](trajii,threshbool);
        }
    }

    t_list->set(syst, n_times, n_trajectories, thresholded, time_conversion);

    delete [] thresholded;

}



template <class valType>
void Value_List<valType>::update_size()
{
  if(values.size()<syst->show_n_trajectories())
  {
    values.resize(syst->show_n_trajectories());
  }
}
}
#endif
