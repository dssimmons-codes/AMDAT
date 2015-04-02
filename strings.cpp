/*Methods for String class - Identifies strings in fluid*/
/*Molecular dynamics analysis toolkit (MDAT)*/
/*Written by David S. Simmons*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include "strings.h"
#include "version.h"

using namespace std;

#include "tokenize.h"
#include "system.h"

/*-----------------------------------------------------------------------------------*/


/*Constructor*/
Strings::Strings(System * sys, int timegap1, float thresh, string sigmatrixname, int maximum_strings, int maximum_stringatoms)
{
	initialize(sys, timegap1, thresh, sigmatrixname, maximum_strings, maximum_stringatoms);
}


/*-----------------------------------------------------------------------------------*/



/*Method to initialize object*/
void Strings::initialize(System * sys, int timegap1, float thresh, string sigmatrixname,int maximum_strings, int maximum_stringatoms)
{
	int timeii, stringii, trajii, lengthii;

	system = sys;

	/*set parameters describing string arrays*/
	maxstrings = maximum_strings;
	maxstringatoms = maximum_stringatoms;
	timegap=timegap1;
	n_trajectories = system->show_n_atoms()+system->show_n_molecules();

	allocate_sig_matrix(sigmatrixname);
	
	threshold=thresh;
	
	n_times=system->timegap_weighting(timegap);
	timetable=system->displacement_times(timegap);

	/*memory allocation for string data.*/
	stringcount = new int [n_times];
	stringID = new int * [n_times];
	stringlist = new int ** [n_times];
	stringsize = new int * [n_times];
	length_distribution = new float [maxstringatoms+1];
	mean_strings=0;
	mean_length=0;
	order_parameter=0;
	trajectories_considered=0;
	
		for(timeii=0;timeii<n_times;timeii++)
		{
			stringID[timeii]=new int[n_trajectories];
			stringcount[timeii]=0;
			stringlist[timeii] = new int * [maxstrings];
			stringsize[timeii] = new int [maxstrings];
			for(trajii=0;trajii<n_trajectories;trajii++)
			{
				stringID[timeii][trajii]=-1;

			}
			for(stringii=0;stringii<maxstrings;stringii++)
			{
				stringsize[timeii][stringii]=0;	//initialize string size to zero
			}
		}
		for(lengthii=0;lengthii<=maxstringatoms;lengthii++)
		{
			length_distribution[lengthii]=0;
		}
	
}


 /*---------------------------------------------------------*/
 /*------Methods to find strings via trajectory lists-------*/
 /*---------------------------------------------------------*/


 void Strings::analyze(Trajectory_List*t_list)
{
	trajectory_list=t_list;
	system->displacement_list(this, timegap);
	postprocess_list();
}


/*-----------------------------------------------------------------------------------*/


 void Strings::list_displacementkernel(int timegapii, int thisii, int nextii)
{
	nexttime = nextii;
	thistime=thisii;
	current_timegap=timegapii;
	trajectories_considered+=(trajectory_list[0]).show_n_trajectories(thistime);
	(trajectory_list[0]).listloop(this,thistime);
}



/*-----------------------------------------------------------------------------------*/


 void Strings::listkernel(Trajectory* trajectory1)
{
	int trajii=0;
	Trajectory * trajectory2;
	int trajectory1ID, trajectory2ID;
	int thisii = int(float(thistime)/float(system->show_n_exponential_steps()));
	int temp_stringID, temp_stringsize, temp_stringID2;
	trajectory1ID=trajectory1->show_trajectory_ID();
	float distance;
	int trajtype1, trajtype2;
	int n_listed_trajectories = trajectory_list->show_n_trajectories(thistime);
	

	//reworking this SECTION
	//while(trajectory2!=0&&trajectory2!=trajectory1)	//only perform analysis if the trajectory list returns a trajectory and if trajectories are distinct
	for(trajii=0;trajii<n_listed_trajectories;trajii++)
	{
		trajectory2 = (*trajectory_list)(thistime, trajii);
		
		if(trajectory2!=trajectory1)
		{
		  trajtype1 = trajectory1->show_type()-1;
		  trajtype2 = trajectory2->show_type()-1;
	  
		  trajectory2ID=trajectory2->show_trajectory_ID();
		  
		  distance = (trajectory1->show_coordinate(thistime)-trajectory2->show_coordinate(nexttime)).length_unwrapped(system->size(thistime));
		  
		  if(distance<threshold*sigmatrix[trajtype1][trajtype2])	//check if pairing satisfies criterion for string
		  {
		    if(stringID[thisii][trajectory1ID]==-1 && stringID[thisii][trajectory2ID]==-1)	//if neither atom is in a string
		    {
		      /*create new string*/
		      create_string(thisii,trajectory1ID,trajectory2ID);
		    }
		    else if(stringID[thisii][trajectory1ID]==-1 && stringID[thisii][trajectory2ID]>=0)	//if atom 2 is in a string but atom 1 is not
		    {
		      /*add atom 1 to string containing atom2*/
		      temp_stringID=stringID[thisii][trajectory2ID];
		      temp_stringsize=stringsize[thisii][temp_stringID];
		      /*check if string is full*/
		      if(temp_stringsize==maxstringatoms)
		      {
			cout << "Error: string size limit reached. Please increase string size limit.\n";
			exit(1);
		      }
		      stringID[thisii][trajectory1ID]=temp_stringID;
		      stringlist[thisii][temp_stringID][temp_stringsize]=trajectory1ID;
		      stringsize[thisii][temp_stringID]++;
		    }
		    else if(stringID[thisii][trajectory1ID]>=0 && stringID[thisii][trajectory2ID]==-1)	//if atom 1 is in a string but atom 2 is not
		    {
		      /*add atom 2 to string containing atom 1*/
		      temp_stringID=stringID[thisii][trajectory1ID];
		      temp_stringsize=stringsize[thisii][temp_stringID];
		      /*check if string is full*/
		      if(temp_stringsize==maxstringatoms)
		      {
			cout << "Error: string size limit reached. Please increase string size limit.\n";
			exit(1);
		      }
		      stringID[thisii][trajectory2ID]=temp_stringID;
		      stringlist[thisii][temp_stringID][temp_stringsize]=trajectory2ID;
		      stringsize[thisii][temp_stringID]++;
		    }
		    else if(stringID[thisii][trajectory1ID]>=0 && stringID[thisii][trajectory2ID]>=0 && stringID[thisii][trajectory1ID]!=stringID[thisii][trajectory2ID])	//if both atoms are in strings but not the same string
		    {
		      /*merge strings*/
		      temp_stringID = stringID[thisii][trajectory1ID];
		      temp_stringID2 = stringID[thisii][trajectory2ID];
		      merge_strings(thisii,temp_stringID,temp_stringID2);
		    }
		  }
		//get next trajectory
		//trajii++;
		//trajectory2 = (trajectory_list[0])(trajii);
		}
	}
}


/*-----------------------------------------------------------------------------------*/


 void Strings::postprocess_list()
{
	defragment_strings();
	calculate_mean_strings();
	calculate_mean_length();
	stringlength_distribution();
	fraction_in_strings();
}


 /*---------------------------------------------------------*/
 /*-------------String array management methods-------------*/
 /*---------------------------------------------------------*/

void Strings::create_string(int thisii, int atom1ID, int atom2ID)
{
	int temp_stringID=stringcount[thisii];	//determine new stringID

	stringID[thisii][atom1ID]=temp_stringID;	//assign this stringID to atom1
	stringID[thisii][atom2ID]=temp_stringID;	//assign this stringID to atom2
	stringlist[thisii][temp_stringID] = new int [maxstringatoms];	//allocate memory for this new string
	stringlist[thisii][temp_stringID][0] = atom1ID;	//put atom1 in string
	stringlist[thisii][temp_stringID][1] = atom2ID;	//put atom2 in string

	stringsize[thisii][temp_stringID]=2;			//set string size to 2

	stringcount[thisii]++;		//increment next open string ID
}



/*-----------------------------------------------------------------------------------*/


/*the problem is here*/
void Strings::merge_strings(int thisii, int stringID1, int stringID2)
{
	int atomii, temp_atomID,temp_stringsize;

	for(atomii=0;atomii<stringsize[thisii][stringID2];atomii++)
	{
		//cout << "\n" << atomii << stringlist[timegapii][thisii][stringID2][atomii]; cout.flush();
		temp_atomID=stringlist[thisii][stringID2][atomii];
		temp_stringsize=stringsize[thisii][stringID1];
		stringID[thisii][temp_atomID]=stringID1;	//switch stringID of atoms in second string to first string
		stringlist[thisii][stringID1][temp_stringsize]=temp_atomID;	//add atom in second string to first string
		stringsize[thisii][stringID1]++;	//increment size of first string
	}

	delete [] stringlist[thisii][stringID2];	//delete second string
	stringsize[thisii][stringID2]=0;		//set second string's size to zero

	/*if deleted string is last string in list, decrement count of strings so that a string list never ends in an empty string*/
	if(stringID2==stringcount[thisii]-1)
	{
		stringcount[thisii]--;
	}
}


/*-------------------------------------------------------*/
/*---------String postprocessing and statistics----------*/
/*-------------------------------------------------------*/


/*method to find empty strings, eliminate them, and consolidate string list accordingly*/
void Strings::defragment_strings()
{
	int timeii, stringii, atomii, temp_stringcount, temp_stringsize, temp_atomID;
	int stringsreduced;		//number of empty strings found so far
	
		for(timeii=0;timeii<n_times;timeii++)
		{
			stringsreduced = 0;		//set number of empty strings found so far to zero
			temp_stringcount = stringcount[timeii];
			for(stringii=0;stringii<temp_stringcount;stringii++)
			{
				if(stringsize[timeii][stringii]==0)	//if string is empty
				{
					stringsreduced++;	//increment number of empty strings found
				}
				else
				{
				  //cout<<"\n"<<stringsize[timegapii][timeii][stringii];
				  if(stringsreduced>0)
				  { 
				    stringlist[timeii][stringii-stringsreduced]=stringlist[timeii][stringii];	//shift string pointers down
				    stringsize[timeii][stringii-stringsreduced]=stringsize[timeii][stringii];	//shift list of string sizes down
				    temp_stringsize=stringsize[timeii][stringii];	//make note of the new string size
				    /*shift stringID of atoms in string downward appropriately*/
				    for(atomii=0;atomii<temp_stringsize;atomii++)
				    {
				      temp_atomID=stringlist[timeii][stringii-stringsreduced][atomii];
				      stringID[timeii][temp_atomID]=stringii-stringsreduced;
				    }
				  }
				}			
			}
			stringcount[timeii]-=stringsreduced;	//reduce count of strings at this timegap and time by number of empty strings found
			

//			/*delete duplicate strings now residing at end of list*/
//			for(stringii=temp_stringcount-stringsreduced;stringii<temp_stringcount;stringii++)
//			{
//				delete [] stringlist[timegapii][timeii][stringii];	//delete string
//				stringsize[timegapii][timeii][stringii]=0;		//set string's size to zero
//			}

			
			
		}
}


/*-----------------------------------------------------------------------------------*/


void Strings::calculate_mean_strings()
{
	int timeii;
	mean_strings=0;

		for(timeii=0;timeii<n_times;timeii++)
		{
			mean_strings+=float(stringcount[timeii])/float(n_times);	
		}
}


/*-----------------------------------------------------------------------------------*/


void Strings::calculate_mean_length()
{
	int timeii,stringii;
	//int * n_strings;

	//n_strings = new int [n_timegaps];

		mean_length=0;
		for(timeii=0;timeii<n_times;timeii++)
		{
			for(stringii=0;stringii<stringcount[timeii];stringii++)
			{
				mean_length+=float(stringsize[timeii][stringii])/float((mean_strings*n_times));	//tally average number of atoms in strings at this timegap
			}
		}
}


/*-----------------------------------------------------------------------------------*/


void Strings::stringlength_distribution()
{
	int timeii,stringii;
	int temp_length;
		for(timeii=0;timeii<n_times;timeii++)
		{
			for(stringii=0;stringii<stringcount[timeii];stringii++)
			{
				temp_length = stringsize[timeii][stringii];
				if(temp_length>=maxstringatoms)
				{
					temp_length = maxstringatoms;
				}
				if(mean_strings*n_times!=0)
				{
					length_distribution[temp_length]+=1/(mean_strings*n_times);
				}
				//cout <<"\n"<<timegapii << "\t" <<temp_length<<"\t"<< length_distribution[timegapii][temp_length];
			}
		}
}


/*-----------------------------------------------------------------------------------*/


/*calculate fractions of particles in strings*/
void Strings::fraction_in_strings()
{

		order_parameter = n_times*mean_strings*mean_length/float(trajectories_considered);

}


/*-----------------------------------------------------------------------------------*/



/*allocate matrix of particle sizes and assign values*/
void Strings::allocate_sig_matrix(string sig_file)
{
  Tokenize tokenize;
  
    string line;
    line = "";
    int sig_tokens=0;
    string * sig_ARGS;
    
    n_atomtypes = system->show_n_atomtypes();
    sig_ARGS =new string [n_atomtypes+1];

    string * species_name;
    species_name = new string [0];

    ifstream file(sig_file.c_str());

    if (file.is_open())
    {
        //get first line of matrix
        getline (file,line);
        sig_tokens = tokenize(line, sig_ARGS);
        species_name = new string [sig_tokens];
        species_name[0] = sig_ARGS[0];

        // initialize matrix
        sigmatrix = new float * [n_atomtypes];
        for (int argii=0; argii<n_atomtypes; argii++)
        {
            sigmatrix[argii] = new float [n_atomtypes];
        }


        //set all elements to 0
        for (int arg1ii=0;arg1ii<n_atomtypes;arg1ii++)
        {
            for (int arg2ii=0;arg2ii<n_atomtypes;arg2ii++)
            {
                sigmatrix[arg1ii][arg2ii] = 0;
            }
        }

        //input first row of the matrix
        for(int argsii=1; argsii<=system->show_n_atomtypes(); argsii++)
        {
            stringstream ss(sig_ARGS[argsii]);
            ss>>sigmatrix[0][argsii-1];
        }

        //input rest of the matrix
        for(int lineii=1;lineii<system->show_n_atomtypes(); lineii++)
        {

            getline (file,line);
            sig_tokens = tokenize(line, sig_ARGS);

            species_name[lineii] = sig_ARGS[0];

            for(int argsii=1; argsii<=system->show_n_atomtypes(); argsii++)
            {
                stringstream ss(sig_ARGS[argsii]);
                ss>>sigmatrix[lineii][argsii-1];
            }
        }
    }
    else
    {
        cout << "\nError: sigma data file not opened succesfully.\n";
        exit(1);
    }
    file.close();
//
//
//    for (int sig1ii=0; sig1ii<n_atomtypes; sig1ii++)
//    {
//        for (int sig2ii=0; sig2ii<n_atomtypes; sig2ii++)
//        {
//            cout << "sigma " << sig1ii+1 << sig2ii+1<< "= "<< sigmatrix[sig1ii][sig2ii]<< endl; cout.flush();
//        }
//    }

}



/*-------------------------------------------------------*/
/*---------------------Data output-----------------------*/
/*-------------------------------------------------------*/

void Strings::write(string filename)const
{
	int lengthii;

	cout << "\nWriting string data to file " << filename <<".";
	ofstream output(filename.c_str());
	output << "String data created by AMDAT v." << VERSION << "\n";

	output << "time\tmean_strings\tmean_length\tmean_length(counting1)\torder_parameter";
	for(lengthii=0;lengthii<=maxstringatoms;lengthii++)
	{
		output << "\t" << lengthii;
	}

	output << "\n";
	
		output << timetable << "\t" << mean_strings << "\t" << mean_length<<"\t"<<order_parameter*mean_length+(1-order_parameter)<<"\t"<<order_parameter;
		
		for(lengthii=0;lengthii<=maxstringatoms;lengthii++)
		{
			output << "\t" << length_distribution[lengthii];
		}
		output << "\n";

}


/*-----------------------------------------------------------------------------------*/


void Strings::write_jump(string filename,int timeindex)const
{
  int stringii, trajii;
  int type;
  int n_traj = 0;
  Coordinate tempcoordinate;
  ofstream output(filename.c_str());

  for(stringii=0;stringii<stringcount[timeindex];stringii++)
  {
    for(trajii=0;trajii<stringsize[timeindex][stringii];trajii++)
    {
      n_traj++;
    }
  }

  output << n_traj << "\nAtoms\n";
  for(stringii=0;stringii<stringcount[timeindex];stringii++)
  {
    for(trajii=0;trajii<stringsize[timeindex][stringii];trajii++)
    {
      tempcoordinate = (system->show_trajectory(stringlist[timeindex][stringii][trajii]))->show_coordinate(timeindex*system->show_n_exponential_steps());
      type = (system->show_trajectory(stringlist[timeindex][stringii][trajii]))->show_type();
      output <<  type  <<  "\t"  << tempcoordinate.show_x() << "\t" << tempcoordinate.show_y() << "\t" << tempcoordinate.show_z() << "\n";
    }
  }

  output << n_traj << "\nAtoms\n";
  for(stringii=0;stringii<stringcount[timeindex];stringii++)
  {
    for(trajii=0;trajii<stringsize[timeindex][stringii];trajii++)
    {
      tempcoordinate = (system->show_trajectory(stringlist[timeindex][stringii][trajii]))->show_coordinate(timeindex*system->show_n_exponential_steps()+timegap);
      type = (system->show_trajectory(stringlist[timeindex][stringii][trajii]))->show_type();
      output <<  type  <<  "\t"  << tempcoordinate.show_x() << "\t" << tempcoordinate.show_y() << "\t" << tempcoordinate.show_z() << "\n";;
    }
  }

}

void Strings::construct_trajectory_list(Trajectory_List* t_list)
{
    Boolean_List* in_string;
    in_string = new Boolean_List [n_times];
    bool threshbool = 0;
    int timeii,stringii, atomii;
    int n_trajectories = system->show_n_trajectories();
    int * time_conversion;
    
    for(timeii=0; timeii<n_times;timeii++)
    {
        in_string[timeii].set(system);
    }

    for(timeii=0;timeii<n_times;timeii++)
    {
      for(stringii=0;stringii<stringcount[timeii];stringii++)
      {
	for(atomii=0;atomii<stringsize[timeii][stringii];atomii++)
	{
	  in_string[timeii](stringlist[timeii][stringii][atomii],1);
	}
      }
    }
 
   time_conversion=new int [system->show_n_timesteps()];
   for(timeii=0;timeii<system->show_n_timesteps();timeii++)
   {   
     time_conversion[timeii]=int(float(timeii-system->show_frt())/float(system->show_n_exponential_steps()));
     if(time_conversion[timeii]>n_times)
     {
       time_conversion[timeii]=n_times;
     }
   }
  
  t_list->set(system, n_times, n_trajectories, in_string, time_conversion);
}
