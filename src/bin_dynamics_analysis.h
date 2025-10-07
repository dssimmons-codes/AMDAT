// Copyright (C) 2010-2025 David S. Simmons

/*Bin_Dynamics_Analysis - performs analysis to binned trajectories
Amorphous Molecular Dynamics Tookit (AMDAT)
Written by Mark Mackura*/

#ifndef BIN_DYNAMICS_ANALYSIS_H
#define BIN_DYNAMICS_ANALYSIS_H

#include "analysis.h"
#include "trajectory_list_bins.h"
#include "system.h"
#include "trajectory_list.h"
#include "boolean_list.h"
#include <string>
#include <sstream>

using namespace std;

template <class Analysis_type>
class Bin_Dynamics_Analysis : public Analysis
{
    Analysis_type*** analysis;
    Trajectory_List_Bins * traj_list_bins;
    Boolean_List * included;
    int * time_conversion;
    int xii,yii,zii;
    
    Trajectory_List last_traj_list;
    
    public:
        //Bin_Dynamics_Analysis();
        Bin_Dynamics_Analysis(Trajectory_List_Bins*,Analysis_type,System*);

        void analyze(Trajectory_List*);
        void list_displacementkernel(int timegap,int time1,int time2);
        void write(string);

};

template <class Analysis_type>
Bin_Dynamics_Analysis<Analysis_type>::Bin_Dynamics_Analysis(Trajectory_List_Bins* t_l_b, Analysis_type analysis_input, System * sys)
{
    /** Create array of analysis objects, one for each bin 
    * @param t_l_b - binning structure for analysis
    * @param analysis_input - analysis object for building the array
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;
    traj_list_bins = t_l_b;
//    Trajectory_List temp(system->show_n_timesteps(),system->show_n_trajectories());
//    last_traj_list = temp;
    analysis = new Analysis_type**[traj_list_bins->show_n_xbins()];
    for(xii=0;xii<traj_list_bins->show_n_xbins();xii++)
    {
        analysis[xii] = new Analysis_type*[traj_list_bins->show_n_ybins()];
        for(yii=0;yii<traj_list_bins->show_n_ybins();yii++)
        {
            analysis[xii][yii] = new Analysis_type[traj_list_bins->show_n_zbins()];
            for(zii=0;zii<traj_list_bins->show_n_zbins();zii++)
            {
                analysis[xii][yii][zii] = analysis_input;
                //analysis[xii][yii][zii].set.(system);
            }
        }
    }
}

template <class Analysis_type>
void Bin_Dynamics_Analysis<Analysis_type>::analyze(Trajectory_List * t_list)
{
    /** Performs analysis
    * @param t_list - list of trajectories for calculation to be done on
    * @author Mark Mackura
    * @date 4/11/2012
    **/
    included = new Boolean_List[system->show_n_timesteps()];
    trajectory_list = t_list;       	//specific subset of particles to analyze
    
    /* set time conv table for last_traj_list */
    time_conversion = new int[system->show_n_timesteps()];
    for(int timeii=0;timeii<system->show_n_timesteps();timeii++)
    {
        time_conversion[timeii]=timeii;
    }
    
    /* set last_traj_list booleans to all true */
    for(int timeii=0; timeii<system->show_n_timesteps(); timeii++)
    {
      included[timeii].set(system);
      for(int trajii=0; trajii<system->show_n_trajectories(); trajii++)
	{
	  included[timeii](trajii,1);
	}
    }
    /* calculations */
    cout<<endl;cout.flush();
    
    system->displacement_list(this);	//calls system to loop over displacement times, in turn calls list_displacementkernel below
    
    cout<<endl;cout.flush();
    /* postprocessing */
    for(xii=0;xii<traj_list_bins->show_n_xbins();xii++)
    {
        for(yii=0;yii<traj_list_bins->show_n_ybins();yii++)
        {
            for(zii=0;zii<traj_list_bins->show_n_zbins();zii++)
            {
                analysis[xii][yii][zii].postprocess_bins();
            }
        }
    }
}

template <class Analysis_type>
void Bin_Dynamics_Analysis<Analysis_type>::list_displacementkernel(int timegapii, int thisii, int nextii)
{
    /** Kernel for loop over time within bin
    * @param timegapii - system time between frames
    * @param thisii - frame # that begins timegap
    * @param nextii - frame # that ends timegap
    * @note Will perform analysis only on trajectories that are in bins for the duration of the timegap
    * @author Mark Mackura
    * @date 4/11/2012
    **/
    
    Trajectory_List temp_traj_list;
    
    for(xii=0;xii<traj_list_bins->show_n_xbins();xii++)
    {
        for(yii=0;yii<traj_list_bins->show_n_ybins();yii++)
        {
            for(zii=0;zii<traj_list_bins->show_n_zbins();zii++)
            {
	      cout<<"Bin ("<<xii<<","<<yii<<","<<zii<<")"<<"\r";cout.flush();
	      temp_traj_list=(*trajectory_list &&((*traj_list_bins)(xii,yii,zii).time_intersection(thisii,nextii)));
	      analysis[xii][yii][zii].bin_hook(&temp_traj_list,timegapii ,thisii ,nextii);
	      //last_traj_list.set(system,system->show_n_timesteps(),system->show_n_trajectories(),included,time_conversion);	//reset last_traj_list to true for all particles for new bin
	      //last_traj_list = *trajectory_list && last_traj_list;
	      
	    }
        }
    } 
    
}

template <class Analysis_type>
void Bin_Dynamics_Analysis<Analysis_type>::write(string filename_template)
{
    /**  Writes output file for each bin
    * @param filename_template - custome filename to be appended to end of bin indicies
    * @author Mark Mackura
    * @date 4/11/2012
    */
    string bin_filename;
    stringstream ssx,ssy,ssz;
    for(xii=0;xii<traj_list_bins->show_n_xbins();xii++)
    {
        for(yii=0;yii<traj_list_bins->show_n_ybins();yii++)
        {
            for(zii=0;zii<traj_list_bins->show_n_zbins();zii++)
            {
                ssx << xii;
		ssx.seekp(0);
		ssy << yii;
		ssy.seekp(0);
		ssz << zii;
		ssz.seekp(0);
                bin_filename = filename_template+"_bin_"+ssx.str()+"."+ssy.str()+"."+ssz.str()+".bindata";
                analysis[xii][yii][zii].write(bin_filename);
            }
        }
    }
}


#endif // BIN_DYNAMICS_ANALYSIS_H
