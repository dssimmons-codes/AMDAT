/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class to calculate n_fold order parameter*/
/*Written by Daniel A. Hunsicker*/


#include <stdlib.h>
#include "n_fold_order_parameter.h"
#include <math.h>
#include "version.h"
#include "trajectory_list.h"
#include "coordinate.h"
#include "value_list.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "progress.h"

using namespace std;
#include "tokenize.h"

N_Fold_Order_Parameter::N_Fold_Order_Parameter(System * sys, float ord, string sig_file, string orientation, float cut, string file_stem, int start, int end)
{
  /** constructor
  * @author Daniel Hunsicker
  * @date 6/16/2012
  **/

    system = syst = sys;
    total_atoms = system->show_n_atoms();
    start_time = start;
    end_time = end;

    n_times = system->show_n_timesteps();
    pdb_stem = file_stem;
    n_atomtypes=system->show_n_atomtypes();
    n_bins = 100;
    sig_cut = cut;
    atomcount_total=0;
    max_param=0;
    n_included=0;
    string plane;
    plane = orientation;
    order = ord;
    n_trajectories=system->show_n_trajectories();

    if(plane=="xy") distancefun = &Coordinate::length_xy;
    else if(plane=="xz") distancefun = &Coordinate::length_xz;
    else if(plane=="yz") distancefun = &Coordinate::length_yz;
    else
    {
    cout<<"Error: plane command "<<orientation<<" not understood.\n";
    exit(1);
    }

    set_time_conv();
    
    included = new Boolean_List[n_times];
    values = new float * [n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
      included[timeii].set(syst);
      values[timeii] = new float [n_trajectories];
      for(int trajii=0;trajii<n_trajectories;trajii++)
      {
	values[timeii][trajii]=0.0;
      }
    }

    param_total_current = new float [n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
        param_total_current[timeii]=0;
    }

    distribution = new float [n_bins+2];

    for (int binii=0;binii<=n_bins+1;binii++)
    {
        distribution[binii]=0;
    }

    //allocate sigma matrix
    string line;
    line = "";
    int sig_tokens=0;
    string * sig_ARGS;
    sig_ARGS =new string [system->show_n_atomtypes()+1];
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


    max_sigma=0;

    for (int sig1ii=0; sig1ii<sig_tokens-1; sig1ii++)
    {
        for (int sig2ii=0; sig2ii<sig_tokens-1; sig2ii++)
        {
            if (max_sigma < sigmatrix[sig1ii][sig2ii])
            {
                max_sigma = sigmatrix[sig1ii][sig2ii];
            }
        }
    }


}

void N_Fold_Order_Parameter::set_time_conv()
{
  /** sets time conversion table
  * @author Daniel Hunsicker
  * @date 6/16/2012
  **/
    time_conversion = new int [n_times];
    defined_times = new bool [n_times];
    
    for (int timeii=0; timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
	defined_times[timeii]=1;

    }
}

void N_Fold_Order_Parameter::analyze(Trajectory_List * t_list)
{
  /** calculates nfold for a system for all times with trajectory lists
  * @author Daniel Hunsicker
  * @date 6/16/2012
  **/

    param_total=0;

    trajectory_list=t_list;

    for(int timeii=0;timeii<n_times;timeii++)
    {

        current_time = timeii;


        atomcount=0;

        trajectory_list->listloop(this,current_time);

        param_total_current[current_time] = param_total_current[current_time]/float(atomcount);

        atomcount_total = atomcount_total + atomcount;

        print_progress(timeii, n_times-1);

    }


    time_average_param = param_total/float(atomcount_total);
    for (int binii=0; binii<n_bins+2;binii++)
    {
        distribution[binii] = distribution[binii]/float(atomcount_total);
    }


}

void N_Fold_Order_Parameter::listkernel(Trajectory * current_trajectory)
{
  /** calculates nfold for a particle with trajectory lists
  * @author Daniel Hunsicker
  * @date 6/16/2012
  **/




    int jtype;
    int ktype;
    int j_ID;
    float sigma_jk;
    Coordinate j_coord;
    Coordinate k_coord;
    Coordinate dist;
    float theta;
    float cutoff;
    float jk_distance;
    float A;
    float B;
    float Atotal;
    float Btotal;
    int binID;

    param_jth=0;

    jtype=0;
    neighborcount=0;
    sigma_jk=0;
    cutoff=0;
    Atotal=0;
    Btotal=0;
    jtype = current_trajectory->show_type()-1;
    j_coord = current_trajectory->show_coordinate(current_time);
    j_ID = current_trajectory->show_trajectory_ID();



    Coordinate const* bounds;
    bounds = system->boundaries(current_time);
    Coordinate box_size;
    box_size = system->size(current_time);

    cutoff = sig_cut * max_sigma;

    atomcount++;
            for (int kii=0;kii<total_atoms;kii++)
            {
                ktype=0;


                Trajectory* kth_trajectory;
                kth_trajectory=system->show_trajectory(kii);
                ktype=kth_trajectory->show_type()-1;
                k_coord = kth_trajectory->show_coordinate(current_time);


                sigma_jk = sigmatrix[jtype][ktype];

                dist = k_coord - j_coord;
                dist -= box_size * ((dist/(box_size*.5)).integer());

                jk_distance = (dist.*distancefun)();
                cutoff = sig_cut * sigma_jk;


                if (jk_distance < cutoff && jk_distance > 0)
                {

                    A=0;
                    B=0;
                    neighborcount = neighborcount++;


                    theta = atan2((dist.show_y()),(dist.show_x()));

                    A = cos(order*theta);
                    B = sin(order*theta);


                    Atotal+=A;
                    Btotal+=B;




                }
            }

            param_jth += sqrt(pow(Atotal,2.0) + pow(Btotal,2.0))/order;



            if (param_jth > max_param)
            {
                max_param=param_jth;
            }

            binID = int(param_jth*n_bins);
            if (binID>=n_bins)
            {
                binID=n_bins+1;
            }
            distribution[binID]++;

            param_total= param_total + param_jth;
            param_total_current[current_time] = param_total_current[current_time] + param_jth;

            set(current_time,j_ID,param_jth);

}


void N_Fold_Order_Parameter::write(string filename)
{
  /** Writes to data file
  * @author Daniel Hunsicker
  * @date 6/16/2012
  **/
    string map_file;
    string thresh_pdb_stem;
    string map_thresh_file;
    cout << "\nWriting hop to file " << filename << "." << endl;

    ofstream output(filename.c_str());

    output << "Hexatic order parameter data created by AMDAT v." << VERSION << endl << endl;

    output << "The time average n_fold order parameter for the system is " << time_average_param << endl << endl;


    for(int timeii=0;timeii<n_times;timeii++)
    {

        if( timeii>=start_time && timeii<=end_time )
        {
            map_file = write_pdb(timeii,pdb_stem, timeii, species_name);
            output << "Map of the n_fold order parameter was written to "<< map_file << endl;
        }

    }

    output << endl;

    float binsize;
    binsize=1/float(n_bins);

    output << "distribution of the average n_fold order parameter over time." << endl;
    output << endl;


    for(int binii=0;binii<n_bins+1;binii++)
    {
        float bin;
        bin = float(binii)*binsize;
        output << bin << "\t" << distribution[binii] << endl;

    }
        output << "overflow\t" << distribution[n_bins+1]<<endl;


        output <<endl;

        output << "max n_fold = " << max_param << endl;

}
