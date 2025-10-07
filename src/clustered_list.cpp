// Copyright (C) 2010-2025 David S. Simmons

#include "clustered_list.h"
#include "tokenize.h"
#include <sstream>
#include "progress.h"
#include <stdlib.h>

using namespace std;


Clustered_List::Clustered_List()
{
	capacity = 0;
    n_atomtypes = 0;
	n_times=0;
	trajectories = new Trajectory ** [n_times];
	n_trajectories = new int [n_times];
	included = new Boolean_List [n_times];
	for(int timeii=0;timeii<n_times;timeii++)
	{
		trajectories[timeii] = new Trajectory * [capacity];
		n_trajectories[timeii]=0;
	}

	//set all system times to convert to 0 internal time
	time_conversion=new int [0];
    neighbor = new int** [n_times];
   for (int timeii=0; timeii<n_times; timeii++)
    {
        neighbor[timeii] = new int* [n_trajectories[timeii]];
        for (int trajii=0; trajii<n_trajectories[timeii];trajii++)
        {
            neighbor[timeii][trajii] = new int  [10];
            for (int ii=0; ii<10; ii++)
            {
                neighbor[timeii][trajii][ii] = 0;
            }

        }
    }

            sigmatrix = new float* [n_atomtypes];
    for (int atomtypeii =0; atomtypeii < n_atomtypes; atomtypeii++)
    {
        sigmatrix[atomtypeii] = new float [n_atomtypes];
        for (int atomtype2ii =0; atomtype2ii < n_atomtypes; atomtype2ii++)
        {
            sigmatrix[atomtypeii][atomtype2ii] = 0;
        }
    }

}


Clustered_List::~Clustered_List()
{
    for (int timeii=0; timeii<n_times; timeii++)
    {
        for (int trajii=0; trajii<n_trajectories[timeii];trajii++)
        {

            delete [] neighbor[timeii][trajii];

        }
        delete [] neighbor[timeii];
    }
    delete [] neighbor;


    for (int atomtypeii =0; atomtypeii < n_atomtypes; atomtypeii++)
    {
        delete [] sigmatrix[atomtypeii];
    }
    delete [] sigmatrix;
}



void Clustered_List::initialize_members()
{
    neighbor = new int** [n_times];

    for (int timeii=0; timeii<n_times; timeii++)
    {
        neighbor[timeii] = new int* [n_system_trajectories()];
        for (int trajii=0; trajii<n_system_trajectories();trajii++)
        {
            neighbor[timeii][trajii] = new int  [10];
            for (int ii=0; ii<10; ii++)
            {
                neighbor[timeii][trajii][ii] = 0;
            }

        }
    }

            sigmatrix = new float* [n_atomtypes];
    for (int atomtypeii =0; atomtypeii < n_atomtypes; atomtypeii++)
    {
        sigmatrix[atomtypeii] = new float [n_atomtypes];
        for (int atomtype2ii =0; atomtype2ii < n_atomtypes; atomtype2ii++)
        {
            sigmatrix[atomtypeii][atomtype2ii] = 0;
        }
    }

}




void Clustered_List::construct_clust_list(Trajectory_List* t_list, string sigma_file, string plane,int primary, int secondary)
{


    allocate_sig_matrix(sigma_file);


    Boolean_List* clustered;
    Boolean_List* previous;
     Boolean_List* previous2;
    clustered = new Boolean_List [sys->show_n_timesteps()];
    previous = new Boolean_List [sys->show_n_timesteps()];
    previous2 = new Boolean_List [sys->show_n_timesteps()];
    bool clustbool;
    int current_time;


    for (int timeii=0; timeii<sys->show_n_timesteps();timeii++)
    {
        current_time  = convert_time(timeii);

        clustered[current_time].set(sys);
        previous[current_time].set(sys);
        previous2[current_time].set(sys);
    }


    cout << "Generating neighbor lists" << endl;cout.flush();
    for (int timeii=0; timeii<sys->show_n_timesteps(); timeii++)
    {current_time  = convert_time(timeii);

        for (int trajii=0; trajii<sys->show_n_trajectories();trajii++)
        {
            generate_neighbor_lists(current_time,trajii,plane);
        }
                print_progress(timeii, sys->show_n_timesteps()-1);
    }

    cout << endl<<"Finding primary members" << endl;cout.flush();

    for (int timeii=0; timeii<sys->show_n_timesteps();timeii++)
    {
        current_time  = convert_time(timeii);


        for (int trajii=0; trajii<n_system_trajectories(); trajii++)
        {   clustbool=0;
            clustbool = cluster(trajii,current_time,primary);
            clustered[current_time](trajii,clustbool);
        }
        print_progress(timeii, sys->show_n_timesteps()-1);
    }



    cout << endl<< "Finding other members" << endl;cout.flush();

    for (int timeii=0; timeii<sys->show_n_timesteps();timeii++)
    {
        current_time  = convert_time(timeii);
        while (previous[current_time] != clustered[current_time])
        {
            previous[current_time] = clustered[current_time];
            for (int trajii=0;trajii<n_system_trajectories();trajii++)
            {
                clustbool=0;
                clustbool = cluster_correction(trajii,current_time, secondary, clustered[current_time]);
                clustered[current_time](trajii,clustbool||clustered[current_time](trajii));
            }
        }
//        while (previous2[current_time] != clustered[current_time])
//        {
//            previous2[current_time] = clustered[current_time];
//            for (int trajii=0;trajii<n_system_trajectories();trajii++)
//            {
//                clustbool=0;
//                clustbool = cluster_correction2(trajii,current_time,clustered[current_time]);
//                clustered[current_time](trajii,clustbool);
//            }
//        }
        print_progress(timeii, sys->show_n_timesteps()-1);
    }

    t_list->set(sys, n_times, n_system_trajectories(), clustered, time_conversion);





}


void Clustered_List::allocate_sig_matrix(string sig_file)
{
    string line;
    line = "";
    int sig_tokens=0;
    string * sig_ARGS;
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
        for(int argsii=1; argsii<=sys->show_n_atomtypes(); argsii++)
        {
            stringstream ss(sig_ARGS[argsii]);
            ss>>sigmatrix[0][argsii-1];
        }

        //input rest of the matrix
        for(int lineii=1;lineii<sys->show_n_atomtypes(); lineii++)
        {

            getline (file,line);
            sig_tokens = tokenize(line, sig_ARGS);

            species_name[lineii] = sig_ARGS[0];

            for(int argsii=1; argsii<=sys->show_n_atomtypes(); argsii++)
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





void Clustered_List::generate_neighbor_lists( int current_time,int trajii,string plane)
{

    if (included[current_time](trajii))
    {

        int neighborii=0;
        float sig_cut = 1.6;
        int jtype=0;
        float cutoff=0;
        Coordinate j_coord;
        Trajectory* current_trajectory;
        current_trajectory = sys->show_trajectory(trajii);
        jtype = current_trajectory->show_type()-1; // type should be mediated by system
        j_coord = current_trajectory->show_coordinate(current_time);

        Coordinate box_size;
        box_size = sys->size();


        typedef float (Coordinate::*length)()const;
        length distancefun;

        if(plane=="xy") distancefun = &Coordinate::length_xy;
        else if(plane=="xz") distancefun = &Coordinate::length_xz;
        else if(plane=="yz") distancefun = &Coordinate::length_yz;
        else
        {
            cout<<"Error: plane command "<<plane<<" not understood.\n";
            exit(1);
        }

        for (int traj2ii=0;traj2ii < n_system_trajectories();traj2ii++)
        {
            if (included[current_time](traj2ii))
            {
                int ktype=0;
                Coordinate k_coord;
                Trajectory* kth_trajectory;
                kth_trajectory=sys->show_trajectory(traj2ii);
                ktype=kth_trajectory->show_type()-1;
                k_coord = kth_trajectory->show_coordinate(current_time);
                Coordinate dist;
                float jk_distance = 0;

                dist = k_coord - j_coord;
                dist -= box_size * ((dist / (box_size*.5)).integer());


                jk_distance = (dist.*distancefun)();

                cutoff = sig_cut * sigmatrix[jtype][ktype];
                if (jk_distance < cutoff && jk_distance > 0)
                {
//                            cout <<"neighbor index = "<< neighborii << endl;cout.flush();
//                         cout << "neighID = "<< traj2ii << endl;cout.flush();
                    neighbor[current_time][trajii][neighborii]=traj2ii+1;
                 //cout << "neighbor["<<current_time<<"][";cout.flush();cout<<trajii;cout.flush();cout<<"]["<<neighborii;cout.flush();cout<<"] =" << neighbor[current_time][trajii][neighborii] << endl;cout.flush();
                    neighborii++;

                }
            }
        }

    }
}







bool Clustered_List::cluster( int trajii, int current_time,int primary)
{
    int neighborcount=0;
    bool clustbool=0;
        if (included[current_time](trajii))
        {
//            cout << "time = "<<current_time<< "\t";cout.flush();
//            cout << "traj = "<<trajii<< endl;cout.flush();
//            cout<< "neighbors of "<< trajii<<" at time "<<current_time<< endl;cout.flush();
//            for (int neighborii=0;neighborii<10;neighborii++)
//            {
//                cout<<neighbor[current_time][trajii][neighborii]<<"\t";cout.flush();
//            }


            for (int traj2ii=1;traj2ii < n_system_trajectories()+1;traj2ii++)
            {
                if (included[current_time](traj2ii-1))
                {
                    for (int neighborii=0;neighborii<10;neighborii++)
                    {
                        if (neighbor[current_time][trajii][neighborii] == traj2ii)
                        {
                            neighborcount++;
                        }
                    }
                }

            }

//cout << "\nneighbors = "<<neighborcount<< endl;cout.flush();
            if (neighborcount>=primary)
            {
                clustbool = 1;

            }
            else
            {
                clustbool = 0;
            }

        }



    return clustbool;

}


bool Clustered_List::cluster_correction( int trajii, int current_time,int secondary, Boolean_List clustered)
{
    int neighborcount=0;
    bool clustbool=0;
    if (included[current_time](trajii))
    {
    for (int traj2ii=1;traj2ii < n_system_trajectories()+1;traj2ii++)
        {
            if (clustered(traj2ii-1))
            {
                for (int neighborii=0;neighborii<10;neighborii++)
                {


                    if (neighbor[current_time][trajii][neighborii] == traj2ii)
                    {
                            neighborcount++;
                    }
                }
            }
        }


    if (neighborcount>=secondary)
    {
        clustbool = 1;
    }
    else
    {
        clustbool = 0;
    }

    }



    return clustbool;

}


bool Clustered_List::cluster_correction2( int trajii, int current_time, Boolean_List clustered)
{
    int neighborcount=0;
    bool clustbool=0;

    if (included[current_time](trajii))
    {
        for (int traj2ii=0;traj2ii < n_system_trajectories();traj2ii++)
        {
            if (clustered(traj2ii))
            {
                for (int neighborii=0;neighborii<10;neighborii++)
                {
                    if (neighbor[current_time][trajii][neighborii] == traj2ii)
                    {
                            neighborcount++;
                    }
                }
            }
        }



            if (neighborcount>=3)
            {
                clustbool = 1;

            }
            else
            {
                clustbool = 0;
            }


    }


    return clustbool;

}
