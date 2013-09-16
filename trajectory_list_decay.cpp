/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for trajectory list decay- determines how a trajectory list decays*/
/*Written by Daniel Hunsicker*/

#include "trajectory_list_decay.h"
#include "analysis.h"
#include <fstream>
#include "system.h"
#include <stdlib.h>
#include "trajectory_list.h"
#include "version.h"
#include "progress.h"




using namespace std;


Trajectory_List_Decay::Trajectory_List_Decay(System* sys)
{

    system = sys;

    n_blocks = system->show_n_exponentials();
    blocksize = system->show_n_exponential_steps();
    current_num = 0;
    initial_num = 0;
    current_time = 0;

    currently_included .set(system);
    initially_included.set(system);


    avg_decay = new float [blocksize];
    total_current = new int [blocksize];
    for (int stepii=0;stepii<blocksize;stepii++)
    {
        total_current[stepii] = 0;
        avg_decay[stepii] = 0;
    }

}


Trajectory_List_Decay::Trajectory_List_Decay()
{

    system = 0;

    n_blocks = 0;
    blocksize = 0;
    current_num = 0;
    initial_num = 0;
    current_time = 0;

//    currently_included = 0;
//    initially_included = 0;



    avg_decay = new float [0];
    total_current = new int [0];

}



Trajectory_List_Decay::Trajectory_List_Decay(const Trajectory_List_Decay & copy)
:Analysis(copy)
{

        system = copy.system;

        n_blocks = copy.n_blocks;
        blocksize = copy.blocksize;
        current_num = copy.current_num;
        initial_num = copy.initial_num;
        current_time = copy.current_time;
        total_initial = copy.total_initial;
        currently_included = copy.currently_included;
        initially_included = copy.initially_included;

        total_current = new int [blocksize];
        avg_decay = new float [blocksize];

        for (int stepii=0;stepii<blocksize;stepii++)
        {
            total_current[stepii] = copy.total_current[stepii];
            avg_decay[stepii] = copy.avg_decay[stepii];
        }

}

Trajectory_List_Decay Trajectory_List_Decay::operator=(const Trajectory_List_Decay & copy)
{
    if (this!=&copy)
    {
        delete [] avg_decay;
        delete [] total_current;

        system = copy.system;

        n_blocks = copy.n_blocks;
        blocksize = copy.blocksize;
        current_num = copy.current_num;
        initial_num = copy.initial_num;
        current_time = copy.current_time;
        total_initial = copy.total_initial;
        currently_included = copy.currently_included;
        initially_included = copy.initially_included;

        total_current = new int [blocksize];
        avg_decay = new float [blocksize];

        for (int stepii=0;stepii<blocksize;stepii++)
        {
            total_current[stepii] = copy.total_current[stepii];
            avg_decay[stepii] = copy.avg_decay[stepii];
        }
    }
      return *this;
}

Trajectory_List_Decay::~Trajectory_List_Decay()
{

        delete [] avg_decay;
        delete [] total_current;

}


void Trajectory_List_Decay::analyze(Trajectory_List* t_list)
{

    for (int blockii=0;blockii<n_blocks;blockii++)
    {
        initially_included = t_list->show_included(blockii*blocksize);
        initial_num = initially_included.show_n_included();
        total_initial += initial_num;
        for (int stepii=0;stepii<blocksize;stepii++)
        {
            current_time = blockii*blocksize+stepii;

            currently_included = t_list->show_included(current_time)&&initially_included;
            current_num = currently_included.show_n_included();

            total_current[stepii] += current_num;
        }
    }



    postprocess();
}


void Trajectory_List_Decay::listkernel(Trajectory* traj)
{

}

void Trajectory_List_Decay::postprocess()
{
    for (int stepii=0;stepii<blocksize;stepii++)
    {
        avg_decay[stepii] = float(total_current[stepii])/float(total_initial);
    }
}



void Trajectory_List_Decay::write(string filename)
{
    cout << "\nWriting trajectory list decay to file " << filename << "." << endl;
    ofstream output(filename.c_str());

    output << "Trajectory List Decay data created by AMDAT v." << VERSION << endl << endl;

    for (int stepii=0;stepii<blocksize;stepii++)
    {
        output << system->show_time(stepii);
        output << "\t" << avg_decay[stepii] << endl;
    }
}
