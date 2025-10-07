/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for class outputting time dependent information about system composition*/
/*Written by David Simmons*/


#include <stdlib.h>
#include "composition.h"
#include "composition_timedependent.h"
#include "version.h"
#include <omp.h>

using namespace std;



void Composition_TimeDependent::write(string filename)
{

    cout << "\nWriting composition to file " << filename << "." << endl;

    ofstream output(filename.c_str());

    output << "Composition data created by AMDAT v." << VERSION << endl;

    output << "time\t" << "number_density\t" << "number";

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        output << "\tspecies"<<typeii+1;
    }

    output << endl << "average\t" << average_density << "\t" << average_density*volume;

    for (int typeii=0;typeii<n_atomtypes;typeii++)
    {
        output << "\t" << time_average_comp[typeii];
    }

    output << endl;

    for (int timeii=0; timeii<n_times; timeii++)
    {
        output << system->show_time(timeii);
        output << "\t" << current_density[timeii];
        output << "\t" << current_density[timeii]*volume;
        for(int typeii=0;typeii<n_atomtypes;typeii++)
        {
            output << "\t" << time_average_comp[typeii];
        }
        output << endl;
    }

}

void Composition_TimeDependent::write(ofstream& output)
{
    output << "Composition data created by AMDAT v." << VERSION << endl;

    output << "time\t" << "number_density";

    for (int typeii=0; typeii<n_atomtypes;typeii++)
    {
        output << "\tspecies"<<typeii+1;
    }

    output << endl << "average\t" << average_density;

    for (int typeii=0;typeii<n_atomtypes;typeii++)
    {
        output << "\t" << time_average_comp[typeii];
    }

    output << endl;

    for(int timeii=0; timeii<n_times; timeii++)
    {
        output << system->show_time(timeii);
        output << "\t" << current_density[timeii];
        for (int typeii=0;typeii<n_atomtypes;typeii++)
        {
            output << "\t" << current_comp[timeii][typeii];
        }
        output << endl;
    }

}
