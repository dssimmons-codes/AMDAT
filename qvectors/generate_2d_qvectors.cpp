
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace std;

int main()
{
    float l;
    int binID;
    int ** x;
    int ** y;
    int n_bins;
    int * vectorID;
    int vector_allocation;
    int half_grid_length;

    half_grid_length = 100;
    vector_allocation =500;


    n_bins = ceil(((sqrt(2*pow(half_grid_length,2))+.25)/.5));

//    cout << "\n\n"<<n_bins<<"\n\n";cout.flush();

    vectorID = new int [n_bins+1];

    for (int binii=0; binii<=n_bins; binii++)
    {
    vectorID[binii]=0;
    }


    x = new int* [n_bins+1];
    y = new int* [n_bins+1];
    for (int binii=0; binii<=n_bins;binii++)
    {
        x[binii]= new int [vector_allocation];

        y[binii]= new int [vector_allocation];
    }

    for (int binii=0; binii<=n_bins;binii++)
    {
        for (int vectorii=0; vectorii<vector_allocation-1; vectorii++)
        {
        x[binii][vectorii] = 0;
        y[binii][vectorii] = 0;
        }
    }


/*-----------calculates and bins x and y values---------------*/


    for (int tempx=-half_grid_length; tempx<=half_grid_length; tempx++)
    {
        for (int tempy=0; tempy<=half_grid_length; tempy++)
        {

        l = sqrt(pow(tempx,2)+pow(tempy,2));

        binID = int(((l+.25)/.5));
//    cout << "\n"<<tempx<<"\t"<<tempy<<"\t"<<l<<"\t\t"<<binID<<"\n";cout.flush();


        x[binID][vectorID[binID]] = tempx;
        y[binID][vectorID[binID]] = tempy;

        vectorID[binID]++;
        }

    }




 /*-----------writes qvector files---------------*/




    for (int binii=0; binii<=n_bins; binii++)
    {
        string filename;
        string bin_str;
        stringstream bin;
        FILE * qvector;

        bin << binii;
        bin_str = bin.str();

        filename = "qvector.";
        if (binii<100)
        {
            filename = filename.append("0");
            if (binii<10)
            {
                filename = filename.append("0");
            }
        }

        filename = filename.append(bin_str);

        qvector = fopen (filename.c_str(),"w");

        for (int vectorii=0; vectorii<vectorID[binii];vectorii++)
        {
            fprintf(qvector, "\t%i\t%i\t0\n",x[binii][vectorii],y[binii][vectorii]);
        }
    }
}

