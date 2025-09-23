/*Trajectoy_List_Bins class - stores a list of trajectories according to bins
Amorphous Molecular Dynamics Tookit (AMDAT)
Written by Mark Mackura*/

#include "system.h"
#include "trajectory_list.h"
#include "trajectory_list_bins.h"
#include "coordinate.h"
#include "progress.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <float.h>
#include <cmath>
#include <algorithm>

using namespace std;

Trajectory_List_Bins::Trajectory_List_Bins()
{
  system=0;
  n_times=0;
  n_trajs=0;
  n_xbins=0;
  n_ybins=0;
  n_zbins=0;

  time_conversion = new int[n_times];
  xlo = new float[n_times];
  ylo = new float[n_times];
  zlo = new float[n_times];
  xhi = new float[n_times];
  yhi = new float[n_times];
  zhi = new float[n_times];


  for(int timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=timeii;
  }

  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }
}

Trajectory_List_Bins::~Trajectory_List_Bins()
{
//  cout << "Destructor called!" << endl;
  for(int xii=0; xii<n_xbins; xii++)
  {
    for(int yii=0; yii<n_ybins; yii++)
    {
      for(int zii=0; zii<n_zbins; zii++)
      {
        for(int timesii=0; timesii<n_times; timesii++)
        {
            delete [] include[xii][yii][zii][timesii];
            //delete [] trajcount[xii][yii][zii][timesii];
        }
        delete [] include[xii][yii][zii];
        delete [] trajcount[xii][yii][zii];
      }
       delete [] include[xii][yii];
       delete [] trajcount[xii][yii];
      //delete [] bins[xii][yii];
    }
    //delete [] bins[xii];
     delete [] include[xii];
     delete [] trajcount[xii];
  }
  delete [] xlo;
  delete [] ylo;
  delete [] zlo;
  delete [] xhi;
  delete [] yhi;
  delete [] zhi;
//  delete [] bins;
  delete [] time_conversion;
  delete [] include;
  delete [] trajcount;
//  cout << "Destructor complete!" << endl;
}

Trajectory_List_Bins::Trajectory_List_Bins(const Trajectory_List_Bins & copy)
{
  system=copy.system;
  n_times=copy.n_times;
  n_trajs=copy.n_trajs;
  n_xbins=copy.n_xbins;
  n_ybins=copy.n_ybins;
  n_zbins=copy.n_zbins;

  time_conversion = new int[n_times];
  xlo = new float[n_times];
  ylo = new float[n_times];
  zlo = new float[n_times];
  xhi = new float[n_times];
  yhi = new float[n_times];
  zhi = new float[n_times];
  //bins = new Trajectory_List ** [n_xbins];
//   include = new Boolean_List***[n_xbins];

  for(int timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=copy.time_conversion[timeii];
  }
  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }
}

Trajectory_List_Bins Trajectory_List_Bins::operator = (const Trajectory_List_Bins & copy)
{
  if(this!=&copy)
  {
    system=copy.system;
    n_times=copy.n_times;
    n_trajs=copy.n_trajs;
    n_xbins=copy.n_xbins;
    n_ybins=copy.n_ybins;
    n_zbins=copy.n_zbins;

//     for(int xii=0; xii<n_xbins; xii++)
//     {
//       for(int yii=0; yii<n_ybins; yii++)
//       {
// 	for(int zii=0; zii<n_zbins; zii++)
// 	{
// 	    delete [] include[xii][yii][zii];
// 	}
// 	delete [] include[xii][yii];
// 	delete [] bins[xii][yii];
//       }
//       delete [] bins[xii];
//       delete [] include[xii];
//     }
    delete [] xlo;
    delete [] ylo;
    delete [] zlo;
    delete [] xhi;
    delete [] yhi;
    delete [] zhi;
    //delete [] bins;
    delete [] time_conversion;
//     delete [] include;

    for(int timeii=0;timeii<n_times;timeii++)
    {
      time_conversion[timeii]=copy.time_conversion[timeii];
    }
  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }
  }
  return *this;
}

Trajectory_List_Bins::Trajectory_List_Bins(System * sys,int xbins,int ybins,int zbins)
{
    /** Bins all system trajectories according to entire system dimensions
    * @param sys - pointer to System object to be binned
    * @param xbins-ybins-zbins - number of bins in each direction
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;                                      //system to bin
    //int xii, yii, zii,
    int timeii;                 	       //indicies
    n_times = system->show_n_timesteps();              //number of system times
    n_trajs = system->show_n_trajectories();           //number of trajectories to bin
    n_xbins = xbins;
    n_ybins = ybins;
    n_zbins = zbins;

    Coordinate const ** boxsize;
    boxsize = new const Coordinate*[n_times];
    xlo = new float[n_times];
    ylo = new float[n_times];
    zlo = new float[n_times];
    xhi = new float[n_times];
    yhi = new float[n_times];
    zhi = new float[n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
      boxsize[timeii] = system->boundaries(timeii);
      /* Add smallest possible value to boundaries so particles that fall directly on boundary are handled properly */
      xlo[timeii] = boxsize[timeii][0](0);
      ylo[timeii] = boxsize[timeii][0](1);
      zlo[timeii] = boxsize[timeii][0](2);
      xhi[timeii] = boxsize[timeii][1](0);
      yhi[timeii] = boxsize[timeii][1](1);
      zhi[timeii] = boxsize[timeii][1](2);

    }

    /* Create time_conversion array with each value same as index*/

    time_conversion = new int[n_times];
    for(timeii=0;timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
    }

    /*Allocate memory for trajectory and boolean arrays for each bin*/


    
  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }
   

   assign_bins();

}
Trajectory_List_Bins::Trajectory_List_Bins(System * sys,int xbins,int ybins,int zbins,float x_lo,float x_hi,float y_lo,float y_hi,float z_lo, float z_hi)
{
    /** Bins all system trajectories only in box specified by high and low positions
    * @param sys - pointer to system object to be binned
    * @param xbins-ybins-zbins - number of bins in each direction
    * @param x_lo-x_hi-y_lo-y_hi-z_lo-z_hi - dimensions of region to bin
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;                                      //system to bin
   // int xii, yii, zii,
    int timeii;                  //indicies
    n_times = system->show_n_timesteps();              //number of system times
    n_trajs = system->show_n_trajectories();           //number of trajectories to bin
    n_xbins = xbins;
    n_ybins = ybins;
    n_zbins = zbins;
    cout<<"\nNumber of bins (x,y,z): ("<<n_xbins<<","<<n_ybins<<","<<n_zbins<<")";

    Coordinate const ** boxsize;
    boxsize = new const Coordinate*[n_times];
    xlo = new float[n_times];
    ylo = new float[n_times];
    zlo = new float[n_times];
    xhi = new float[n_times];
    yhi = new float[n_times];
    zhi = new float[n_times];

    for (int timeii=0; timeii<n_times;timeii++)
    {
      boxsize[timeii] = system->boundaries(timeii);
      /* Add smallest possible value to boundaries if necessary so particles that fall directly on boundary are included */
      /* Check which dimensions are full box sizes, and assign region dimensions if not full size */
      if ( x_lo == 0 || x_lo<boxsize[timeii][0](0))
      {
	xlo[timeii] = boxsize[timeii][0](0); //binning boundary set to system boundary @ timeii
      }
      else
      {
	xlo[timeii] = x_lo; //binning boundary set to the user defined value @ timeii
      }
      if ( y_lo == 0 || y_lo<boxsize[timeii][0](1))
      {
	ylo[timeii] = boxsize[timeii][0](1);
      }
      else
      {
	ylo[timeii] = y_lo;
      }
      if ( z_lo == 0 || z_lo<boxsize[timeii][0](2))
      {
	zlo[timeii] = boxsize[timeii][0](2);
      }
      else
      {
	zlo[timeii] = z_lo;
      }
      if ( x_hi == 0 || x_hi>boxsize[timeii][1](0))
      {
	xhi[timeii] = boxsize[timeii][1](0);
      }
      else
      {
	xhi[timeii] = x_hi;
      }
      if ( y_hi == 0 || y_hi>boxsize[timeii][1](1))
      {
	yhi[timeii] = boxsize[timeii][1](1);
      }
      else
      {
	yhi[timeii] = y_hi;
      }
      if ( z_hi == 0 || z_hi>boxsize[timeii][1](2))
      {
	zhi[timeii] = boxsize[timeii][1](2);
      }
      else
      {
	zhi[timeii] = z_hi;
      }
    }

    /* Create time_conversion array with each value same as index*/

    time_conversion = new int[n_times];
    for(timeii=0;timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
    }
    //cout<<"\nTimes: "<<n_times<<"\t Trajectories: "<<n_trajs;cout.flush();

    /*Allocate memory for trajectory and boolean arrays for each bin*/

  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

   assign_bins();

}
Trajectory_List_Bins::Trajectory_List_Bins(System * sys,int xbins,int ybins,int zbins,Coordinate los,Coordinate his)
{
    /** Bins all system trajectories only in box specified two coordinates
    * @NOTE METHOD IS OUT OF DATE - this method is not current and needs to be revised to function as the 'float' based boundary constructor, with varying boxsize
    * @param los - Coordinate object holding lower bounds for box dimensions
    * @param his - Coordinate object holding upper bounds for box dimensions
    * @param xbins-ybins-zbins - number of bins in each direction
    * @param sys - pointer to System object to be binned
    * @author Mark Mackura
    * @date 4/12/2012
    **/
//     system = sys;                                      //system to bin
//     int xii, yii, zii, timeii;                  //indicies
//     n_times = system->show_n_timesteps();              //number of system times
//     n_trajs = system->show_n_trajectories();           //number of trajectories to bin
//     n_xbins = xbins;
//     n_ybins = ybins;
//     n_zbins = zbins;
    cout<<"\n Error: Coordinate constructor for binning is incomplete, see javadoc comment."<<endl;
    cout.flush();
    exit(1);


    /* Define boundaries of bins*/
//     xlo = los(0);
//     ylo = los(1);
//     zlo = los(2);
//     xhi = his(0);
//     yhi = his(1);
//     zhi = his(2);
//
//     /* Create time_conversion array with each value same as index*/
//
//     time_conversion = new int[n_times];
//     for(timeii=0;timeii<n_times;timeii++)
//     {
//         time_conversion[timeii]=timeii;
//     }
//
//     /*Allocate memory for trajectory and boolean arrays for each bin*/
//
  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

    assign_bins();


}



Trajectory_List_Bins::Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, Trajectory_List* clustered_list)
{
    /** Bins all system trajectories according to entire system dimensions
    * @param sys - pointer to System object to be binned
    * @param xbins-ybins-zbins - number of bins in each direction
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;                                      //system to bin
    n_times = system->show_n_timesteps();              //number of system times
    n_trajs = system->show_n_trajectories();           //number of trajectories to bin
    n_xbins = 1;
    n_ybins = 1;

    Coordinate const ** boxsize;
    boxsize = new const Coordinate*[n_times];
    xlo = new float[n_times];
    ylo = new float[n_times];
    zlo = new float[n_times];
    xhi = new float[n_times];
    yhi = new float[n_times];
    zhi = new float[n_times];
    float length =0;
    for (int timeii=0; timeii<n_times;timeii++)
    {
        float x_length,y_length,z_length;
      boxsize[timeii] = system->boundaries(timeii);
      /* Add smallest possible value to boundaries so particles that fall directly on boundary are handled properly */
      xlo[timeii] = boxsize[timeii][0](0);
      ylo[timeii] = boxsize[timeii][0](1);
      zlo[timeii] = boxsize[timeii][0](2);
      xhi[timeii] = boxsize[timeii][1](0);
      yhi[timeii] = boxsize[timeii][1](1);
      zhi[timeii] = boxsize[timeii][1](2);
      x_length = (xhi[timeii] - xlo[timeii])/2;
      y_length = (yhi[timeii] - ylo[timeii])/2;
      z_length = (zhi[timeii] - zlo[timeii])/2;
      if (pow(pow(x_length,2)+pow(y_length,2)+pow(z_length,2),.5)>length)
      {
          length = pow(pow(x_length,2)+pow(y_length,2)+pow(z_length,2),.5);
      }
    }


    int max_zbins = int(length/bin_thickness)+1;
    if (n_bins < max_zbins)
    {
        n_zbins = n_bins;
    }
    else
    {
        cout<< endl<< "system too small for "<< n_bins <<" bins. using " <<max_zbins <<" instead."<<endl; cout.flush();
        n_zbins = max_zbins;
    }

 for (int timeii=0; timeii<n_times;timeii++)
    {
      /* set boundaries to zero to avoid wrong volume calculations from methods (composition)*/
      xlo[timeii] = 0;
      ylo[timeii] = 0;
      zlo[timeii] = 0;
      xhi[timeii] = 0;
      yhi[timeii] = 0;
      zhi[timeii] = 0;
    }

    /* Create time_conversion array with each value same as index*/

    time_conversion = new int[n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
    }

    /*Allocate memory for trajectory and boolean arrays for each bin*/


  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

   assign_bins_distance_clusters(binning_list,clustered_list,bin_thickness);

}

Trajectory_List_Bins::Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, string pln, float posit, string dir)
{
    /** Bins all system trajectories according to entire system dimensions
    * @param sys - pointer to System object to be binned
    * @param xbins-ybins-zbins - number of bins in each direction
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;                                      //system to bin
    n_times = system->show_n_timesteps();              //number of system times
    n_trajs = system->show_n_trajectories();           //number of trajectories to bin


    Coordinate const ** boxsize;
    boxsize = new const Coordinate*[n_times];
    xlo = new float[n_times];
    ylo = new float[n_times];
    zlo = new float[n_times];
    xhi = new float[n_times];
    yhi = new float[n_times];
    zhi = new float[n_times];
    float x_length,y_length,z_length;
    for (int timeii=0; timeii<n_times;timeii++)
    {

      boxsize[timeii] = system->boundaries(timeii);
      /* Add smallest possible value to boundaries so particles that fall directly on boundary are handled properly */

          xlo[timeii] = boxsize[timeii][0](0);
          ylo[timeii] = boxsize[timeii][0](1);
          zlo[timeii] = boxsize[timeii][0](2);
          xhi[timeii] = boxsize[timeii][1](0);
          yhi[timeii] = boxsize[timeii][1](1);
          zhi[timeii] = boxsize[timeii][1](2);
      if (dir=="average")
      {
            x_length = (xhi[timeii] - xlo[timeii])/2;
            y_length = (yhi[timeii] - ylo[timeii])/2;
            z_length = (zhi[timeii] - zlo[timeii])/2;
      }
      else if (dir=="above" || dir=="below")
      {
            x_length = (xhi[timeii] - xlo[timeii]);
            y_length = (yhi[timeii] - ylo[timeii]);
            z_length = (zhi[timeii] - zlo[timeii]);
      }
        else
        {
           cout << "could not recognize direction " << dir << ". valid choices are above, below or average.";cout.flush();
        }
    }



    Coordinate pos;
    if (pln == "x")
    {
        int max_bins = int(x_length/bin_thickness)+1;
        if (n_bins < max_bins)
        {
            n_xbins = n_bins;
        }
        else
        {
            cout<< endl<< "system too small for "<< n_bins <<" bins. Using " <<max_bins<< " instead."<<endl; cout.flush();
            n_xbins = max_bins;
        }
        n_ybins = 1;
        n_zbins = 1;
	    for (int timeii=0; timeii<n_times;timeii++)
    {

      if (dir=="average")
      {
            xlo[timeii] = posit - bin_thickness*n_xbins;
            xhi[timeii] = posit + bin_thickness*n_xbins;
      }
      else if (dir=="above")
      {
            xlo[timeii] = posit;
            xhi[timeii] = posit + bin_thickness*n_xbins;
      }
      else if (dir=="below")
      {
            xlo[timeii] = posit - bin_thickness*n_xbins;
            xhi[timeii] = posit;
      }
    }
    }
    else if (pln == "y")
    {
        n_xbins = 1;
        int max_bins = int(y_length/bin_thickness)+1;
        if (n_bins < max_bins)
        {
            n_ybins = n_bins;
        }
        else
        {
            cout<< endl<< "system too small for " <<n_bins <<" bins. Using "<< max_bins <<" instead."<<endl; cout.flush();
            n_ybins = max_bins;
        }
        n_zbins = 1;
        pos.set(0,posit,0);
	    	    for (int timeii=0; timeii<n_times;timeii++)
    {

      if (dir=="average")
      {
            ylo[timeii] = posit - bin_thickness*n_ybins;
            yhi[timeii] = posit + bin_thickness*n_ybins;
      }
      else if (dir=="above")
      {
            ylo[timeii] = posit;
            yhi[timeii] = posit + bin_thickness*n_ybins;
      }
      else if (dir=="below")
      {
            ylo[timeii] = posit - bin_thickness*n_ybins;
            yhi[timeii] = posit;
      }
    }
    }
    else if (pln == "z")
    {
        n_xbins = 1;
        n_ybins = 1;
        int max_bins = int(z_length/bin_thickness)+1;
            if (n_bins < max_bins)
        {
            n_zbins = n_bins;
        }
        else
        {
            cout<< endl<<  "system too small for " <<n_bins<< " bins. Using " <<max_bins <<" instead."<<endl; cout.flush();
            n_zbins = max_bins;
        }
        pos.set(0,0,posit);
	    	    for (int timeii=0; timeii<n_times;timeii++)
    {

      if (dir=="average")
      {
            zlo[timeii] = posit - bin_thickness*n_zbins;
            zhi[timeii] = posit + bin_thickness*n_zbins;
      }
      else if (dir=="above")
      {
            zlo[timeii] = posit;
            zhi[timeii] = posit + bin_thickness*n_zbins;
      }
      else if (dir=="below")
      {
            zlo[timeii] = posit - bin_thickness*n_zbins;
            zhi[timeii] = posit;
      }
    }
    }
    else
    {
        cout << "plane "<<pln<<" not recognized, plane must be x, y, or z." <<endl;cout.flush();
        exit(1);
    }





//cout << n_xbins<< "\t"<< n_ybins<< "\t"<< n_zbins<< endl; cout.flush();

    /* Create time_conversion array with each value same as index*/

    time_conversion = new int[n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
    }

    /*Allocate memory for trajectory and boolean arrays for each bin*/


  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

   assign_bins_distance_plane(binning_list,pln,posit,dir,bin_thickness);

}



Trajectory_List_Bins::Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, Coordinate pnt)
{
    /** Bins all system trajectories according to entire system dimensions
    * @param sys - pointer to System object to be binned
    * @param xbins-ybins-zbins - number of bins in each direction
    * @author Mark Mackura
    * @date 4/12/2012
    **/
    system = sys;                                      //system to bin
    n_times = system->show_n_timesteps();              //number of system times
    n_trajs = system->show_n_trajectories();           //number of trajectories to bin
    n_xbins = 1;
    n_ybins = 1;

    Coordinate const ** boxsize;
    boxsize = new const Coordinate*[n_times];
    xlo = new float[n_times];
    ylo = new float[n_times];
    zlo = new float[n_times];
    xhi = new float[n_times];
    yhi = new float[n_times];
    zhi = new float[n_times];
    float length =0;
    for (int timeii=0; timeii<n_times;timeii++)
    {
        float x_length,y_length,z_length;
      boxsize[timeii] = system->boundaries(timeii);
      /* Add smallest possible value to boundaries so particles that fall directly on boundary are handled properly */
      xlo[timeii] = boxsize[timeii][0](0);
      ylo[timeii] = boxsize[timeii][0](1);
      zlo[timeii] = boxsize[timeii][0](2);
      xhi[timeii] = boxsize[timeii][1](0);
      yhi[timeii] = boxsize[timeii][1](1);
      zhi[timeii] = boxsize[timeii][1](2);
      x_length = (xhi[timeii] - xlo[timeii])/2;
      y_length = (yhi[timeii] - ylo[timeii])/2;
      z_length = (zhi[timeii] - zlo[timeii])/2;
      if (pow(pow(x_length,2)+pow(y_length,2)+pow(z_length,2),.5)>length)
      {
          length = pow(pow(x_length,2)+pow(y_length,2)+pow(z_length,2),.5);
      }
    }


    int max_zbins = int(length/bin_thickness)+1;
    if (n_bins < max_zbins)
    {
        n_zbins = n_bins;
    }
    else
    {
        cout<< endl<< "system too small for "<< n_bins <<" bins. using " <<max_zbins <<" instead."<<endl; cout.flush();
        n_zbins = max_zbins;
    }
for (int timeii=0; timeii<n_times;timeii++)
    {
      /* set boundaries to zero to avoid wrong volume calculations from methods (composition)*/
      xlo[timeii] = 0;
      ylo[timeii] = 0;
      zlo[timeii] = 0;
      xhi[timeii] = 0;
      yhi[timeii] = 0;
      zhi[timeii] = 0;
    }

    /* Create time_conversion array with each value same as index*/

    time_conversion = new int[n_times];
    for(int timeii=0;timeii<n_times;timeii++)
    {
        time_conversion[timeii]=timeii;
    }

    /*Allocate memory for trajectory and boolean arrays for each bin*/


  include = new int****[n_xbins];
   trajcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     include[xii] = new int***[n_ybins];
     trajcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
 	include[xii][yii] = new int**[n_zbins];
	trajcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
 	    include[xii][yii][zii] = new int*[n_times];
	    trajcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

   assign_bins_distance_point(binning_list,pnt,bin_thickness);

}
Trajectory_List Trajectory_List_Bins::persistant_trajectories(int xii,int yii,int zii,int starttime,int endtime)
{
    /** Returns list of particles in bin at single timestep (start == end) or for last two timesteps (end-1 && end)
    * @param xii-yii-zii - bin indicies
    * @param starttime-endtime - time interval in terms of frames
    * @author Mark Mackura
    * @date 4/12/2012
    **/

    Boolean_List * bins;

    bins = new Boolean_List[n_times];

    bin_boolean_calculation(xii,yii,zii,bins);

    Trajectory_List bin = Trajectory_List(system,n_times,n_trajs,bins,time_conversion);
    if (starttime == endtime || starttime == 0)
    {
	return bin.time_intersection(starttime,endtime); //return Trajectory_List of particles in bin only at starttime or from starttime=0 to endtime
    }
    else
    {
        return bin.time_intersection(endtime-1,endtime); //return Trajectory_List of particles who are in bin for last two times in interval
    }
}


void Trajectory_List_Bins::bin_boolean_calculation(int xii,int yii,int zii, Boolean_List * bins)
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
    bins[timeii].set(system,include[xii][yii][zii][timeii],trajcount[xii][yii][zii][timeii]);
  }
}



Trajectory_List Trajectory_List_Bins::operator()(int xii,int yii,int zii)
{
    /** Returns Trajectory_List for specific bin
    *@param xii-yii-zii - bin indicies
    */
    Boolean_List * bins;

    bins = new Boolean_List[n_times];

    bin_boolean_calculation(xii,yii,zii,bins);

    Trajectory_List templist(system,n_times,n_trajs,bins,time_conversion);

    delete [] bins;

    return templist;
}





int Trajectory_List_Bins::show_n_xbins()
{
  return n_xbins;
}
int Trajectory_List_Bins::show_n_ybins()
{
  return n_ybins;
}
int Trajectory_List_Bins::show_n_zbins()
{
   return n_zbins;
}
float Trajectory_List_Bins::show_lx(int timeii)
{
  return xhi[timeii]-xlo[timeii];
}
float Trajectory_List_Bins::show_ly(int timeii)
{
  return yhi[timeii]-ylo[timeii];
}
float Trajectory_List_Bins::show_lz(int timeii)
{
  return zhi[timeii]-zlo[timeii];
}



void Trajectory_List_Bins::assign_bins()
{
  /**Rectangularly bins all trajs at all times within given dimensions (each bin has it's own Boolean_List)
     * @author Mark Mackura and David Simmons
     * @date 3/1/2013
     */
  int timeii,trajii,xii,yii,zii;
  int **** tempcount;
  float xcoord,ycoord,zcoord;
  int *** temp_bins;

  temp_bins = new int ** [n_trajs];
  for(trajii=0; trajii<n_trajs; trajii++)
  {
    temp_bins[trajii]=new int * [n_times];
    for(timeii=0; timeii<n_times; timeii++)
    {
      temp_bins[trajii][timeii]=new int [3];
    }
  }
  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
	{
	  for(timeii=0; timeii<n_times; timeii++)
	  {
	    trajcount[xii][yii][zii][timeii]=0;
	  }
	}
     }
   }

  for(trajii=0; trajii<n_trajs; trajii++)
  {
    for(timeii=0; timeii<n_times; timeii++)
    {
      xcoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_x();
      ycoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_y();
      zcoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_z();
      if((xcoord<xhi[timeii] && xcoord>xlo[timeii]) && (ycoord<yhi[timeii] && ycoord>ylo[timeii]) && (zcoord<zhi[timeii] && zcoord>zlo[timeii]))
      {
	xii = int((xcoord-xlo[timeii]+FLT_EPSILON*10.0)/(xhi[timeii]-xlo[timeii]+FLT_EPSILON*20.0)*float(n_xbins));
	yii = int((ycoord-ylo[timeii]+FLT_EPSILON*10.0)/(yhi[timeii]-ylo[timeii]+FLT_EPSILON*20.0)*float(n_ybins));
	zii = int((zcoord-zlo[timeii]+FLT_EPSILON*10.0)/(zhi[timeii]-zlo[timeii]+FLT_EPSILON*20.0)*float(n_zbins));

	temp_bins[trajii][timeii][0]=xii;
	temp_bins[trajii][timeii][1]=yii;
	temp_bins[trajii][timeii][2]=zii;

// 	                               cout <<xii<<"\t"; cout.flush();
//                    cout <<yii<<"\t"; cout.flush();
//                    cout <<zii<<"\t"; cout.flush();
//                    cout << timeii<<"\t"; cout.flush();
                   
	
	trajcount[xii][yii][zii][timeii]++;
      }
    }
  }

   tempcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     tempcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
	tempcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	    tempcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	  for(timeii=0; timeii<n_times; timeii++)
	  {
	    include[xii][yii][zii][timeii]=new int [trajcount[xii][yii][zii][timeii]];
	    tempcount[xii][yii][zii][timeii]=0;
	  }
	}
     }
   }


  for(trajii=0; trajii<n_trajs; trajii++)
  {
    for(timeii=0; timeii<n_times; timeii++)
    {
      xcoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_x();
      ycoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_y();
      zcoord = system->show_trajectory(trajii)->show_coordinate(timeii).show_z();
      if((xcoord<xhi[timeii] && xcoord>xlo[timeii]) && (ycoord<yhi[timeii] && ycoord>ylo[timeii]) && (zcoord<zhi[timeii] && zcoord>zlo[timeii]))
      {
      xii=temp_bins[trajii][timeii][0];
      yii=temp_bins[trajii][timeii][1];
      zii=temp_bins[trajii][timeii][2];

      if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii])
      {
	cout<<"\n"<<xii<<"\t"<<yii<<"\t"<<zii<<"\t"<<timeii<<"\t"<<trajii<<"\t"<<n_times<<"\t";cout.flush();
      }
      if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii])
      {
	cout<<tempcount[xii][yii][zii][timeii]<<"\t"<<include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]<<"\t"<<trajcount[xii][yii][zii][timeii]<<"\t";cout.flush();
	
      }

      include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]=trajii;
      tempcount[xii][yii][zii][timeii]++;
      }
    }
  }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	  delete [] tempcount[xii][yii][zii];
	}
	delete [] tempcount[xii][yii];
     }
     delete [] tempcount[xii];
   }
   delete [] tempcount;


  for(trajii=0; trajii<n_trajs; trajii++)
  {
    for(timeii=0; timeii<n_times; timeii++)
    {
      delete [] temp_bins[trajii][timeii];
    }
    delete [] temp_bins[trajii];
  }
  delete [] temp_bins;

}

void Trajectory_List_Bins::assign_bins_distance_clusters(Trajectory_List * binned_list,Trajectory_List * cluster_list,float bin_thickness)
{
  /**bin based on distance to a second trajectory list (each bin has it's own Boolean_List)
     * @author Mark Mackura and David Simmons
     * @date 3/1/2013
     */
  int zii;

  float length;
  float temp_length;
  
  for(int xii=0; xii<n_xbins; xii++)
  {
    for(int yii=0; yii<n_ybins; yii++)
    {
      for(int zii=0; zii<n_zbins; zii++)
      {
        for(int timeii=0; timeii<n_times; timeii++)
        {
          trajcount[xii][yii][zii][timeii]=0;
        }
      }
    }
  }

  vectorstoragething.resize(n_xbins,vector<vector<vector<vector<int>>>>(n_ybins,vector<vector<vector<int>>>(n_zbins,vector<vector<int>>(n_times,vector<int>(0)))));
  for(int timeii=0; timeii<n_times; timeii++)
  {
    for(int trajii=0; trajii<binned_list->show_n_trajectories(timeii); trajii++)
    {
      length = system->size(timeii).max();
      Coordinate coordii = (*binned_list)(timeii,trajii)->show_coordinate(timeii);
      int n_traj2s = cluster_list->show_n_trajectories(timeii);
      for(int traj2ii=0; traj2ii<n_traj2s; traj2ii++)
      {
        temp_length=
                (
                  (
                    coordii -
                    (*cluster_list)(timeii,traj2ii)->show_coordinate(timeii)
                  ).vector_unwrapped(system->size())
                ).length();

        length=min(length,temp_length);
      }
      zii = int(length/bin_thickness);
      if (zii<n_zbins)
      {
        trajcount[0][0][zii][timeii]++;
        vectorstoragething[0][0][zii][timeii].push_back((*binned_list)(timeii,trajii)->show_trajectory_ID());
      }
    }
  }

  for(int xii=0; xii<n_xbins; xii++)
  {
    for(int yii=0; yii<n_ybins; yii++)
    {
      for(int zii=0; zii<n_zbins; zii++)
      {
        for(int timeii=0; timeii<n_times; timeii++)
        {
          /*allocate memory for include IDs at each x,y,z,t; then loop over included IDs in vectorstoragething[xii,yii,zii,timeii] and then copy them to include[xii][yii][zii][timeii][ii]*/
	        include[xii][yii][zii][timeii] = new int [vectorstoragething[xii][yii][zii][timeii].size()];
//          include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]=trajii;
          for (int trajii=0; trajii < vectorstoragething[xii][yii][zii][timeii].size(); trajii++)
          {
            include[xii][yii][zii][timeii][trajii] = vectorstoragething[xii][yii][zii][timeii][trajii];
          }
        }
      }
    }
  }

//   /*loops after here are old - delete*/
//   
//   
//  for(int trajii=0; trajii<n_trajs; trajii++)
//  {
//    for(int timeii=0; timeii<n_times; timeii++)
//    {
////cout << binned_list->is_included(timeii,trajii)<<"\t";cout.flush();
//          Coordinate box_size(system->boundaries(timeii)[1](0)-system->boundaries(timeii)[0](0),system->boundaries(timeii)[1](1)-system->boundaries(timeii)[0](1),system->boundaries(timeii)[1](2)-system->boundaries(timeii)[0](2));
//          if (binned_list->is_included(timeii,trajii))
//          {
////              cout<<"traj "<< trajii <<" at time "<< timeii<<" is included in binned list"<<endl;cout.flush();
////              cout << cluster_list->show_n_trajectories(timeii)<< endl;cout.flush();
//
//           Coordinate dist;
//           float length=100000;
//            for(int traj2ii=0; traj2ii<n_trajs; traj2ii++)
//            {//cout << cluster_list->is_included(timeii,trajii)<<"\t";cout.flush();
//                if (cluster_list->is_included(timeii,traj2ii))
//                {
//                    //cout <<"traj "<< traj2ii<<" is included in clustered list"<<endl;cout.flush();
//                    float temp_length;
//                    dist = system->show_trajectory(trajii)->show_coordinate(timeii) - system->show_trajectory(traj2ii)->show_coordinate(timeii);
//                    dist -= box_size * ((dist/(box_size*.5)).integer());
//
//                    //cout <<dist.show_x()<<"\t"<<dist.show_y()<<"\t"<<dist.show_z()<<endl; cout.flush();
//                    temp_length = dist.length();
//
//                    if (temp_length <length)
//                    {
//                        length=temp_length;
//                    }
//                }
//            }
//
//
//
//
//            xii = 0;
//            yii =0;
//            zii = int(length/bin_thickness);
//
//                    temp_bins[trajii][timeii][0]=xii;
//                    temp_bins[trajii][timeii][1]=yii;
//                    temp_bins[trajii][timeii][2]=zii;
//                      if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
//                {
//                    //            cout <<xii<<"\t"; cout.flush();
//        //            cout <<yii<<"\t"; cout.flush();
//        //            cout <<zii<<"\t"; cout.flush();
//        //            cout << timeii<<"\t"; cout.flush();
//        //            cout <<trajcount[xii][yii][zii][timeii]<<endl; cout.flush();
//
//
//                    trajcount[xii][yii][zii][timeii]++;
//                }
//            }
//    }
//}
//
//   tempcount = new int***[n_xbins];
//   for(int xii=0; xii<n_xbins; xii++)
//   {
//     tempcount[xii] = new int**[n_ybins];
//     for(int yii=0; yii<n_ybins; yii++)
//     {
//	tempcount[xii][yii] = new int*[n_zbins];
// 	for(int zii=0; zii<n_zbins; zii++)
// 	{
//	    tempcount[xii][yii][zii] = new int[n_times];
// 	}
//     }
//   }
//
//  for(int xii=0; xii<n_xbins; xii++)
//   {
//     for(int yii=0; yii<n_ybins; yii++)
//     {
// 	for(int zii=0; zii<n_zbins; zii++)
// 	{
//	  for(int timeii=0; timeii<n_times; timeii++)
//	  {
//	    include[xii][yii][zii][timeii]=new int [trajcount[xii][yii][zii][timeii]];
//	    tempcount[xii][yii][zii][timeii]=0;
//	  }
//	}
//     }
//   }
//
//
//  for(int trajii=0; trajii<n_trajs; trajii++)
//  {
//
//    for(int timeii=0; timeii<n_times; timeii++)
//    {
//
//          if (binned_list->is_included(timeii,trajii))
//          {
//              xii=temp_bins[trajii][timeii][0];
//              yii=temp_bins[trajii][timeii][1];
//              zii=temp_bins[trajii][timeii][2];
//                if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
//                {
//                     if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii]){cout<<"\n"<<xii<<"\t"<<yii<<"\t"<<zii<<"\t"<<timeii<<"\t"<<trajii<<"\t"<<n_times<<"\t";cout.flush();}
//
//                     if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii])
//                     {cout<<tempcount[xii][yii][zii][timeii]<<"\t"<<include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]<<"\t"<<trajcount[xii][yii][zii][timeii]<<"\t";cout.flush();}
//
//                     include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]=trajii;
//                     
//                     tempcount[xii][yii][zii][timeii]++;
//                }
//
//          }
//    }
//  }
//
//  for(int xii=0; xii<n_xbins; xii++)
//   {
//     for(int yii=0; yii<n_ybins; yii++)
//     {
// 	for(int zii=0; zii<n_zbins; zii++)
// 	{
//	  delete [] tempcount[xii][yii][zii];
//	}
//	delete [] tempcount[xii][yii];
//     }
//     delete [] tempcount[xii];
//   }
//   delete [] tempcount;
//  
//  for(int trajii=0; trajii<n_trajs; trajii++)
//  {
//    for(int timeii=0; timeii<n_times; timeii++)
//    {
//      delete [] temp_bins[trajii][timeii];
//    }
//    delete [] temp_bins[trajii];
//  }
//  delete [] temp_bins;

}


void Trajectory_List_Bins::assign_bins_distance_point(Trajectory_List * binned_list,Coordinate point,float bin_thickness)
{
  /**Rectangularly bins all trajs at all times within given dimensions (each bin has it's own Boolean_List)
     * @author Mark Mackura and David Simmons
     * @date 3/1/2013
     */
  int xii,yii,zii;
  int **** tempcount;
  int *** temp_bins;

  temp_bins = new int ** [n_trajs];
  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    temp_bins[trajii]=new int * [n_times];
    for(int timeii=0; timeii<n_times; timeii++)
    {
      temp_bins[trajii][timeii]=new int [3];
    }
  }
  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
	{
	  for(int timeii=0; timeii<n_times; timeii++)
	  {
	    trajcount[xii][yii][zii][timeii]=0;
	  }
	}
     }
   }



  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {
          Coordinate box_size(system->boundaries(timeii)[1](0)-system->boundaries(timeii)[0](0),system->boundaries(timeii)[1](1)-system->boundaries(timeii)[0](1),system->boundaries(timeii)[1](2)-system->boundaries(timeii)[0](2));
          if (point > box_size)
          {
              cout << "point chosen does not lie withing the box dimmensions." << endl;cout.flush();
              exit(1);
          }
          if (binned_list->is_included(timeii,trajii))
          {

           Coordinate dist;
           float length;
                    dist = system->show_trajectory(trajii)->show_coordinate(timeii) - point;
                    dist -= box_size * ((dist/(box_size*.5)).integer());

                    length = dist.length();







            xii = 0;
            yii =0;
            zii = int(length/bin_thickness);





                    temp_bins[trajii][timeii][0]=xii;
                    temp_bins[trajii][timeii][1]=yii;
                    temp_bins[trajii][timeii][2]=zii;
                if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
                {
                    //            cout <<xii<<"\t"; cout.flush();
        //            cout <<yii<<"\t"; cout.flush();
        //            cout <<zii<<"\t"; cout.flush();
        //            cout << timeii<<"\t"; cout.flush();
        //            cout <<trajcount[xii][yii][zii][timeii]<<endl; cout.flush();


                    trajcount[xii][yii][zii][timeii]++;
                }

            }
    }
}

   tempcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     tempcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
	tempcount[xii][yii] = new int*[n_zbins];
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	    tempcount[xii][yii][zii] = new int[n_times];
 	}
     }
   }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	  for(int timeii=0; timeii<n_times; timeii++)
	  {
	    include[xii][yii][zii][timeii]=new int [trajcount[xii][yii][zii][timeii]];
	    tempcount[xii][yii][zii][timeii]=0;
	  }
	}
     }
   }


  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {

          if (binned_list->is_included(timeii,trajii))
          {
              xii=temp_bins[trajii][timeii][0];
              yii=temp_bins[trajii][timeii][1];
              zii=temp_bins[trajii][timeii][2];
                if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
                {
                  if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii]){cout<<"\n"<<xii<<"\t"<<yii<<"\t"<<zii<<"\t"<<timeii<<"\t"<<trajii<<"\t"<<n_times<<"\t";cout.flush();}
                  if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii]){cout<<tempcount[xii][yii][zii][timeii]<<"\t"<<include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]<<"\t"<<trajcount[xii][yii][zii][timeii]<<"\t";cout.flush();}

                  include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]=trajii;
                  tempcount[xii][yii][zii][timeii]++;
                }
          }

    }
  }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	  delete [] tempcount[xii][yii][zii];
	}
	delete [] tempcount[xii][yii];
     }
     delete [] tempcount[xii];
   }
   delete [] tempcount;


  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {
      delete [] temp_bins[trajii][timeii];
    }
    delete [] temp_bins[trajii];
  }
  delete [] temp_bins;

}

void Trajectory_List_Bins::assign_bins_distance_plane(Trajectory_List * binned_list,string plane,float position, string direction, float bin_thickness)
{
  /**Rectangularly bins all trajs at all times within given dimensions (each bin has it's own Boolean_List)
     * @author Mark Mackura and David Simmons
     * @date 3/1/2013
     */
  int xii,yii,zii;
  int **** tempcount;
  //float xcoord,ycoord,zcoord;
  int *** temp_bins;

  temp_bins = new int ** [n_trajs];
  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    temp_bins[trajii]=new int * [n_times];
    for(int timeii=0; timeii<n_times; timeii++)
    {
      temp_bins[trajii][timeii]=new int [3];
    }
  }
  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
	{
	  for(int timeii=0; timeii<n_times; timeii++)
	  {
	    trajcount[xii][yii][zii][timeii]=0;
	  }
	}
     }
   }

  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {
        //cout << timeii<<"\t"<<trajii<<endl;cout.flush();
        //  cout << binned_list->is_included(timeii,trajii)<<"\n";cout.flush();
          Coordinate box_size(system->boundaries(timeii)[1](0)-system->boundaries(timeii)[0](0),system->boundaries(timeii)[1](1)-system->boundaries(timeii)[0](1),system->boundaries(timeii)[1](2)-system->boundaries(timeii)[0](2));

          if (binned_list->is_included(timeii,trajii))
          {

               float dist;
               float length;


                    if (plane == "x")
                    {
                        if (position > box_size.show_x())
                        {
                            cout << "plane chosen does not lie within the box dimmensions." << endl;cout.flush();
                            exit(1);
                        }
                            dist = system->show_trajectory(trajii)->show_coordinate(timeii).show_x() - position;

			    
                        if (direction == "average")
                        {
                            length = abs((box_size.show_x() * int((dist/(box_size.show_x()*.5))))-dist);
                        }
                        else if (direction == "above")
                        {
                            if (dist<0)
                            {
                                length =  box_size.show_x()+dist;
                            }
                            else
                            {
                                length = dist;
                            }
                        }
                        else if (direction == "below")
                        {

                            if (dist <0)
                            {
                                length = abs(dist);
                            }
                            else
                            {
                                length =  box_size.show_x()- dist;
                            }
                        }
                        else
                        {
                           cout << "could not recognize direction " << direction << ". valid choices are above, below or average.";cout.flush();
                        }
                        xii = int(length/bin_thickness);
                        yii =0;
                        zii = 0;
                    }
                    else if (plane == "y")
                    {
                        if (position > box_size.show_y())
                        {
                            cout << "plane chosen does not lie within the box dimmensions." << endl;cout.flush();
                            exit(1);
                        }
                        dist = system->show_trajectory(trajii)->show_coordinate(timeii).show_y() - position;

                        if (direction == "average")
                        {
                            length = abs((box_size.show_y() * int((dist/(box_size.show_y()*.5))))-dist);
                        }
                        else if (direction == "above")
                        {
                            if (dist<0)
                            {
                                length = box_size.show_y()+dist;
                            }
                            else
                            {
                                length = dist;
                            }
                        }
                        else if (direction == "below")
                        {

                            if (dist <0)
                            {
                                length = abs(dist);
                            }
                            else
                            {
                                length = box_size.show_y()-dist;
                            }
                        }
                        else
                        {
                           cout << "could not recognize direction " << direction << ". valid choices are above, below or average.";cout.flush();
                        }
                        xii = 0;
                        yii = int(length/bin_thickness);
                        zii = 0;
                    }
                    else if (plane == "z")
                    {
                        if (position > box_size.show_z())
                        {
                            cout << "plane chosen does not lie within the box dimmensions." << endl;cout.flush();
                            exit(1);
                        }
                        dist = system->show_trajectory(trajii)->show_coordinate(timeii).show_z() - position;

                        if (direction == "average")
                        {
                            length = abs((box_size.show_z() * int((dist/(box_size.show_z()*.5))))-dist);
                        }
                        else if (direction == "above")
                        {
                            if (dist<0)
                            {
                                length = box_size.show_z()+dist;
                            }
                            else
                            {
                                length = dist;
                            }
                        }
                        else if (direction == "below")
                        {

                            if (dist <0)
                            {
                                length = abs(dist);
                            }
                            else
                            {
                                length = box_size.show_z()-dist;
                            }
                        }
                        else
                        {
                           cout << "could not recognize direction " << direction << ". valid choices are above, below or average.";cout.flush();
                        }
                        xii = 0;
                        yii =0;
                        zii = int(length/bin_thickness);
                    }
                   else
                   {
                       cout << "could not recognize direction " << direction << ". valid choices are above, below or average.";cout.flush();
                   }



                        temp_bins[trajii][timeii][0]=xii;
                        temp_bins[trajii][timeii][1]=yii;
                        temp_bins[trajii][timeii][2]=zii;


                    if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
                    {
//                      cout <<xii<<"\t"; cout.flush();
//                        cout <<yii<<"\t"; cout.flush();
//                        cout <<zii<<"\t"; cout.flush();
//                        cout << timeii<<"\t"; cout.flush();
//                        cout <<trajcount[xii][yii][zii][timeii]<<endl; cout.flush();
                        trajcount[xii][yii][zii][timeii]++;
                    }
            }
    }
  }

   tempcount = new int***[n_xbins];
   for(int xii=0; xii<n_xbins; xii++)
   {
     tempcount[xii] = new int**[n_ybins];
     for(int yii=0; yii<n_ybins; yii++)
     {
        tempcount[xii][yii] = new int*[n_zbins];
        for(int zii=0; zii<n_zbins; zii++)
        {
            tempcount[xii][yii][zii] = new int[n_times];
        }
     }
   }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
        for(int zii=0; zii<n_zbins; zii++)
        {
          for(int timeii=0; timeii<n_times; timeii++)
          {
            include[xii][yii][zii][timeii]=new int[trajcount[xii][yii][zii][timeii]];
            tempcount[xii][yii][zii][timeii]=0;
          }
        }
     }
   }

  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {

          if (binned_list->is_included(timeii,trajii))
          {
              xii=temp_bins[trajii][timeii][0];
              yii=temp_bins[trajii][timeii][1];
              zii=temp_bins[trajii][timeii][2];
            if (xii<=n_xbins-1 && yii<=n_ybins-1 && zii<=n_zbins-1)
            {
              if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii]){cout<<"\n"<<xii<<"\t"<<yii<<"\t"<<zii<<"\t"<<timeii<<"\t"<<trajii<<"\t"<<n_times<<"\t";cout.flush();}
              if(tempcount[xii][yii][zii][timeii]==trajcount[xii][yii][zii][timeii]){cout<<tempcount[xii][yii][zii][timeii]<<"\t"<<include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]<<"\t"<<trajcount[xii][yii][zii][timeii]<<"\t";cout.flush();}

              include[xii][yii][zii][timeii][tempcount[xii][yii][zii][timeii]]=trajii;
              tempcount[xii][yii][zii][timeii]++;
            }
          }
    }
  }

  for(int xii=0; xii<n_xbins; xii++)
   {
     for(int yii=0; yii<n_ybins; yii++)
     {
 	for(int zii=0; zii<n_zbins; zii++)
 	{
	  delete [] tempcount[xii][yii][zii];
	}
	delete [] tempcount[xii][yii];
     }
     delete [] tempcount[xii];
   }
   delete [] tempcount;


  for(int trajii=0; trajii<n_trajs; trajii++)
  {
    for(int timeii=0; timeii<n_times; timeii++)
    {
      delete [] temp_bins[trajii][timeii];
    }
    delete [] temp_bins[trajii];
  }
  delete [] temp_bins;

}


void Trajectory_List_Bins::write_bins_xyz(string file)
{
  /** Writes an xyz traj file for each bin (all times)
  * @param file filename template to have bin indicies prepended to
  * @author Mark Mackura
  * @date 5/22/2012
  **/

  string bin_filename;
  int total_ii;
  total_ii = n_xbins*n_ybins*n_zbins;
  int done_ii=1;
  stringstream ssx,ssy,ssz;
  Trajectory_List bin;
  for(int xii = 0; xii<n_xbins; xii++)
  {
      for(int yii = 0; yii < n_ybins; yii++)
      {
        for(int zii = 0; zii < n_zbins; zii++)
        {
	  print_progress(done_ii++,total_ii);
	  ssx << (xii+1);
	  ssx.seekp(0);
	  ssy << (yii+1);
	  ssy.seekp(0);
	  ssz << (zii+1);
	  ssz.seekp(0);
	  bin_filename = file+"."+ssx.str()+"_"+ssy.str()+"_"+ssz.str()+".xyz";

	  bin.write_xyz(bin_filename);
        }
     }
  }
}
void Trajectory_List_Bins::write_single_bin_xyz(string file,int xii, int yii, int zii)
{
  /** Writes an xyz traj file for one bin (all times)
  * @param file filename template to have bin indicies prepended to
  * @param xii-yii-zii bin indicies
  * @author Mark Mackura
  * @date 5/22/2012
  **/

    Boolean_List * bins;

    bins = new Boolean_List[n_times];

    bin_boolean_calculation(xii,yii,zii,bins);

  cout<<"This method is broken."<<endl;
  exit(1);
  string bin_filename;
Trajectory_List bin;
  char buffer[MAXFILENAMECHAR];
  bin_filename = sprintf(buffer,"xbinii_%d_ybinii_%d_zbinii_%d_%s",xii,yii,zii,file.c_str());
  bin = Trajectory_List(system,n_times,n_trajs,bins,time_conversion);
  bin.write_xyz(bin_filename);
}
