
/*Multibody_List_Bins class - stores a list of multibodies according to bins
Amorphous Molecular Dynamics Tookit (AMDAT)
Written by David Simmons based off of early code by Mark Mackura*/
#ifndef MULTIBODY_LIST_BINS
#define MULTIBODY_LIST_BINS

#include "multibody_list.h"
#include "system.h"
#include "boolean_list.h"
#include "trajectory_list.h"
#include "coordinate.h"
#include "progress.h"

#define MAXFILENAMECHAR	255 //maximum allowable filename size by linux 2.6.32-5

namespace std{

class Multibody_List_Bins
{
    System * system;
    int n_times;                                        //number of system times;
    int n_trajs;                                        //number of total binned trajectories
    int n_xbins, n_ybins, n_zbins;                      //number of bins in each direction
    float *xlo,*xhi,*ylo,*yhi,*zlo,*zhi;                //box dimensions
    Multibody_List***include;
    //int ***** include;
    //int **** trajcount;

    int * time_conversion;			        //time_conversion array pointer (note this is the same for all bins)

    void bin_boolean_calculation(int,int,int,Boolean_List*);
    //void set_trajectory_lists();

    void assign_bins();
    void assign_bins_distance_clusters(Trajectory_List*,Trajectory_List*,float);
    void assign_bins_distance_point(Trajectory_List*,Coordinate,float);
    void assign_bins_distance_plane(Trajectory_List*,string,float,string,float);

    public:
      Trajectory_List_Bins(); //default constructor
      ~Trajectory_List_Bins(); //destructor
      Trajectory_List_Bins(const Trajectory_List_Bins &); //copy constructor
      Trajectory_List_Bins(System * sys,int xbins = 1,int ybins = 1,int zbins = 1);                                                                  //bins entire system rectangularly
      Trajectory_List_Bins(System * sys,int,int,int,float,float,float,float,float,float);                              //bins specified region rectangularly using absolute dimensions
      Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, Trajectory_List* clustered_list);     //bins a list as a distance from particles in another list
      Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, string pln, float posit, string dir);     //bins a list as a distance from particles in another list
      Trajectory_List_Bins(System * sys, float bin_thickness, int n_bins, Trajectory_List* binning_list, Coordinate pnt);     //bins a list as a distance from particles in another list
      /*Don't use*/Trajectory_List_Bins(System * sys,int,int,int,Coordinate,Coordinate);                                        //bins region between two coordinates rectangularly
      Trajectory_List_Bins operator=(const Trajectory_List_Bins &);	//assignment

      /* Methods to return Trajectory lists for specific bins */
      Trajectory_List persistant_trajectories(int xii,int yii,int zii,int starttime,int endtime);       //returns trajectory_list of particles in specified bin for time duration

      /* Methods to give information about bin structure, sizes */
      int show_n_xbins();
      int show_n_ybins();
      int show_n_zbins();
      float show_lx(int timeii=0);
      float show_ly(int timeii=0);
      float show_lz(int timeii=0);

      /* Methods to give information about Trajectories in bins */

      Trajectory_List operator() (int,int,int);     	//returns trajectory_List of bin xii,yii,zii
      void write_bins_xyz(string);		//writes xyz file for each bin, prepending bin index to filename
      void write_single_bin_xyz(string,int,int,int);	//writes xyz file for single bin
};
}
#endif
