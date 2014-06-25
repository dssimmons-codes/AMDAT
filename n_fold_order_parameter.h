/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class to calculate n_fold order parameter*/
/*Written by Daniel A. Hunsicker*/


#ifndef N_FOLD_ORDER_PARAMETER
#define N_FOLD_ORDER_PARAMETER

#include "analysis.h"
#include "value_list.h"
#include <string>
#include <sstream>

namespace std{


class N_Fold_Order_Parameter : public Value_List<float>, public Analysis
{
    int total_atoms;
    int start_time;
    float order;
    int end_time;
    int n_bins;
    int n_included;
    float sig_cut;
    float threshold;
    float * distribution;
    string pdb_stem;
    int n_atomtypes;
    int current_time;
    float max_param;
    bool want_map;
    float ** sigmatrix;
    float * timetable;

    string * species_name;
    Boolean_List * n_fold_thresh;
    float time_average_param;
    float max_sigma;
//    Value_List * orientation
    int neighborcount;
    float param_jth;
    float param_total;
    float * param_total_current;
    int atomcount_total;      //temp atomcount variable
    int atomcount;

    typedef float (Coordinate::*length)()const;
    length distancefun;

    public:
    N_Fold_Order_Parameter(System * sys, float ord, string sig_file, string orientation = "xy", float cut = 1.25, string file_stem = "map_hop", int start=0, int end=-1);

    Analysis_Type what_are_you(){Analysis_Type type = n_fold_order_parameter; return type;};		//virtual method to report the type of analysis
    
    void analyze(Trajectory_List * t_list);
    void listkernel(Trajectory *);
    void write(string);
    void find_ordered();
    void set_time_conv();


};
}
#endif
