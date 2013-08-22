#ifndef CLUSTERED_LIST_H
#define CLUSTERED_LIST_H

#include <iostream>
#include "trajectory_list.h"
#include "analysis.h"
#include "system.h"
namespace std{

class Clustered_List: public Trajectory_List
{   protected:
    float ** sigmatrix;
    int *** neighbor;


    public:

        Clustered_List();
        ~Clustered_List();
        Clustered_List(const Trajectory_List & t_list):Trajectory_List(t_list){initialize_members();};
        void initialize_members();
        void allocate_sig_matrix(string);

        void construct_clust_list(Trajectory_List*,string,string,int,int);
        void generate_neighbor_lists(int, int,string);
        bool cluster(int,int,int);
        bool cluster_correction(int,int,int,Boolean_List);
        bool cluster_correction2(int,int,Boolean_List);

};

}
#endif // CLUSTERED_LIST_H
