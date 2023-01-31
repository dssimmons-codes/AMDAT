#ifndef LINKED_TRAJECTORIES
#define LINKED_TRAJECTORIES


namespace std{


class Linked_Trajectories
{
        Trajectory* base_trajectory;
        vector <Trajectory*> linked_trajectories;

    public:
        Linked_Trajectories();
        Linked_Trajectories(Trajectory*);
        Linked_Trajectories(Trajectory*,vector<Trajectory*>);
        Linked_Trajectories(Trajectory*,Trajectory**,int);

        Linked_Trajectories operator &&(const Linked_Trajectories& compare) const;    //return new Linked_Trajectories object containing only those trajectories in both lists
        Linked_Trajectories operator ||(const Linked_Trajectories& compare) const;    //return new Linked_Trajectories list containing all trajectories in both
        Linked_Trajectories operator-(const Linked_Trajectories& compare) const;    //return new Linked_Trajectories object equal to this one but with any trajectories removed that are also contained in compare
        bool operator <(int);
        bool operator >(int);
        bool operator <=(int);
        bool operator >=(int);
        bool operator <(const Linked_Trajectories &);
        bool operator >(const Linked_Trajectories &);
        bool operator <=(const Linked_Trajectories &);
        bool operator >=(const Linked_Trajectories &);

        void operator +=(Trajectory*);      //add a trajectory

        operator float() const{return float(size());};
        operator int() const{return size();};

        int size()const{return trajectories.size();};     //return number of trajectories in object



};


}


#endif // LINKED_TRAJECTORIES
