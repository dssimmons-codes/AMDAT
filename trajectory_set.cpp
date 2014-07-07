/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Trajectory_set class - Class to store sets of user-defined trajectories
/*Written by David S. Simmons*/


Trajectory_Set::Trajectory_Set()
{
    n_trajectories=0;
    trajectories=new Trajectory [n_trajectories];
}

Trajectory_Set::Trajectory_Set(const Trajectory_Set & copy)
{
    int trajii;

    n_trajectories = copy.n_trajectories;
    trajectories=new Trajeectory [n_trajectories];

    for(trajii=0;trajii<n_trajectories;trajii++)
    {
        trajectories[trajii]=copy.trajectories[trajii];
    }
}

Trajectory_Set::~Trajectory_Set()
{
    delete [] trajectories;
}


Trajectory_Set::operator=(const Trajectory_Set & copy)
{
    int trajii;

    if(this!=&copy)
    {
        delete [] trajectories;

        n_trajectories = copy.n_trajectories;
        trajectories=new Trajeectory [n_trajectories];

        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=copy.trajectories[trajii];
        }
    }
}



Trajectory_Set (const Multibody_Set & multibodies, bool centertype)
{
    int trajii;
    n_trajectories = multibodies.show_n_multibodies();

    trajectories = new Trajectory [n_trajectories];

    if(centertype) //compute and save centroid trajectory
    {
        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=multibodies.show_multibody[trajii]->centroid_trajectory;
        }
    }
    else    //compute and save COM trajectory
    {
        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=multibodies.show_multibody[trajii]->center_of_mass_trajectory;
        }
    }
}
