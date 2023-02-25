//Add standard methods up here. Contructors, equality operator, etc.

void ParticlesBetween::timekernel2(int timeii)
{
   trajectory_list->listloop(this,0, timeii, 0);  
}


void Partices_Between::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Particles_Between::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii) //need to check back on which of these time indices is the right one.
{
   Trajectory* traj3;   
   bool isbetween;    //boolean to store whether the particle is in the 'between' set at this time
   
   int n_trajs=trajectory_list2->show_n_trajectories(thisii);  //record how many particles are in the second trajectory list at this time
   
   for(trajii=0;trajii<n_trajs;trajii++)  //here is the second loop over particles
   {
      traj3 = (*trajectory_list2)(thisii,trajii)  //this is the syntax to pull the second particle in set 2
      //we now have pointers to all 3 particles here.
      //traj1 is the central particle you are trying to sort
      //traj2 is the first particle in the group you are looking between
      //traj3 is the second particle in the group you are looking between
         
      //do some math here to determine whether particle is in between
      //if it is, set isbetween to true.
      //if it is not, set isbetween to false
      
         if(isbetween)
         {
            addtrajectory(thisii,traj1);        //this line will add the trajectory to the trajectory list
         }
   }
      
   float distance;
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    bin(thisii,distance);
  }
   
}
