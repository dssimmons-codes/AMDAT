#include <math>


void Radial_Distribution_Function::analyze(Trajectory_List * t_list, Trajectory_List * t_list2)
{
  int timeii;
  trajectory_list=t_list;
  trajectory_list2=tlist2;

  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms_1[timeii]=trajectory_list->show_n_trajectories(timeii);
    n_atoms_2[timeii]=trajectory_list2->show_n_trajectories(timeii);
    trajectory_list->listloop(this,0, timeii, 0);
  }
  
  postprocess_list();
}


void Radial_Distribution_Function::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  trajectory_list2->listloop(this,0, thisii, 0);
}


void Radial_Distribution_Function::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(thisiiii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    bin(timeii,distance);
  }
}


void Radial_Distribution_Function::bin(int timestep, float distance)
 {
  int binindex;
  binindex = int((distance)/bin_size);
    
  if(binindex>=0)
  {
    if(binindex<n_bins)
    {
      (time_rdf[timestep][binindex])++;
    }
  }
 }
 
 
 void Radial_Distribution_Function::postprocess_list()
 {
   int timeii,binii;
   float boxvolume;
   Coordinate boxsize;
   float rhoj;
   float rshell,shellvolume;
   int n_i_total=0;
   
   for(timeii=0;timeii<n_times;timeii++)
   {
     boxsize = system->size(timeii);
     boxvolume=boxsize.show_x()*boxsize.show_y()*boxsize.show_z();
     rhoj=float(n_atoms_j[timeii])/boxvolume;
     n_i_total+=n_atoms_i[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       rshell = float(binii)*bin_size;						//determine inner radius of bin
       shellvolume = (4.0/3.0)*PI*(pow(rshell+bin_size,3.0)-pow(rshell,3.0));		//calculate volume of bin
       time_rdf[timeii][binii]/=(float(n_atoms_i[timeii])*rhoj*shellvolume);
       mean_rdf[binii]+=time_rdf[timeii]*float(n_atoms_i[timeii]);
     }
   }
   
   for(binii=0;binii<n_bins;binii++)
   {
     mean_rdf[binii]/=float(n_i_total);
   }
 }
 
 #ifdef NEVER
 
 Radial_Distribution_Function operator+ (const Radial_Distribution_Function & increment)
 {
   
   Radial_Distribution_Function() temp;
   
   
   //compute time_rdf
   for(timeii=0;timeii<n_times;timeii++)
   {
     temp.n_atoms_i[timeii]=n_atoms_i+increment[timeii].i_atoms_j[timeii];
     temp.n_atoms_j[timeii]=n_atoms_j+increment[timeii].n_atoms_j[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       time_rdf
     }
   }
 
   
   //compute and normalize mean_rdf
   int n_i_total = 0;
   for(timeii=0;timeii<n_times;timeii++)
   {
     n_i_total+=n_atoms_i[timeii];
     for(binii=0;binii<n_bins;binii++)
     {
       temp.mean_rdf[binii]+=time_rdf[timeii]*float(n_atoms_i[timeii]);
     }
   }
   
   for(binii=0;binii<n_bins;binii++)
   {
     temp/ mean_rdf[binii]/=float(n_i_total);
   }
 }
 
 #endif
 
 void Radial_Distribution_Function::write(string filename)
 {
   
 }