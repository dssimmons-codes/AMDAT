<h1>delete\_trajectories</h1>

Deletes an existing trajectory_list, freeing up the associated memory.

_delete_trajectory_list \<name of existing trajectory_list to delete\>_

In most cases, trajectory\_lists require only a modest amount of memory to store, such that this is not strictly necessary. However, if an analysis tool defines a trajectory list that changes at every timestep in the trajectory, it required memory storage may be considerably higher, and deletion may be worthwhile if you are memory limited. 
