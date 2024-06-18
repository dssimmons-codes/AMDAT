<h1>invert_list</h1>

Creates a new trajectory list by inverting an existing trajectory list (meaning it unincludes all included trajectories and includes all unincluded trajectories) and then intersecting it with a second list (meaning it includes only those trajectories in both lists). The net effect of this is to include only trajectories that are not in trajectory\_list 1 but are in trajectory\_list 2.

invert\_list \<initial trajectory\_list\> \<second trajectory\_list\> \<name of new inverted list to be created\>_
