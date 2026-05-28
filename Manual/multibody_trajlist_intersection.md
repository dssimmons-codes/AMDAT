<h1>multibody_trajlist_intersection</h1>

Creates a new time-dependent multibody list based on an existing multibody list, by selecting at each time only those multibodies in the existing list that have at least \<threshold\> constituent beads within a specified `trajectory_list` at that time. This check is done separately at each time in the trajectory such that the resulting multibody list will differ by time.

multibody\_trajlist\_intersection \<name of new multibody\_list to create\> \<name of existing multibody\_list to downselect geometrically\> \<name of trajectory\_list to check against\> \<threshold\>

Note that this command does not actually create any new multibodies - it merely creates a new list downselected from an existing list of multibodies. 
