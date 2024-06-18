<h1>thresholded_list</h1>

Creates a thresholded list of trajectories based on a specified analysis\_value\_list.

_thresholded\_list \<analysis\_value\_list name\> \<thresholded\_list name\>\<threshold type\>\<first threshold\>\<optional:second threshold\>_

An analysis\_value\_list is an object that specifies a value (for example, a float) for some or all of the particle trajectories in the system. For example, this could be used to store some measure of local order around each particle. Thresholded list could then be employed to construct a trajectory\_list containing the most ordered particles at each time.

Allowable values for _threshold type_ are

_greater_ – selects trajectories with values greater than the first threshold

_less_ – selects trajectories with values less than the first threshold

_between_ – selects trajectories with values between and including the first and second threshold
