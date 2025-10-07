<h1>region_multibody_list</h1>

Creates a new time-dependent multibody list based on an existing multibody list, by selecting at each time only those multibodies in the existing list that have centroid or center of mass (as originally selected) within a specified rectangular prism within the box. This check is done separately at each time in the trajectory such that the resulting multibody list will differ by time.

_region\_multibody\_list \<name of new multibody\_list to create\> \<name of existing multibody\_list to downselect geometrically\> \<xlo\> \<ylo\> \<zlo\> \<xhi\> \<yhi\> \<zhi\>_

\<xlo\> \<ylo\> \<zlo\> \<xhi\> \<yhi\> \<zhi\> set the lower and upper bounds of the domain to be selected in the x,y, and z directions. Note that this command does not actually create any new multibodies - it merely creates a new list downselected from an existing list of multibodies.
