Creates a list of all trajectories based on a binning scheme (specific subsets are called in analysis tools by taking the intersection of this list with a trajectory list defined with the create_list command).

_create\_bin\_list \<listname\>_

_\<type\> \<args\>_

Formatting of second line is:

|type	|	args | effect |
|----------|----------|----------|
|all	|		\<n\_xbins\> \<n\_ybins\> \<n\_zbins\>| All trajectories are binned into a (n_xbins,n_ybins,n_zbins) structure where the bins are evenly distributed across the dimensions of the entire system |
| region	|		\<n_xbins\> \<n_ybins\> \<n_zbins\> \<xlo\> \<xhi\> \<ylo\> \<yhi\> \<zlo\> \<zhi\> |All trajectories are binned similarly to <all>, but the bins are evenly split across the specified region.|
| distance	| trajectory \<list\_to\_bin\> \<list\_to\_take\_distance\_from\> \<bin\_thickness\> \<n\_bins\> | All trajectories are binned as a function of distance from the nearest particle in a second list. |
| distance |	point \<list\_to\_bin\> \<x\_coordinate\> \<y\_coordinate\> \<z\_coordinate\> \<bin\_thickness\> \<n\_bins\> | All trajectories are binned as a function of radial distance from the specified point.  |
| distance |	plane \<direction\> \<list\_to\_bin\> \<axes perpendicular to plane\> \<position of plane\> \<bin\_thickness\> \<n\_bins\> | All trajectories are binned as a function of perpendicular distance from the specified plane. \<direction\> can be “above” or “below”, which will only include particles above or below the plane, or “average” which will include particles on both sides of the plane. |
