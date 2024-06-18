<h1>structure_factor</h1>

Calculates structure factor.

_structure\_factor \<output file\> \<symmetry\> \<geometry\> \<max\_length\_scale\> \<(optional) fullblock\>_

_\<target\>_

_\<target 2(optional, if symmetry = asymmetric only)_

Options for _\<geometry\>_ are _xyz, xy, xz, yz, x, y,_ and _z_. This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function. _xyz_ computes the full radial three dimensional isf, _xy, yz, and xz_ calculate two-dimensional in-plane radial isf's, and _x, y, and z_ compute one-dimensional isf's. _\<symmetry\>_ is either _symmetric_ or _asymmetric_. If _symmetric_, the analysis calculates the structure factor between the set of particles specified in target and itself. If _asymmetric_, then a second target must be specified, and the partial structure factor describing correlations of _\<target 1\>_ with only _\<target 2\>_ is calculated. _\<first frame\>_ and _\<last frame\>_ are the indices of the limits on time spacings to be calculated.

_\<max\_length\_scale\>_ determines the longest distance which will be decomposed into inverse space, if a distance of 0 is given, the full box size is used_. \<(optional) fullblock\>_is 0 if structure factor should only be calculated at the start of each block, 1 if it should be calculated at all times
