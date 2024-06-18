<h1>isfs</h1>

Calculates self-part of the intermediate scattering function. 

_isfs \<output file\> \>\<first wavevector index\> \<last wavevector index\> \<geometry\> \<max\_length\_scale\> \<(optional) fullblock\>_

_\<target\>_

_\<first wavevector index\> and \<lastwavevector index\>_ give the limits wave vectors to be calculated. Options for _\<geometry\>_ are _xyz, xy, xz, yz, x, y,_ and _z_. This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function. _xyz_ computes the full radial three dimensional isf, _xy, yz, and xz_ calculate two-dimensional in-plane radial isf's, and _x, y, and z_ compute one-dimensional isf's. _\<max\_length\_scale\>_ determines the longest distance which will be decomposed into inverse space. If a distance of 0 is given, the full box size is used_. \<(optional) fullblock\>_ is an optional argument that can be either 0 or 1. The default is 0, in which case time spacings spanning multiple blocks use only the first time of each block. A setting of 1 specifies that it should use all block times for times spanning blocks. This may result in substantial computational time increases but offers the possibility of modestly increased data quality at very long times.
