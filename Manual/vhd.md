<h1>vhd</h1>

Calculates distinct part of the Van Hove function.

_vhd \<output file\> \<min cell size\> \<maximum range binned\> \<# of bins\>_

_\<target\>_

This algorithm uses a spatial decomposition of the box into cells, with cells guaranteed to be of size _\<min cell size\>_ or larger (the algorithm will select the minimum cell size that is an integer divisor of the box size and is larger than _\<min cell size\>_.
