<h1>radial_count</h1>

Computes a non-normalized radial distribution function - essentially the mean number of particles in shells a distance $r$ to $r+\Delta r$ from a central particle. 

Syntax is as follows.

_rdf \<output filename\> \<symmetry keyword\> \<number of distance bins\> \<timescheme\> \<max length scale to compute to\>_

_\<symmetry keyword\>_ is either "symmetric" or "asymmetric". If "symmetric", then only a single target trajectory list line is provided. The RDF is then computed between all particles of this list in a symmetric manner. If "asymmetric", then two target trajectory list lines are provided. The RDF is then constructed by binning members of the second list according to distance from the first list, leading to an asymmetric radial distribution function.

_\<number of distance bins\>_ specifies at how many points (distances) to compute the radial distibution funcion

_\<timescheme\>_ determines what times to loop over. If timescheme = -1, loop over all times. If timescheme is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block.

_\<max length scale to compute to\>_ sets the maximum distance to which the rdf is computed. If this is set to a value of 0, the calculation defaults to something like half the smallest box dimension.

