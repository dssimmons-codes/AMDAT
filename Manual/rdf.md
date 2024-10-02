<h1>rdf</h1>

Computes a radial distribution function. The radial distribution may in general be asymmetric, computed between two distinct sets of atoms $a$ and $b$ ($a$ and $b$ may be the same). When asymmetric, $a$ is the central particle in the g(r) calculation and $b$ is the set of particles for which densities around the central particles are computed.

This is computed as follows.

$$  g(r) = \frac{1}{\sum_{i=1}^{T} N^a(t_i)} \sum_{i=1}^{T} \sum_{k=1}^{N^a(t_i)}  \frac{n_{k}^{b}(r \pm \Delta r,t_i)}{V(r \pm \Delta r) \rho^b(t_i)} $$

Here T is the number of times in the trajectory, $N^a(t_i)$ is the number of particles in group $a$ at time $t_i$.


_rdf \<output filename\> \<symmetry keyword\> \<number of distance bins\> \<timescheme\> \<max length scale to compute to\>_

_\<symmetry keyword\>_ is either "symmetric" or "asymmetric". If "symmetric", then only a single target trajectory list line is provided. The RDF is then computed between all particles of this list in a symmetric manner. If "asymmetric", then two target trajectory list lines are provided. The RDF is then constructed by binning members of the second list according to distance from the first list, leading to an asymmetric radial distribution function.

_\<number of distance bins\>_ specifies at how many points (distances) to compute the radial distibution funcion

_\<timescheme\>_ determines what times to loop over. If timescheme = -1, loop over all times. If timescheme is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block.

_\<max length scale to compute to\>_ sets the maximum distance to which the rdf is computed. If this is set to a value of 0, the calculation defaults to something like half the smallest box dimension.
