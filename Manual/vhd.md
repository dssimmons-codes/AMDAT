<h1>vhd</h1>

<h2>Function</h2>

Calculates distinct part of the Van Hove function. The distinct van Hove is by construction asymmetric, computed between two distinct sets of atoms $a$ and $b$ ($a$ and $b$ may be the same). When asymmetric, $a$ is the central particle in the $G(r,\Delta t)$ calculation and $b$ is the set of particles for which densities around the central particles are computed.

$$  G(r,\Delta t) = \frac{1}{\sum_{i=1}^{T} N^a(s_i)} \sum_{i=1}^{T} \sum_{k=1}^{N^a(t_i)}  \frac{n_{k}^{b}(r \pm \frac{\Delta r}{2},s_i,s_i + \Delta t)}{V_{shell}(r \pm \frac{\Delta r}{2})} $$

Here T is the number of starting times to be employed in the dynamical analysis (in a blocked trajectory this is equal to the number of blocks), $N^a(t_i)$ is the number of particles in group $a$ at time $t_i$, $n_{k}^{b}(r \pm \frac{\Delta r}{2},t_i)$ is the number of particles of type $b$, at time $s_i+\Delta t$ within a shell extending from $r-\frac{\Delta r}{2}$ to $r+\frac{\Delta r}{2}$ around the location of a particle $k$ of type $a$ at an earlier time $s_i$, and $V_{shell}(r \pm \frac{\Delta r}{2})$ is the volume of this shell.

<h2>Syntax</h2>

_vhd \<output file\> \<min cell size\> \<maximum range binned\> \<# of bins\>_

_\<target\>_

This algorithm uses a spatial decomposition of the box into cells, with cells guaranteed to be of size _\<min cell size\>_ or larger (the algorithm will select the minimum cell size that is an integer divisor of the box size and is larger than _\<min cell size\>_.
