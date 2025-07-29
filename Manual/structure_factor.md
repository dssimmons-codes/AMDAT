<h1>structure_factor</h1>
<h2>function</h2>

Calculates the structure factor of the system via the Fourier transform of the density, as per the following equations.

$$ S\left( {\vec{k}} \right) = \left\langle \frac{1}{N}\rho _{{\vec{k}}}^{I}\rho _{-\vec{k}}^{II} \right\rangle $$

$$ {{\rho }^{I}}\left( {\vec{k}} \right)=\sum\limits_{i=1}^{{{N}^{I}}}{\exp \left( -i{{{\vec{k}}}_{j}}\cdot {{{\vec{r}}}_{i}} \right)} $$

(For scientific background and discussion, see Chapter 4 of Hansen and McDonald’s “Theory of Simple Liquids). Here the superscripts I and II denote the set of particles under consideration. When I and II refer to the same set of particles, this is the symmetric structure factor. When they are different sets of particles, this is an asymmetric structure factor reporting on the cross-correlations between distinct sets of particles.

Mechanically, AMDAT first computes the Fourier transform of the density for each required wavevector as follows.

$$ \[\rho \left( {{{\vec{k}}}_{j}},t \right)=\sum\limits_{i=1}^{N}{\cos \left( {{{\vec{k}}}_{j}}\cdot {{{\vec{r}}}_{i}}\left( t \right) \right)}+i\sum\limits_{i=1}^{N}{\sin \left( {{{\vec{k}}}_{j}}\cdot {{{\vec{r}}}_{i}}\left( t \right) \right)}\] $$

If the computation being performed is a symmetric structure factor, this calculation is performed once. If it is an asymmetric structure factor, it is performed distinctly for each set of particles. These calculations are separately performed at each timestep being analyzed.

AMDAT ultimately outputs the structure factor as a function of scalar wavenumber rather an as a function of wavevector. To do so, it averages over many wavevectors with approximately equivalent wavenumber (more on this below). It also averages over many times. AMDAT thus computes the structure factor at a wavenumber k as follows.

$$ \[S\left( k \right)=\frac{1}{SH}\sum\limits_{l=1}^{S}{\sum\limits_{j=1}^{H}{\left[ \rho _{real}^{I}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)\rho _{real}^{II}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)+\rho _{imag}^{I}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)\rho _{imag}^{II}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right) \right]}}\] $$

where S is the number of times over with S is averaged and H is the number of wavevectors corresponding to the wavenumber k. 



_structure\_factor \<output file\> \<symmetry\> \<geometry\> \<max\_length\_scale\> \<(optional) timescheme\>_

_\<target\>_

_\<target 2> (optional, if symmetry = asymmetric only)_

Options for _\<geometry\>_ are _xyz, xy, xz, yz, x, y,_ and _z_. This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function. _xyz_ computes the full radial three dimensional isf, _xy, yz, and xz_ calculate two-dimensional in-plane radial isf's, and _x, y, and z_ compute one-dimensional isf's. _\<symmetry\>_ is either _symmetric_ or _asymmetric_. If _symmetric_, the analysis calculates the structure factor between the set of particles specified in target and itself. If _asymmetric_, then a second target must be specified, and the partial structure factor describing correlations of _\<target 1\>_ with only _\<target 2\>_ is calculated. _\<first frame\>_ and _\<last frame\>_ are the indices of the limits on time spacings to be calculated.

_\<max\_length\_scale\>_ determines the longest distance which will be decomposed into inverse space, if a distance of 0 is given, the full box size is used. 

_\<(optional) timescheme\>_ determines what times to loop over. If timescheme = -1, loop over all times. If timescheme is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block. 
In most cases this should be set to either 0 or -1, with -1 giving somewhat improved statistical strength at the cost of much longer compute times.
