<h1>structure_factor</h1>
<h2>function</h2>

Calculates the structure factor of the system via the Fourier transform of the density, as per the following equations.

$$ S\left( {\vec{k}} \right) = \left\langle \frac{\tilde{\rho}^{I} \left( \vec{k} \right) \tilde{\rho}^{II} \left( -\vec{k} \right)}{N}  \right\rangle $$

$$ {\tilde{\rho}}^{I} \left( \vec{k} \right) = \sum_{i=1}^{{{N}^{I}}} \exp \left( -i \vec{k}_{j} \cdot \vec{r}\_{i} \right)$$

(For scientific background and discussion, see Chapter 4 of Hansen and McDonald’s “Theory of Simple Liquids). Here the superscripts I and II denote the set of particles under consideration. When I and II refer to the same set of particles, this is the symmetric structure factor. When they are different sets of particles, this is an asymmetric structure factor reporting on the cross-correlations between distinct sets of particles.

Mechanically, AMDAT first computes the Fourier transform of the density for each required wavevector as follows.

$$ \tilde{\rho} \left( \vec{k}\_j, t \right) = \sum_{i=1}^{N} \cos \left( \vec{k}\_{j} \cdot \vec{r}\_{i} \left( t \right) \right) + i \sum_{i=1}^{N} \sin \left( \vec{k}\_j \cdot \vec{r}\_i \left( t \right) \right) $$
<!-- ORIGINAL EQUATION: $$ \[\tilde{\rho} \left( {{{\vec{k}}}_{j}},t \right)=\sum\limits_{i=1}^{N}{\cos \left( {{{\vec{k}}}_{j}}\cdot {{{\vec{r}}}_{i}}\left( t \right) \right)}+i\sum\limits_{i=1}^{N}{\sin \left( {{{\vec{k}}}_{j}}\cdot {{{\vec{r}}}_{i}}\left( t \right) \right)}\] $$ -->

If the computation being performed is a symmetric structure factor, this calculation is performed once. If it is an asymmetric structure factor, it is performed distinctly for each set of particles. These calculations are separately performed at each timestep being analyzed.

AMDAT ultimately outputs the structure factor as a function of scalar wavenumber rather an as a function of wavevector. To do so, it averages over many wavevectors with approximately equivalent wavenumber (more on this below). It also averages over many times. AMDAT thus computes the structure factor at a wavenumber k as follows.

$$ S\left( k \right) = \frac{1}{SH} \sum \limits_{l=1}^{S}{\sum\limits\_{j=1}^{H}{\left[ \tilde{\rho}\_{real}^{I}\left( {{{\vec{k}}}\_{j}},{{t}\_{l}} \right)\tilde{\rho}\_{real}^{II}\left( {{{\vec{k}}}\_{j}},{{t}\_{l}} \right)+\tilde{\rho}\_{imag}^{I}\left( {{{\vec{k}}}\_{j}},{{t}\_{l}} \right)\tilde{\rho}\_{imag}^{II}\left( {{{\vec{k}}}\_{j}},{{t}\_{l}} \right) \right]}} $$
<!-- $$ S\left( k \right)=\frac{1}{SH}\sum\limits_{l=1}^{S}{\sum\limits_{j=1}^{H}{\left[ \tilde{\rho} _{real}^{I}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)\tilde{\rho} _{real}^{II}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)+\tilde{\rho} _{imag}^{I}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right)\tilde{\rho} _{imag}^{II}\left( {{{\vec{k}}}_{j}},{{t}_{l}} \right) \right]}} $$
 -->

where $S$ is the number of times over which $S$ is averaged and $H$ is the number of wavevectors corresponding to the wavenumber $k$. 

```
structure_factor <output file> <symmetry> <geometry> <max_length_scale> <(optional) timescheme>

<target>

<target 2> (optional, if symmetry = asymmetric only)
```

Options for `<geometry>` are `xyz`, `xy`, `xz`, `yz`, `x`, `y`, and `z`. This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function. `xyz` computes the full radial three dimensional `isf`, `xy`, `yz`, and `xz` calculate two-dimensional in-plane radial `isf`'s, and `x`, `y`, and `z` compute one-dimensional `isf`'s. `<symmetry>` is either `symmetric` or `asymmetric`. If `symmetric`, the analysis calculates the structure factor between the set of particles specified in target and itself. If `asymmetric`, then a second target must be specified, and the partial structure factor describing correlations of `<target 1>` with only `<target 2>` is calculated. `<first frame>` and `<last frame>` are the indices of the limits on time spacings to be calculated.

`<max_length_scale>` determines the longest distance which will be decomposed into inverse space, if a distance of 0 is given, the full box size is used. 

`<(optional) timescheme>` determines what times to loop over. If `timescheme` is -1, loop over all times. If `timescheme` is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block. 
In most cases this should be set to either 0 or -1, with -1 giving somewhat improved statistical strength at the cost of much longer compute times.
