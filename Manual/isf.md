<h1>isf</h1>
<h2>function</h2>

Calculates the full intermediate scattering function of the system via the Fourier transform of the density.

<h2>Syntax</h2>

```
structure_factor <output file> <method> <symmetry> <plane> <timescheme> <minimum wavenumber index> <maximum wavenumber index>
<target>
<target 2> (if symmetry = asymmetric only)
```

`<method>` is either `auto` or `manual`. The `manual` functionality is not currently documented; employ auto.

`<symmetry>` is either `symmetric` or `asymmetric`.
If `symmetric`, the analysis calculates the structure factor between the set of particles specified in target and itself.
If `asymmetric`, then a second target must be specified, and the partial structure factor describing correlations of `<target>` with only `<target 2>` is calculated.
`<first frame>` and `<last frame>` are the indices of the limits on time spacings to be calculated.
Options for ``<plane>`` are ``xyz``, ``xy``, ``xz``, ``yz``, ``x``, ``y``, and ``z``.
This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function.
``xyz`` computes the full radial three dimensional isf, ``xy``, ``yz``, and ``xz`` calculate two-dimensional in-plane radial isf's, and ``x``, ``y``, and ``z`` compute one-dimensional isf's.
``<max_length_scale>`` determines the longest distance which will be decomposed into inverse space.
If a distance of 0 is given, the full box size is used.
Any deviation from '0' will in general produce incorrect results for the structure factor, especially at low k.

`<timescheme>` determines what times to loop over.
If `timescheme` is -1, loop over all times.
If `timescheme` is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block.
In most cases this should be set to either 0 or -1, with -1 giving improved statistical strength at the cost of much longer compute times.
The level of improvement in statistical strength will depend on the timescale for structural decorrelation in comparison to the length of a trajectory block.
