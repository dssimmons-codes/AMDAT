<h1>composition_vs_time</h1>

Calculates the composition and number density.

_composition\_vs\_time \<output file\>_

_\<target\>_

This method is most commonly employed with spatially binned analysis to compute gradients in particle density, number, or composition with spatial position in a simulation box.

A header section reports a time-average number density and particle fraction of each species (per bin in the case of binned analysis). This is followed by a section that reports, for each time in the trajectory, a number density, total number of particles, and then particle fractions for each species.

For the purposes of computing a density, this method currently assumes the system (or a bin) to have a time-invariant volume. The volume is based on the first timestep. For trajectory/_list analysis, the volume employed is the entire box. For binned analysis using any rectilinear bin shape, the volume of the bin at the first timestep is employed. For trajectory/_bin/_lists constructed based on distance to a second trajectory/_list, the bin volume is not well defined and is treated as -1 so that density outputs are negative and clearly indicated as not meaningful to the user. For these situations, bin counts and compositions can be employed instead.
