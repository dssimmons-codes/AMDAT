<h1>Index of Commands</h1>

| Command | Effect | Input | Output |
|----------|----------|----------|----------|
| [clustered\_list](clustered_list_decay.md) | Generates a new trajectory list by finding clustered particles in the original list. | trajectory\_list | New trajectory\_list |
| [compare_gaussian](compare_gaussian.md) | Finds points where the self Van Hove crosses the Gaussian approximation to the self Van Hove. | trajectory\_list | Tab demarkated data file |
| [composition](composition.md) | Calculates the composition and number density. | trajectory\_list | Tab demarkated data file |
| create_distance_neighborlist | Constructs a neighborlist object based on a distance cutoff criterion. | trajectory\_list | neighbor_list |
| create_voronoi_neighborlist | Constructs a neighborlist object based on the voronoi tesselation | trajectory\_list | neighbor_list |
| displacement_list | Computes the displacement of each particle at some time separation and stores the results in a value_list | trajectory\_list | Tab demarkated data file and value_list |
| [displacement_map](displacement_map.md) | Generates maps of particle displacement magnitude and writes to pdb files with value beta specifying displacement. | trajectory\_list | PDB file |
| [isf\_list](isf_list.md) | Calculate full intermediate scattering function | trajectory\_list | Tab demarkated data file |
| [isfs](isfs.md) | Calculates the self-part of the intermediate scattering function. | trajectory\_list | Tab demarkated data file |
| [mean\_closest\_distance](mean_closest_distance.md) | Calculates the average distance between particles in two distinct trajectory\_lists. | trajectory\_list | Tab demarkated data file |
| [msd](msd.md) | Calculates the mean-square displacement. | trajectory\_list | Tab demarkated data file |
| [msd\_2d](msd_2d.md) | Calculates the two-dimensional mean-square displacement. | trajectory\_list | Tab demarkated data file |
| [n\_fold](n_fold.md) | Calculates mean and time dependence of n-fold order parameter and writes maps to files. | trajectory\_list | Tab demarkated data file and a PDB file |
| [ngp](ngp.md) | Calculates non-Gaussian parameter of the mean squared displacement. | trajectory\_list | Tab demarkated data file |
| [rdf](rdf.md) | Computes the radial distribution function. | trajectory\_list | Tab demarkated data file |
| [rnf](radial_count.md) | Computes the mean number of particles in shells around a set of central particles | Tab demarcated data file |
| [rgtensor\_stats](rgtensor_stats.md) | Calculates statistics describing the radius of gyration tensor characterizing particle trajectories as a function of time. | trajectory\_list | Tab demarkated data file |
| [stiffness\_dist](stiffness_dist.md) | Calculates distribution of inverse Debye-Waller factor values 1/u2. | trajectory\_list | Tab demarkated data file |
| [strings](strings.md) | Analyzes 'stringlike' collective motion. Planned for deprecation in favor of string_multibodies. | trajectory\_list | Tab demarkated data file. |
| string_multibodies | Builds multibodies comprised of 'stringlike' collective displacements | trajectory\_list | multibody_list |
| [structure\_factor](structure_factor.md) | Computes the structure factor. | trajectory\_list | Tab demarkated data file |
| [trajectory\_list\_decay](trajectory_list_decay.md) | Computes an autocorrelation function for inclusion in the trajectory\_list. | trajectory\_list | Tab demarkated data file |
| [u2dist](u2dist.md) | Calculates distribution of square displacements at a specified time. | trajectory\_list | Tab demarkated data file |
| [vac\_fourier](vac_fourier.md) | Calculates fourier transform of velocity autocorrelation function. | trajectory\_list | Tab demarkated data file |
| [vac\_function](vac_function.md) | Calculates velocity autocorrelation function by taking two derivatives of the mean-square-displacement. | trajectory\_list | Tab demarkated data file |
| [vhd](vhd.md) | Calculates distinct part of the Van Hove function. | trajectory\_list |  Tab demarkated data file |
| [vhs](vhs.md) | Calculates self-part of the Van Hove function.| trajectory\_list | Tab demarkated data file |
| [vht](vht.md) | Adds together previously computed vhd and vhs.| trajectory\_list | Tab demarkated data file |
| [write\_list\_trajectory](write_list_trajectory.md) | Writes to file an xyz trajectory from a trajectory list.| trajectory\_list| Trajectory file in wrapped xyz format |
| [write\_list\_trajectory\_full](write_list_trajectory_full.md) | Writes to file an xyz trajectory from a trajectory list.| trajectory\_list| Trajectory file in wrapped xyz format |
