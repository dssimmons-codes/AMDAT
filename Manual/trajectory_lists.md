<h1>trajectory_list</h1>

`trajectory_list`s are the foundational data object in AMDAT. A `trajectory_list` stores a list of particles for analysis, which may (or may not) vary with the time frame in the trajectory. Some `trajectory_list`s are static and contain a fixed list of particles for analysis. Some are dynamic and may contain different particles at each time step or in each time block (in exponential or blocked time schemes). 

<h2>Creating and manipulating trajectory lists</h2>

There are many methods of creating `trajectory_list` objects. The simplest way to create a static trajectory list is the [`create_list`](create_list.md) command, which allow the user to select particles based on molecular species, atom types, or other characteristics within the chemical structure of the system. More complex methods may perform some analysis on a trajectory list or some other data object and generate a dynamical `trajectory_list` based on some analysis of the system. A table of commands that create or destroy `trajectory_list` objects follows.

**Commands that create or destroy `trajectory_list` objects**

| Command | Effect |
|----------|----------|
| [`combine_trajectories`](combine_trajectories.md)    | Creates a new list of trajectories by merging multiple existing of existing `trajectory_list` objects. |
| [`create_list`](create_list.md)    | Creates a static list of particle trajectories. |
| [`delete_trajectory_list`](delete_trajectories.md) | Deletes an existing `trajectory_list`, freeing up the associated memory. |
| [`flatten_multibodies`](flatten_multibodies.md)    | Creates a trajectory list comprised of all trajectories in a set of multibodies. |
| [`invert_list`](invert_list.md)    | Creates a new trajectory list by inverting an existing trajectory list and then intersecting it with a second list. |
| [`thresholded_list`](thresholded_list.md)    | Creates a  list of trajectories by applying a value threshold to a specified `value_list`. |
| [`traj_list_from_bin_list`](traj_list_from_bin_list.md) | Creates a `trajectory_list` object from a bin of a [`trajectory_bin_list`](trajectory_bin_list.md) |

In addition to the above, all commands that create multibodies also creates an associated `trajectory_list` comprised of either the centroids or centers of mass of the multibodies, as specified by the user. See the [multibodies](multibodies.md) page for more information.

<h2>Analyzing trajectory lists</h2>

Most trajectory analysis methods consist of two lines, where the first line (the "type" line) defines the type of analysis to be performed and the output file, and the second line (the "target" line) specifies the atoms or molecules on which the analysis is to be performed. The target line generally has the syntax

`list <name of trajectory_list to be analyzed>`

Some trajectory analysis methods must be provided two targets, which are specified via two consecutive lines below the command itself. The overall structure is then as follows.

```
<Analysis command> <keywords>
list <name of trajectory_list 1 for analysis>
list <name of trajectory_list 2 for analysis>
```

A large number of analysis tools target `trajectory_list` objects. 

**Analysis commands that target `trajectory_list` objects**

| Command | Effect | Output |
|----------|----------|----------|
| [`clustered_list`](clustered_list_decay.md) | Generates a new trajectory list by finding clustered particles in the original list. | New `trajectory_list` |
| [`compare_gaussian`](compare_gaussian.md) | Finds points where the self Van Hove crosses the Gaussian approximation to the self Van Hove. | Tab demarkated data file |
| [`composition`](composition.md) | Calculates the composition and number density. | Tab demarkated data file |
| [`composition_vs_time`](composition_vs_time.md) | Calculates and outputs the composition and number density as a function of time. | Tab demarkated data file |
| [`create_distance_neighborlist`](create_distance_neighborlist.md) | Constructs a neighborlist object based on a distance cutoff criterion | `neighbor_list` |
| [`create_voronoi_neighborlist`](create_voronoi_neighborlist.md) | Constructs a neighborlist object based on the voronoi tesselation | `neighbor_list` |
| [`displacement_list`](displacement_list.md) | Computes the displacement of each particle at some time separation and stores the results in a `value_list` | Tab demarkated data file and `value_list` |
| [`displacement_dist`](displacement_dist.md) | Computes the displacement of each particle at some time separation and reports a distribution | Tab demarcated data file. |
| [`displacement_map`](displacement_map.md) | Generates maps of particle displacement magnitude and writes to pdb files with value beta specifying displacement.| PDB file |
| [`incremental_mean_displacement`](incremental_mean_displacement.md) | Calculates mean displacement of particles as a function of time. | Tab demarcated data file. |
| [`isf_list`](isf_list.md) | Calculate full intermediate scattering function | Tab demarkated data file |
| [`isfs`](isfs.md) | Calculates the self-part of the intermediate scattering function. | Tab demarkated data file |
| [`mean_closest_distance`](mean_closest_distance.md) | Calculates the average distance between particles in two distinct `trajectory_list`s. | Tab demarkated data file |
| [`mean_displacement`](mean_displacement.md) | Calculates the mean (not mean square) displacement vector of all particles as a function of time. | Tab demarcated data file |
| [`msd`](msd.md) | Calculates the mean-square displacement. | Tab demarkated data file |
| [`msd_2d`](msd_2d.md) | Calculates the two-dimensional mean-square displacement. | Tab demarkated data file |
| [`n_fold`](n_fold.md) | Calculates mean and time dependence of n-fold order parameter and writes maps to files. | Tab demarkated data file and a PDB file |
| [`ngp`](ngp.md) | Calculates non-Gaussian parameter of the mean squared displacement. | Tab demarkated data file |
| [`rnf`](radial_count.md) | Computes the mean number of particles in shells around a set of central particles | Tab demarcated data file |
| [`rdf`](rdf.md) | Computes the radial distribution function. | Tab demarkated data file |
| [`rgtensor_stats`](rgtensor_stats.md) | Calculates statistics describing the radius of gyration tensor characterizing particle trajectories as a function of time. | Tab demarkated data file |
| [`stiffness_dist`](stiffness_dist.md) | Calculates distribution of inverse Debye-Waller factor values 1/u2. | Tab demarkated data file |
| [`strings`](strings.md) | Analyzes 'stringlike' collective motion. Planned for deprecation in favor of string_multibodies. | Tab demarkated data file. |
| `streamlined_strings` | Functionality unknown. Requires additional investigation. | Tab demarkated data file. |
| `string_multibodies` | Builds multibodies comprised of 'stringlike' collective displacements | `multibody_list` |
| [`structure_factor`](structure_factor.md) | Computes the structure factor. | Tab demarkated data file |
| [`trajectory_list_decay`](trajectory_list_decay.md) | Computes an autocorrelation function for inclusion in the `trajectory_list`. | Tab demarkated data file |
| [`u2dist`](u2dist.md) | Calculates distribution of square displacements at a specified time. | Tab demarkated data file |
| [`unsteady_velocity`](unsteady_velocity.md) | Calculates mean velocity as a function of time, provided that velocities were read in at file read time. | Tab demarcated data file. |
| [`vac_fourier`](vac_fourier.md) | Calculates fourier transform of velocity autocorrelation function. | Tab demarkated data file |
| [`vac_function`](vac_function.md) | Calculates velocity autocorrelation function by taking two derivatives of the mean-square-displacement. | Tab demarkated data file |
| [`vhd`](vhd.md) | Calculates distinct part of the Van Hove function.| Tab demarkated data file |
| [`vhs`](vhs.md) | Calculates self-part of the Van Hove function. | Tab demarkated data file |
| [`vht`](vht.md) | Adds together previously computed vhd and vhs. | Tab demarkated data file |
| [`write_list_trajectory`](write_list_trajectory.md) | Writes to file an xyz trajectory from a trajectory list.| Trajectory file in wrapped xyz format |
| [`write_list_trajectory_full`](write_list_trajectory_full.md) | Writes to file an xyz trajectory from a trajectory list.| Trajectory file in wrapped xyz format |
| [`write_single_particle`](write_single_particle.md) | Writes an xyz file from a single particle, in which its locations at distinct times are written out at a single time. | xyz trajectory file |
| `find_fast` | Deprecated; does not work properly. In principle find particles that are 'fast' based on the above comparison between the self Van Hove and Gaussian approximation.| - |
| `find_fast_fixedthreshold` | Finds fast particles based on displacements higher than some value. Likely to become deprecated in favor of valuelist approach. | new `trajectory_list` |
| `find_between` | Finds particles with displacements in some value range. Likely to become deprecated in favor of valuelist approach. | new `trajectory_list` |
| `gyration_radius` | Deprecated; to be removed. | - |
| `isf` | Deprecated; to be removed. | - |
| `radial_debye_waller` | Computes dwf as a function of distance from (0,0,0). Planned for deprecation.|
| `vector_autocorrelation_function` | Computes autocorrelation function for a specified set of vectors. Essentially replaced by [`raf`](raf.md) [`multibody`](multibodies.md) functionality and planned for deprecation. | Tab-demarcated data file. |


Finally, there are a small number of analysis tools that do not directly run on `trajectory_list`s, but output some transformation of a prior analysis performed on a trajectory list.

| Command | Effect | Output |
|----------|----------|----------|
| `isfd` | Deprecated and likely incorrect; to be removed. | - |
| `structure_factor_from_rdf` | Performs a fourier transform on data from a prior rdf calculation to arrive at a structure factor. | Tab demarcated data file |

