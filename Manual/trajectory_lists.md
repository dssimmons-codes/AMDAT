<h1>trajectory_list</h1>

trajectory\_lists are the foundational data object in AMDAT. A trajectory\_list stores a list of particles for analysis, which may (or may not) vary with the time frame in the trajectory. Some trajectory\_lists are static and contain a fixed list of particles for analysis. Some are dynamic and may contain different particles at each time step or in each time block (in exponential or blocked time schemes). 

<h2>Creating and manipulating trajectory_lists</h2>

There are many methods of creating _trajectory\_list_ objects. The simplest way to create a static trajectory list is the _[create\_list](create_list.md)_ command, which allow the user to select particles based on molecular species, atom types, or other characteristics within the chemical structure of the system. More complex methods may perform some analysis on a trajectory list or some other data object and generate a dynamical trajectory\_list based on some analysis of the system. A table of commands that create or destroy trajectory\_list objects follows.

**Commands that create or destroy trajectory\_list objects**

| Command | Effect |
|----------|----------|
| [combine\_trajectories](combine_trajectories.md)    | Creates a new list of trajectories by merging multiple existing of existing trajectory\_lists. |
| [create\_list](create_list.md)    | Creates a static list of particle trajectories. |
| [delete\_trajectory_list](delete_trajectories.md) | Deletes an existing trajectory\_list, freeing up the associated memory. |
| [flatten\_multibodies](flatten_multibodies.md)    | Creates a trajectory list comprised of all trajectories in a set of multibodies. |
| [invert\_list](invert_list.md)    | Creates a new trajectory list by inverting an existing trajectory list and then intersecting it with a second list. |
| [thresholded\_list](thresholded_list.md)    | Creates a  list of trajectories by applying a value threshold to a specified value\_list. |

<h2>Analyzing trajectory_lists</h2>

Most trajectory analysis methods consist of two lines, where the first line (the "type" line) defines the type of analysis to be performed and the output file, and the second line (the "target" line) specifies the atoms or molecules on which the analysis is to be performed. The target line generally has the syntax

_list \<name of trajectory\_list to be analyzed\>_

Some trajectory analysis methods must be provided two targets, which are specified via two consecutive lines below the command itself. The overall structure is then as follows.

_\<Analysis command\> \<keywords\>_ <br>
_list <name of trajectory\_list 1 for analysis>_ <br>
_list <name of trajectory\_list 2 for analysis>_

A large number of analysis tools target trajectory\_list objects. 

**Analysis commands that target trajectory\_list objects**

| Command | Effect | Output |
|----------|----------|----------|
| [clustered\_list](clustered_list_decay.md) | Generates a new trajectory list by finding clustered particles in the original list. | New trajectory\_list |
| [compare_gaussian](compare_gaussian.md) | Finds points where the self Van Hove crosses the Gaussian approximation to the self Van Hove. | Tab demarkated data file |
| [composition](composition.md) | Calculates the composition and number density. | Tab demarkated data file |
| create_distance_neighborlist | Constructs a neighborlist object based on a distance cutoff criterion | neighbor_list |
| create_voronoi_neighborlist | Constructs a neighborlist object based on the voronoi tesselation | neighbor_list |
| displacement_list | Computes the displacement of each particle at some time separation and stores the results in a value_list | Tab demarkated data file and value_list |
| [displacement_map](displacement_map.md) | Generates maps of particle displacement magnitude and writes to pdb files with value beta specifying displacement.| PDB file |
| [isf\_list](isf_list.md) | Calculate full intermediate scattering function | Tab demarkated data file |
| [isfs](isfs.md) | Calculates the self-part of the intermediate scattering function. | Tab demarkated data file |
| [mean\_closest\_distance](mean_closest_distance.md) | Calculates the average distance between particles in two distinct trajectory\_lists. | Tab demarkated data file |
| [msd](msd.md) | Calculates the mean-square displacement. | Tab demarkated data file |
| [msd\_2d](msd_2d.md) | Calculates the two-dimensional mean-square displacement. | Tab demarkated data file |
| [n\_fold](n_fold.md) | Calculates mean and time dependence of n-fold order parameter and writes maps to files. | Tab demarkated data file and a PDB file |
| [ngp](ngp.md) | Calculates non-Gaussian parameter of the mean squared displacement. | Tab demarkated data file |
| [rdf](rdf.md) | Computes the radial distribution function. | Tab demarkated data file |
| [rgtensor\_stats](rgtensor_stats.md) | Calculates statistics describing the radius of gyration tensor characterizing particle trajectories as a function of time. | Tab demarkated data file |
| [stiffness\_dist](stiffness_dist.md) | Calculates distribution of inverse Debye-Waller factor values 1/u2. | Tab demarkated data file |
| [strings](strings.md) | Analyzes 'stringlike' collective motion. Planned for deprecation in favor of string_multibodies. | Tab demarkated data file. |
| string_multibodies | Builds multibodies comprised of 'stringlike' collective displacements | multibody_list |
| [structure\_factor](structure_factor.md) | Computes the structure factor. | Tab demarkated data file |
| [trajectory\_list\_decay](trajectory_list_decay.md) | Computes an autocorrelation function for inclusion in the trajectory\_list. | Tab demarkated data file |
| [u2dist](u2dist.md) | Calculates distribution of square displacements at a specified time. | Tab demarkated data file |
| [vac\_fourier](vac_fourier.md) | Calculates fourier transform of velocity autocorrelation function. | Tab demarkated data file |
| [vac\_function](vac_function.md) | Calculates velocity autocorrelation function by taking two derivatives of the mean-square-displacement. | Tab demarkated data file |
| [vhd](vhd.md) | Calculates distinct part of the Van Hove function.| Tab demarkated data file |
| [vhs](vhs.md) | Calculates self-part of the Van Hove function. | Tab demarkated data file |
| [vht](vht.md) | Adds together previously computed vhd and vhs. | Tab demarkated data file |
| [write\_list\_trajectory](write_list_trajectory.md) | Writes to file an xyz trajectory from a trajectory list.| Trajectory file in wrapped xyz format |
| [write\_list\_trajectory\_full](write_list_trajectory_full.md) | Writes to file an xyz trajectory from a trajectory list.| Trajectory file in wrapped xyz format |
| find_fast | Deprecated; does not work properly. In principle find particles that are 'fast' based on the above comparison between the self Van Hove and Gaussian approximation.| - |
| gyration_radius | Deprecated; to be removed. | - |
| isf | Deprecated; to be removed. | - |

