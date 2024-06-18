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

| Command | Effect |
|----------|----------|
| composition | Calculates the composition and number density. |
| isfs | Calculates the self-part of the intermediate scattering function. |
| msf | Calculates the mean-square displacement. |
| structure\_factor | Computes the structure factor. |
| trajectory\_list\_decay | Computes an autocorrelation function for inclusion in the trajectory\_list. |
| clustered\_list | Generates a new trajectory list by finding clustered particles in the original list. |
| compare_gaussian | Finds points where the self Van Hove crosses the Gaussian approximation to the self Van Hove. |
| displacement_map | Generates maps of particle displacement magnitude and writes to pdb files with value beta specifying displacement.|
| find_fast | Find particles that are 'fast' based on the above comparison between the self Van Hove and Gaussian approximation.|
| gyration_radius | Deprecated; to be removed. |
| isf | Deprecated; to be removed. |
| isfs\_list | Calculate full intermediate scattering function |
| mean\_closest\_distance | Calculates the average distance between particles in two distinct trajectory\_lists. |
| n\_fold | Calculates mean and time dependence of n-fold order parameter and writes maps to files. |
| ngp | Calculates non-Gaussian parameter of the mean squared displacement. |
| rdf | Computes the radial distribution function. |
| rgtensor\_stats | Calculates statistics describing the radius of gyration tensor characterizing particle trajectories as a function of time. |
| stiffness\_dist | Calculates distribution of inverse Debye-Waller factor values 1/u2. |
| strings | Analyzes 'stringlike' collective motion |
| u2dist | Calculates distribution of square displacements at a specified time. |
| vac\_fourier | Calculates fourier transform of velocity autocorrelation. |
| vac\_function | Calculates velocity autocorrelation function by taking two derivatives of the mean-square-displacement. |
| vhd | Calculates distinct part of the Van Hove function.|
| vhs | Calculates self-part of the Van Hove function. |
| vht | Adds together previously computed vhd and vhs. |
| write\_list\_trajectory | Writes to file an xyz trajectory from a trajectory list.|
| write\_list\_trajectory\_full | Writes to file an xyz trajectory from a trajectory list.|

