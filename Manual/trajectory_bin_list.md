<h1>trajectory_bin_list</h1>

<h2>creating and deleting trajectory_bin_list objects</h2>

| Command | Effect |
|----------|----------|
| create\_bin\_list | Creates a list of all trajectories based on a binning scheme |
| remove\_bin\_list | Deletes a bin list from memory |

<h2>Analyzing trajectory_bin_list objects</h2>

**Analysis commands that target trajectory\_list objects**

The target for these analysis tools is specified in a second line that follows the line providing the analysis command line. This target line involves specifying a defined trajectory\_bin\_list and a [trajectory\_list](trajectory_list.md). Analysis will occur on a bin-by-bin basis as determined by the trajectory\_bin\_list specified. Only those particles included in the specified trajectory\_list at each time will be employed in the analysis. The syntax of the target line is as follows.

_bin\_list \<name of trajectory\_bin\_list\> \<name of trajectory_list\> <persistence check? 0 or 1 (optional)>_

Put another way, this takes the intersection of the trajectory_list and (separatedly_) each bin of the trajectory_bin_list, and performs the calculation on each bin individually. The analysis will in general output a separate tab-demarcated file for each bin, which may produce a large number of files. An ending will be appended to the name of each output file with the bin indicies (e.g. <output file>.2.1.2.bindata)

The optional argument _\<persistence check\>_ specifies whether or not, for dynamic analyses, the algorithm retains only particles that are in the bin at start and end times of a given time gap. The default is zero, meaning this check is not performed.

| Command | Effect | Output |
|----------|----------|----------|
| [composition](composition.md) | Calculates the composition and number density. | Tab demarcated data file |
| [find_edge](find_edge.md) | Finds the leading edge of a set of particles along a vector. | Tab demarcated data file |
| [isfs](isfs.md) | Calculates the self-part of the intermediate scattering function. | Tab demarcated data file |
| [mean\_displacement](mean_displacement.md) | Calculates the mean (not mean square) displacement vector of all particles as a function of time. | Tab demarcated data file |
| [msd](msd.md) | Calculates the mean-square displacement. | Tab demarcated data file |
| [msd\_2d](msd_2d.md) | Calculates the two-dimensional mean-square displacement. | Tab demarcated data file |
| [structure\_factor](structure_factor.md) | Computes the structure factor. | Tab demarcated data file |
| [trajectory\_list\_decay](trajectory_list_decay.md) | Computes an autocorrelation function for inclusion in the trajectory\_list. | Tab demarcated data file |
| [unsteady\_velocity](unsteady_velocity.md) | Calculates mean velocity as a function of time, provided that velocities were read in at file read time. | Tab demarcated data file. |
| [write\_bin\_xyz](write_bin_xyz.md) | Writes xyz files for all bins. | Series of xyz trajectory files |]
