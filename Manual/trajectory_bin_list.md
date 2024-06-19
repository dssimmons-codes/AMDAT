<h1>trajectory_bin_list</h1>

<h2>creating and deleting trajectory_bin_list objects</h2>

| Command | Effect |
|----------|----------|
| create\_bin\_list | Creates a list of all trajectories based on a binning scheme |
| remove\_bin\_list | Deletes a bin list from memory |

<h2>Analyzing trajectory_bin_list objects</h2>

**Analysis commands that target trajectory\_list objects**

The target for these analysis tools is specified in a second line that follows the line providing the analysis command line. This target line involves specifying a defined trajectory\_bin\_list and a [trajectory\_list](trajectory_list.md). Analysis will occur on a bin-by-bin basis as determined by the trajectory\_bin\_list specified. Only those particles included in the specified trajectory\_list at each time will be employed in the analysis. The syntax of the target line is as follows.

_bin\_list \<name of trajectory\_bin\_list\> \<name of trajectory_list\>_

The analysis will in general output a separate tab-demarcated file for each bin, which may produce a large number of files.  

| Command | Effect | Output |
|----------|----------|----------|
| [composition](composition.md) | Calculates the composition and number density. | Tab demarcated data file |
| [isfs](isfs.md) | Calculates the self-part of the intermediate scattering function. | Tab demarcated data file |
| [msd](msd.md) | Calculates the mean-square displacement. | Tab demarcated data file |
| [msd\_2d](msd_2d.md) | Calculates the two-dimensional mean-square displacement. | Tab demarcated data file |
| [structure\_factor](structure_factor.md) | Computes the structure factor. | Tab demarcated data file |
| [trajectory\_list\_decay](trajectory_list_decay.md) | Computes an autocorrelation function for inclusion in the trajectory\_list. | Tab demarcated data file |
| find_edge | Finds the leading edge of a set of particles along a vector. | Tab demarcated data file |
