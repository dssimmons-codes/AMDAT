<h1>Value_List</h1>

A value\_list object stores a list of numbers for some set of particles as a function of time. Value\_lists can generally be created in two ways - during read-in of a custom trajectory file using the [custom_manual](read_custom_manual.md) read in method, or by an analysis performed during the LAMMPS run. For information on creation of value_lists at trajectory readin, see the page on the [custom_manual](read_custom_manual.md) read in method.

<h2>Analysis tools that create value_list objects</h2>

To be compiled. For now, see the page on [trajectory_lists](trajectory_lists.md), where most or all of these analysis tools are documented.

<h2>Analysis tools that operate on value_list objects</h2>

| Command | Effect | Output |
|----------|----------|----------|
| autocorrelate_value_list | Calculates the time autocorrelation of a value_list. | Tab-demarcated data file|
| crosscorrelate_value_lists | Determines the mean time cross-correlation of two value lists. | Tab-demarcated data file}
| value_list (keyword threshold_value) | Constructs a new trajectory list, by selecting only those particles with values in the selected value list that are in a specified value range. | a new trajectory\_list |
| value_list (keyword threshold_percentile) | Constructs a new trajectory list, by selecting only those particles whose values in the selected value list are in a specified percentile range. | a new trajectory\_list |
| value_list (keyword write_pdb) | Writes out a .pdb file in which the value in the beta column is drawn from a value list.  | PDF file |
| value_statistics | Outputs distribution statistics on a value_list | Tab-demarcated data file |

