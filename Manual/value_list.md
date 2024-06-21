<h1>Value_List</h1>

A value\_list object stores a list of numbers for some set of particles as a function of time. Value\_lists can generally be created in two ways - during read-in of a custom trajectory file using the [custom_manual](read_custom_manual.md) read in method, or by an analysis performed during the LAMMPS run. For information on creation of value_lists at trajectory readin, see the page on the [custom_manual](read_custom_manual.md) read in method.

<h2>Analysis tools that create or remove value_list objects</h2>

To be compiled. For now, see the page on [trajectory_lists](trajectory_lists.md), where most or all of these analysis tools are documented.

| Command | Effect |
|----------|----------|
|[delete\_valuelist](delete_valuelist.md)| Deletes a value\_list from memory.|
|[displacement\_list](displacement_list.md)| Creates a value list containing particle displacements at a specified time gap.|

<h2>Analysis tools that operate on value_list objects</h2>

| Command | Effect | Output |
|----------|----------|----------|
| [autocorrelate_value_list](autocorrelate_value_list.md) | Calculates the time autocorrelation of a value_list. | Tab-demarcated data file|
| [crosscorrelate_value_lists](crosscorrelate_value_lists.md) | Determines the mean time cross-correlation of two value lists. | Tab-demarcated data file}
| [value_list](value_list_command.md) (keyword _threshold\_value_) | Constructs a new trajectory list, by selecting only those particles with values in the selected value list that are in a specified value range. | a new trajectory\_list |
| [value_list](value_list_command.md) (keyword _threshold\_percentile_) | Constructs a new trajectory list, by selecting only those particles whose values in the selected value list are in a specified percentile range. | a new trajectory\_list |
| [value_list](value_list_command.md) (keyword _write\_pdb_) | Writes out a .pdb file in which the value in the beta column is drawn from a value list.  | PDF file |
| [value_statistics](value_statistics.md) | Outputs distribution statistics on a value_list | Tab-demarcated data file |
| [value_statistics_pertime](value_statistics_pertime.md) | Outputs value_list statistics in a time-resolved manner |Tab-demarcated data file|

