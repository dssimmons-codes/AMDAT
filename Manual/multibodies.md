<h1>multibodies</h1>

A multibody is an entity that consists of multiple particles (each with its own trajectory). Multibodies are tracked and analyzed via `multibody_list`s, each of which contains a list of many multibodies that may change with time. These objects allow for analysis of internal degrees of freedom of multiparticle entities. In general, no multibodies are defined at the beginning of an instance of AMDAT. This is distinct from the case of [`trajectory_list`](trajectory_list.md), where the underlying trajectories in many cases already exist at the beginning of an AMDAT run (since they consist of particles created during trajectory file read-in). Therefore, some methods create multibodies themselves, which are organized into `multibody_list`s for analysis. Others merely create new `multibody_list`s storing existing multibodies.

<h2>Commands that create multibodies and multibody lists</h2>

| Command | Effect |
|----------|----------|
| [`create_multibodies`](create_multibodies.md) | Creates multibodies based upon specified locations within molecules |
| [`comover_multibodies`](comover_multibodies.md) | Creates multibodies comprised of particles that remain neighbors based on some distance cutoff. (In beta, not ready for general use, may be removed.) |
| [`relative_displacement_strings`](relative_displacement_strings.md) | Creates multibodies comprised of particles that move together. (In beta, not ready for general use, may be removed.) |
| [`string_multibodies`](string_multibodies.md) | Builds multibodies comprised of 'stringlike' collective displacements | 


<h2>Commands that create new multibody lists from existing multibodies or delete multibody lists</h2>

| Command | Effect |
|----------|----------|
| [`combine_multibody_lists`](combine_multibody_lists.md) | Creates a new `multibody_list` by merging a set of existing `multibody_list`s |
| [`region_multibody_list`](region_multibody_list.md) | Creates a new  `multibody_list` by sorting an existing `multibody_list` into spatial bins based on their centroid or center of mass location |
| [`threshold_multibody_list`](threshold_multibody_list.md) | Creates a new `multibody_list` by downselecting a `multibody_list` based on the number of bodies in each multibody. |
| [`delete_multibody_list`](delete_multibody_list.md) | Deletes an existing `multibody_list` and frees the associated memory. Does not delete the underlying multibodies themselves. |

<h2>Analyzing multibody lists</h2>

The target `multibody_list` for these analysis tools is generally provided as an argument within the command line itself.

| Command | Effect | Output |
|----------|----------|----------|
| [`flatten_multibodies`](flatten_multibodies.md) | Generates a `trajectory_list` containing all particles in multibodies in the specified `multibody_list`. | new `trajectory_list` |
| [`gyration_radius`](gyration_radius.md) | Computes the mean multibody gyration radius of for the multibodies in a specified `multibody_list`. | Tab-demarcated data file |
| [`orientational_correlation`](orientational_correlation.md) | Calculates the orientational correlation of 2-body multibodies with an externally defined vector. | Tab-demarcated data file |
| [`raf`](raf.md) | Computes the reorientation autocorrelation function for a list of 2-body multibodies. | Tab-demarcated data file |
| [`size_statistics`](size_statistics.md) | Reports statistics on the distribution of the number of bodies in multibodies in a given `multibody_list`. | Tab-demarcated data file |
| `baf` | Deprecated and planned for removal. | - |
