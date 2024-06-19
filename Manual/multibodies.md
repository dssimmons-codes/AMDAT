<h1>multibodies.md</h1>

A multibody is an entity that consists of multiple particles (each with its own trajectory). Multibodies are tracked and analyzed via multibody_lists, each of which contains a list of many multibodies that may change with time. These objects allow for analysis of internal degrees of freedom of multiparticle entities. In general, no multibodies are defined at the beginning of an instance of AMDAT. This is distinct from the case of [trajectory_lists], where the underlying trajectories in many cases already exist at the beginning of an AMDAT run (since they consist of particles created during trajectory file read-in). Therefore, some methods create multibodies themselves, which are organized into multibody_lists for analysis. Others merely create new multibody_lists storing existing multibodies.

<h2>Commands that create multibodies and multibody_lists</h2>


| Command | Effect |
|----------|----------|
| create\_multibodies | Creates multibodies based upon specified locations within molecules |
| string\_multibodies | Builds multibodies comprised of 'stringlike' collective displacements | 


<h2>Commands that create new multibody_lists from existing multibodies or delete multibody_lists</h2>

| Command | Effect |
|----------|----------|
| combine\_multibody\_lists | Creates a new multibody_list by merging a set of existing multibody\_lists |
| region\_multibody\_list | Creates a new  multibody list by sorting an existing multibody_list into spatial bins based on their centroid or center of mass location |
| threshold\_multibody\_list| Creates a new multibody list by downselecting a multibody list based on the number of bodies in each multibody. |
| delete\_multibody\_list | Deletes an existing multibody\_list and frees the associated memory. Does not delete the underlying multibodies themselves. |

<h2>Analyzing multibody_lists</h2>

The target multibody_list for these analysis tools is generally provided as an argument within the command line itself.

| Command | Effect | Output |
|----------|----------|----------|
| size\_statistics | Reports statistics on the distribution of the number of bodies in multibodies in a given multibody\_list. | Tab-demarcated data file |
| raf | Computes the reorientation autocorrelation function for a list of 2-body multibodies. | Tab-demarcated data file |
| baf | Deprecated and planned for removal. | - |
| orientational\_correlation | Calculates the orientational correlation of 2-body multibodies with an externally defined vector. | Tab-demarcated data file |
| gyration\_radius | Computes the mean multibody gyration radius of for the multibodies in a specified multibody_list. | Tab-demarcated data file |
