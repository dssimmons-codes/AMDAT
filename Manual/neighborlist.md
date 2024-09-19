<h1>neighborlist</h1>

neighborlist objects store a list of neighbors at each time of each particle in a list of particles. Every neighborlist is also a [value_list](value_list.md), with the stored values corresponding to the number of neighbors of each particle at each time. neighborlists are thus also subject to all analysis tools that can be applied to value_lists.

<h2>Creating neighborlists</h2>

Presently there are two commands for neighborlist creation, which differ in their criteria for idenfying neighbors. For most applications, it is suggested that neighborlists be computed using all particles in the system. Use of a subset of particles may lead to strange behavior (exe. in voroni tesselation). This implies that one should typically employ a static trajectory_list in constructing neighbor lists, not one that evolves over time. However, this is not strictly required, and use of particle subsets may perform reasonably with distance-based neighborlist building. There are also cases in which partial particle sets may be desired for voronoi neighbor list building (e.g. if one wishes to analyze the neighbor structures of composite objects, such as atomistic nanoparticles, at a center of mass level) Neighborlist analysis tools (discussed further below) allow for downselection of which particles to be analyzed at a later time, so even if neighborlists are build with all particles in the system, they can be downselected later.

**Commands that create or remove neighborlists**

| Command | Effect |
|----------|----------|
| [create\_distance\_neighborlist](create_distance_neighborlist.md) | Constructs a neighborlist object based on a distance cutoff criterion. | 
| [create\_voronoi\_neighborlist](create_voronoi_neighborlist.md) | Constructs a neighborlist object based on the voronoi tesselation. |
| [delete\_neighborlist](delete_neighborlist.md) | Removes a neighborlist from memory. |

<h2>Analyzing neighborlists</h2>

In general, the target neighborlist for analysis is provided as an argument within the analysis command line as specified in each analysis command page. In some cases a trajectory list for downselection of particles to be included may be input as an additional target line or in the command line. See the page for each analysis tool for details.

**Commands that analyze neighborlists**

| Command | Effect | Output |
|----------|----------|----------|
| [neighbor_decorrelation_function](neighbor_decorrelation_function.md) | compute an autocorrelation function for the decay in number of conserved neighbors | Tab-demarcated data file |
| [persistent\_neighbors](persistent_neighbors.md) | Constructs multibodies comprised of particles that remain neighbors after a specified time gap | Multibody_list|
| [value\_statistics](value_statistics.md) | Computes statistics on the distribution of neighbors per particle. | Tab-demarcated data file |
