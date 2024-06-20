<h1>create_distance_neighborlist</h1>

Create a neighbor list, for all particles, based on a distance cutoff scheme

_create\distance\_neighborlist \<name of neighborlist to be created\> \<distance cutoff multiplier\> \<path to sigma matrix file\> \<optional: index of first timestep to build neighborlists\> <optional: index of last timestep to build neighborlists\>_

The sigma matrix file specifies a contact distance for each pair of particle types. The cutoff distance for neighborlist building is specific to each pair of particle types and is their contact distance times the cutoff distance multipler. The file should be tab demarcated in the following form.

||||||
|----------|----------|----------|----------|----------|
| Particle type name 1 | $\sigma_{1-1}$ | $\sigma_{1-2}$ | ... | $\sigma_{1-n}$ |
| Particle type name 2 | $\sigma_{2-1}$ | $\sigma_{2-2}$ | ... | $\sigma_{2-n}$ |
| ... | ... |... | ... |... |
| Particle type name n | $\sigma_{n-1}$ | $\sigma_{n-2}$ | ... | $\sigma_{n-n}$ |
