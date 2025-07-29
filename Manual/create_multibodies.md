<h1>create_multibodies</h1>

Creates multibodies, and an associated `multibody_list` for later access, based on specified structural location within molecules. Additionally creates a new `trajectory_list` of the same name, which stores a list of trajectories of either the centroids or centers of mass of the multibodies, as specified by the user.

`create_multibodies <name of multibody list to create> <name of atom type to be assigned to center of mass or centroid> <either "centroid" or "com"> <keywords (see below)>`

The `centroid` or `com` keywords specify whether to employ the centroid or the center of mass to compute a trajectory characteristic of the collective translational motion of each multibody. These trajectories are accesible later via a trajectory list of the same name specified for the name of the multibody list itself.

Allowable keywords for atom selection, and the options required for each keyword, are as follows.

`all_molecule`

Creates a list of multibodies in which each molecule in the system is converted to a multibody containing all atoms in that molecule.

`species_molecule <species_name>`

Creates a list of multibodies in which each molecule of species `<species_name>` in the system is converted to a multibody containing all atoms in that molecule.

`species_type <species_name> <type_name>`

Creates a list of multibodies in which each multibody is generated from each molecule of species `<species_name>`, with each of these multibodies containing all atoms of type `<type_name>` in a given molecule of that species.

`species_atomlist <species_name> <type_name1> <index1> <type_name2> <index2> ... <type_namelast> <indexlast>`

Creates a list of multibodies in which each multibody is generated from each molecule of species `<species_name>`, with each of these multibodies containing specified by a series of `atom_type`, `atom_index` pairs. Each atom is specified by the type of the atom, and then the index of the atom within that atomtype within that species. These indices are determined by original order of read in from the trajectory file. For a given molecule, AMDAT separately indices atoms of distinct atom types. Indexing starts at 0. Therefore for an atomtype "H", for example, an index of 4 would specify the "H" atom within that molecule that was read in after four prior "H" atoms were read into that molecule.
