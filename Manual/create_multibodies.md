<h1>create_multibodies</h1>

Creates multibodies, and an associated multibody\_list for later access, based on specified structural location within molecules. Additionally creates a new trajectory\_list of the same name, which stores a list of trajectories of either the centroids or centers of mass of the multibodies, as specified by the user.

_create\_multibodies \<name of multibody list to create\> \<name of atom type to be assigned to center of mass or centroid\> \<either "centroid" or "com"\> <keywords (see below)_

The "centroid" or "com" keywords specify whether to employ the centroid or the center of mass to compute a trajectory characteristic of the collective translational motion of each multibody. These trajectories are accesible later via a trajectory list of the same name specified for the name of the multibody list itself.

Allowable keywords for atom selection, and the options required for each keyword, are as follows.

_all\_molecule_

Creates a list of multibodies in which each molecule in the system is converted to a multibody containing all atoms in that molecule.

_species\_molecule \<species\_name\>_

Creates a list of multibodies in which each molecule of species _<species\_name>_ in the system is converted to a multibody containing all atoms in that molecule.

_species\_type \<species\_name\> \<type\_name\>_

Creates a list of multibodies in which each multibody is generated from each molecule of species _<species\_name>_, with each of these multibodies containing all atoms of type _<type\_name>_ in a given molecule of that species.

_species\_atomlist \<species\_name\> \<type\_name1\> \<index1\> \<type\_name2\> \<index2\> ... \<type\_namelast\> \<indexlast\>_

Creates a list of multibodies in which each multibody is generated from each molecule of species _<species\_name>_, with each of these multibodies containing specified by a series of atom_\type atom_\index pairs. Each atom is specified by the type of the atom, and then the index of the atom within that atomtype within that species. These indices are determined by original order of read in from the trajectory file. For a given molecule, AMDAT separately indices atoms of distinct atom types. Indexing starts at 0. Therefore for an atomtype "H", for example, an index of 4 would specify the "H" atom within that molecule that was read in after four prior "H" atoms were read into that molecule.
