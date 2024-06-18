<h1>create_list</h1>
Creates a static list of particle trajectories for analysis based upon set features such as the particle's type and location within molecules of a given species.

_create\_list \<listname\>_

_\<keyword\> \<arguments\>_

Available keywords and associated arguments are listed below. It is possible to combine multiple keywords and associated arguments to enable selection of complex sets, as follows:

_create\_list \<listname\>_

_\<keyword1 \> \<arguments1\> \<keyword2 \> \<arguments2\>_

keywordargs

_all none_

Selects all atoms in the system

_type\_system \<atom\_type\>_

Selects atoms of type _atom\_type_

_type\_species \<species\_name\> \<atom\_type\>_

Selects atoms of type _atom\_type_ in species _species\_name_

_species \<species\_name\>_

Selects all atoms within species _species\_name_

_molecule \<species\_name\> \<molecule\_index\>_

Selects atoms within a single molecule of species _\<species\_name\>_, where the molecule is specified by _molecule\_index_, which corresponds to the position of the molecule in the trajectory file within the set of molecules of that species (with the initial molecule of a species indexed as 0.

_type\_molecule \<species\_name\> \<molecule\_index\> \<atom\_type\>_

Selects atoms of type _\<atom\_type\>_ within a single molecule of species _\<species\_name\>_, where the molecule is specified by _molecule\_index_, which corresponds to the position of the molecule in the trajectory file within the set of molecules of that species

_atom\_species \<species\_name\> \<atom\_type\> \<atom\_index\>_

Selects a single atom of type _\<atom\_type\>_ within every molecule species _\<species\_name\>_, with the atom being specified by _atom\_index_, which corresponds to the order in which the atoms of that type were first read from the trajectory file for that species (again, indexed from 0). To be more specific, each atom is specified by the type of the atom, and then the index of the atom within that atomtype within that species. These indices are determined by original order of read in from the trajectory file. For a given molecule, AMDAT separately indices atoms of distinct atom types. Indexing starts at 0. Therefore for an atomtype "H", for example, an index of 4 would specify the "H" atom within that molecule that was read in after four prior "H" atoms were read into that molecule.
