<h1>create_list</h1>
Creates a static list of particle trajectories for analysis based upon set features such as the particle's type and location within molecules of a given species.

```
create_list <listname>

<keyword> <arguments>
```

Available keywords and associated arguments are listed below. It is possible to combine multiple keywords and associated arguments to enable selection of complex sets, as follows:

```
create_list <listname>

<keyword1> <arguments1> <keyword2> <arguments2>
```

| Keyword             | Arguments                                                   | Description |
|---------------------|-------------------------------------------------------------|-------------|
| ``all``             | ``none``                                                    | Selects all atoms in the system. |
| ``type_system``     | ``<atom_type>``                                              | Selects atoms of type ``<atom_type>``. |
| ``type_species``    | ``<species_name> <atom_type>``                               | Selects atoms of type ``<atom_type>`` in species ``<species_name>``. |
| ``species``         | ``<species_name>``                                           | Selects all atoms within species ``<species_name>``. |
| ``molecule``        | ``<species_name> <molecule_index>``                          | Selects atoms within a single molecule of species ``<species_name>``, where the molecule is specified by ``<molecule_index>`` (index within the set of molecules of that species in the trajectory file, starting at 0). |
| ``type_molecule``   | ``<species_name> <molecule_index> <atom_type>``              | Selects atoms of type ``<atom_type>`` within a single molecule of species ``<species_name>``, specified by ``<molecule_index>`` (index within the set of molecules of that species in the trajectory file, starting at 0). |
| ``atom_species``    | ``<species_name> <atom_type> <atom_index>``                  | Selects a single atom of type ``<atom_type>`` within every molecule of species ``<species_name>``, where ``<atom_index>`` is the index of the atom of that type within each molecule (starting from 0 based on the order in the trajectory file). To be more specific, each atom is specified by the type of the atom, and then the index of the atom within that atomtype within that species. These indices are determined by original order of read in from the trajectory file. For a given molecule, AMDAT separately indices atoms of distinct atom types. Indexing starts at 0. Therefore for an atomtype "H", for example, an index of 4 would specify the "H" atom within that molecule that was read in after four prior "H" atoms were read into that molecule. |
