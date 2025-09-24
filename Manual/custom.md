<h1>custom</h1>

<h2>Overview</h2>

Custom files are a type of trajectory file created by LAMMPS. Currently, AMDAT can only read custom trajectory files that are sorted:  _AMDAT assumes that the atoms are in the same order in each frame_; if the trajectory does not conform to this rule, results will generally be incorrect. As with [xyz](xyz.md) read-in,  a template file may also be specified. 

Coordinate types understood by AMDAT included wrapped (x, y, and z), scaled wrapped (xs, ys, and zs), unwrapped (xu, yu, and zu), and scaled unwrapped (xsu, ysu, and zsu). At least one complete set of these coordinate types must be present in order for AMDAT to proceed. AMDAT handles this information in the following way.

1.	If both scaled and unscaled coordinates of either wrapped or unwrapped type are provided, AMDAT defaults to reading the unscaled coordinates and ignores the scaled coordinates.
2.	If only wrapped coordinates are provided, they are read, and the code later attempts to infer unwrapped coordinates from the particle displacements between frames. However, if time between frames sufficiently large that particle displacements are exceeding 40% of the box size in any direction, this may lead to incorrect results due to unwrapping ambiguity.
3.	If only unwrapped coordinates are provided, they are read, and AMDAT determines the corresponding wrapped coordinates.
4.	If wrapped and unwrapped coordinates are provided, AMDAT reads both.
5.	If wrapped coordinates are provided, unwrapped coordinates are not provided, and the ‘box image index’ is provided, then AMDAT reads the wrapped coordinates and determines the unwrapped coordinates by adding to x, y, and z the box image index times the box size.
6.	If the box image index is provided by unwrapped coordinates are also provided, the box image index is ignored.

AMDAT also understands velocities (vx, vy, and vz), and masses (m), although these categories are not required and are only used for certain analysis techniques. An additional data category that must be included in the custom trajectory file is “type”. These columns may be in any order in the trajectory file. Custom files include the box bounds information in the header of each frame, so this information need not be provided by the user in the input file.

<h2>syntax</h2>

```
custom
<input_file_path> <optional: template_file> 
<first_species_name> <number_of_molecules> <second_species_name> <number_of_molecules> ... <last_species_name> <number_of_molecules> 
<first_atom_type>                     <second_atom_type>                      ...    <last_atom_type>
<# of first type in first species>    <# of second type in first species>     ...    <# of last type in first species> 
<# of first type in second species>   <# of second type in second species>    ...    <# of last type in second species>
...
...
<# of first type in last species>     <# of second type in last species>      ...    <# of last type in last species>
```

Minimally, ``<filenames>`` must include the path to the custom trajectory file. If this is the only file provided, AMDAT will assume that the trajectory is ordered such that all molecules of a given species are grouped together, with the species appearing in the `custom` file in the order in which they are provided in the list of species above.

Alternatively, the user may include an optional ‘template’ file, which provides the order of molecules by species type. The format is:

```
<species_k_name> <# of continuous molecules of species k>
<species_j_name> <# of continuous molecules of species j>
...
```
 
Species in the template file may generally be repeated (ie, molecules of a given species need not be continuous). This permits more complex orderings of molecule types within the file.

Since box size data is included in LAMMPS custom file output at each time step, time-dependent box dimensions are automatically read in from the custom file.
