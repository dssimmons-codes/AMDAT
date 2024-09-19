<h1>custom_manual</h1>

<h2>Overview</h2>

Custom files are a type of trajectory file created by LAMMPS. Currently, AMDAT can only read custom trajectory files that are sorted:  _AMDAT assumes that the atoms are in the same order in each frame_; if the trajectory does not conform to this rule, results will generally be incorrect. As with [xyz](xyz.md) read-in,  a template file may also be specified. 

In comparison to the [custom](custom.md) read-in method, custom_manual allows the user increased control over the type of trajectory data that is read in, and also allows read-in of arbitrary data columns to be stored [value_lists](value_list.md) for additional analysis. For example, this could allow read-in of data produced by per-atom computes in LAMMPS, such as per-atom stresses.

read_manual can presently process two types of corrdinate sets: unwrapped (xu, yu, and zu), and wrapped_indexed (x, y, and z, and ix, iy, and iz)  The choice of which of these coordinate types to employ, as well as specification of any other data columns to read in as value_lists, is specified in an extra header_file that is provided to AMDAT and read in as specified in the syntax section below.

<h2>syntax</h2>

The syntax to be employed is as follows.

_custom_   
_<input\_file\_path> <header\_file\_path> <optional: template\_file>_  
_<first\_species\_name> <number\_of\_molecules> <second\_species\_name> <number\_of\_molecules> …_  
_<first\_atom\_type> <second\_atom\_type> …_  
_<# of first type in first species> <# of second type in first species> …_  
_<# of first type in second species> <# of second type in second species>…_  
_..._  
_..._  
_<# of first type in last species>…<# of last type in last species>_  

Minimally, <filenames> must include the path to the custom trajectory file and the path to a header file. The header file is structured as follows.

The first line must begin with a coordinate style keyword, which is either _unwrapped_ or _wrapped\_indexed_. This is followed, on the same line, by the column headers (the headers at the tops of the columns in the custom trajectory file) for the associated coordinate data. For style _unwrapped_, three headers must be provided: those corresponding to the columns containing the unwrapped x-coordinate, y-coordinate, and z-coordinate (in this order). For style _wrapped\_indexed_, six headers must be provided: those corresponding to the columns containing the wrapped x-coordinate, y-coordinate, and z-coordinate, and those corresponding to the image index flags in the x dimension, y dimension, and z dimension (in this order).

Following this initial line, the user may include additional lines, each of which specifies the header of an additional data column for read-in as a value_list, followed by the name of the value\_list to be created. The overall syntax for this file is thus as follows

_<coordinate\_style\_keyword> <coordinate\_headers>_  
_<optional: additional\_column\_header> <value\_list name>_  
...

Where an unlimited number of optional lines can be included to read in additional data columns.

If the user provides only thes two files, AMDAT will assume that the trajectory is ordered such that all molecules of a given species are grouped together, with the species appearing in the custom file in the order in which they are provided in the list of species above.

Alternatively, the user may include an optional ‘template’ file, which provides the order of molecules by species type. The format is:

<species_k_name> <# of continuous molecules of species k>  
<species_j_name> <# of continuous molecules of species j>  
…     

Species in the template file may generally be repeated (ie, molecules of a given species need not be continuous). This permits more complex orderings of molecule types within the file.

Since box size data is included in LAMMPS custom file output at each time step, time-dependent box dimensions are automatically read in from the custom file.
