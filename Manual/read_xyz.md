<h1>read_xyz</h1>

<h2>Overview</h2>

xyz is a common trajectory file format produced by codes such as LAMMPS. It is a text file containing the positions of each particle at each time frame sequentially. Each frame is headed by two lines. The first reads “atoms” and the second contains the number of atoms in that frame. These are followed by a list of all atoms, where the first column is the atom type index, and the second, third, and fourth are the x, y, and z coordinates of that atom at that time step. _Note that AMDAT assumes that the atoms are in the same order in each frame_; if the trajectory does not conform to this rule, results will generally be incorrect. 

Since xyz files employ ‘wrapped’ coordinates, they generally face the challenge that they do not always facilitate unambiguous reconstruction of unwrapped particle trajectories. AMDAT will attempt to infer unwrapped trajectories by assuming that the shortest distance travelled by a given particle between two timesteps, considering all possible boundary crossings, is the correct one. If AMDAT encounters trajectories in which displacements approach half the box size in any dimension, it will output a warning indicating that unambiguous particle unwrapping is not possible. If you encounter this warning, it is strongly suggested that you not rely upon any dynamical analysis resulting from this trajectory. This problem can often be solved by employing more closely spaced frames or more generally employing a lammps custom trajectory with either unwrapped coordinates or image flags provided.

<h2>syntax</h2>

This command consists of multiple lines. Syntax is as follows.

_xyz_  
_<input\_file\_path> <optional: template\_file>_  
_<first\_species\_name> <number\_of\_molecules> <second\_species\_name> <number\_of\_molecules>…_  
_<first\_atom type> <second\_atom type> …_  
_<# of first type in first species> <# of second type in first species> …_  
_<# of first type in second species> <# of second type in second species>…_  
_..._  
_..._  
_<# of first type in last species>…<# of last type in last species>_  
_<Lx> <Ly> <Lz>_  
_<xlo> <xhi> <ylo> <yhi> <zlo> <zhi>_  

Minimally, <filenames> for xyz must include the path to the xyz trajectory file. If this is the only file provided, AMDAT will assume that the trajectory is ordered such that all molecules of a given species are grouped together, with the species appearing in the xyz file in the order in which they are provided in the list of species above.

Alternatively, the user may include an optional ‘template’ file, which provides the order of molecules by species type. The format is:

<species_k_name> <# of continuous molecules of species k>
<species_j_name> <# of continuous molecules of species j>
…   
Species in the template file may generally be repeated (ie, molecules of a given species need not be continuous). This permits more complex orderings of molecule types within the file.

In the last two input lines above, Li is the total length of the box in the i’th dimension; ilo is the position of the lower-valued boundary in the i’th direction, and ihi is the position of the upper-valued boundary in the i’th direction.

