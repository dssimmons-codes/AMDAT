<h1>Overview of AMDAT</h1>

AMDAT is a toolkit of molecular dynamics analysis methods. It currently has a focus on analysis of non-crystalline materials, although it contains some analyses relevant to crystalline systems (and this focus may evolve moving forward). It is designed to be a high-memory fast-compute analysis package. This means it loads the entire trajectory into active memory for analysis. This makes it extremely fast, because it does not need to read hard drive storage during operation. However, it also means that in practice memory requirements are typically 2-3 times (or more, depending on the types of analyses performed) the trajectory file size. This may become prohibitive for analysis of extremely large trajectory files.

AMDAT can read multiple file formats produced by the LAMMPS molecular dynamics simulation software, and it also can read .xtc format trajectory produced by Gromacs (although this functionality was last tested many version of GROMACS ago). The type of trajectory file employed affects the type of analyses that can be performed in AMDAT. The most versatile functionality is presently achievable with LAMMPS custom trajectory files. AMDAT enables [versatile selection of sets of particles](trajectory_lists.md) from within a simulation and contains methods of characterizing both the structure and dynamics of the particles. It enables user-definition of '[multibodies](multibodies.md)' consisting of multiple particles within a molecule and analysis of either the center of mass motion of these multibodies or of multibody correlations wotjom these multibodies. It allows for definition of [neighborlists](neighborlist.md) that track the neighbors of each particle. It allows read-in-from-trajectory-file or live calculation of [ancillary sets of values](value_list.md) on a per-particle-per-time basis and their analysis. AMDAT is currently in a beta testing phase, and functionality is being continuously improved or updated. That being said, AMDAT is written in an object oriented manner in C++, and it is designed to enable straightforward extension with new analysis techniques. We welcome submissions of new additions to this code. Some guidance on how to go about this will be provided in the developer section of this manual (to be added in the near future).

A central feature of AMDAT is its ability to read in and analyze blocked, exponentially-spaced trajectories. This is essential for analysis for dynamics over many orders of magnitude in time. The handling of this is currently documented in the [input_files](input_file.md) documentation. Both the documentation and the underlying functionality for this are expected to be expanded in the near future.

AMDAT includes no built-in tools for visualization of either trajectories or properties of trajectories. All results are output in a text-based form and must be viewed in another software (exe python, julia, excel, matlab,...). AMDAT does not include the ability to aggregate statistics acquired based on multiple simulations; this requires use of post-analysis software, such as python, julia, matlab, igor, excel, or similar. AMDAT presently does not include multithreading, although this capability is planned for introduction in a future update.

<h2>Making AMDAT</h2>

For instructions on making AMDAT, see the page on [Making AMDAT](making_amdat.md)

<h2>General concepts and terminology</h2>

AMDAT is run from the command line (see section below on Running AMDAT). The user provides an [input file](input_file.md), which provides information on the trajectory file and system metadata and provides a script specifying analyses to be performed. AMDAT reads in the specified simulation trajectory, and based on this constructs a database of particle trajectories (corresponding the the particles in the simulation trajectory). In order to perform analyses on these trajectories, the user selects one or more lists of trajectories ([trajectory_lists](trajectory_lists.md)) specifying a set of particle trajectories for analysis. Trajectories can also be selected and analyzed in a spatially-resolved manner using [trajectory_bin_lists](trajectory_bin_list.md) The user can also define [multibodies](multibodies.md), which are entities comprised of multiple particle trajectories (for example, one could define a multibody corresponding to a whole molecule or a single sidegroup). They may also define [neighbor_lists](neighbor_list.md) and [value_lists](value_list.md). Distinct sets of analysis tools can be applied to each of these data objects and are described in their respective pages. Many of these tools interact. For example, a number of analysis methods themselves generate trajectory lists, and the set of particle trajectories included in these lists may in general vary over the duration of the simulation trajectory (for example, a trajectory list of the most mobile particles – see analysis method find\_fast – will vary over the course of the simulation as particles become more or less mobile with time).  In a similar manner to trajectories, the user may define lists of these multibodies for analysis (for example, calculation of the multibody gyration radius).

**Simulation Trajectory** : The time-series of configurations being analyzed

**Frame** : The index of a given configuration in a trajectory: frame n is the nth configuration in the trajectory being analyzed. In keeping with C++ conventions, indices throughout are understood to begin at 0.

**Time** : refers to the time corresponding to a given frame

**Lag time or displacement time** or **timegap** or **time spacing** : refers to the time in displacement and time correlation calculations, in which time refers to the time from some arbitrary time origin that may be chosen at many points along the trajectory.

<h2>Running AMDAT</h2>

AMDAT is run from the command line, specifying the filename of a user [input file](input_file.md):

_./AMDAT \<FLAGS\>_

\<FLAGS\> refers to a series of options built in to AMDAT. These flags can given in any combination in any order. Below is a list of flags and their respective syntax.

_-i \<INPUT FILE\>_

\<INPUTFILE\> is the path to a file directing AMDAT to a trajectory file for analysis, providing some metadata on the trajectory, and specifying the analysis to be performed. A more detailed description can be found in section 1. \<INPUTFILE\> is required to run the program.

_-n \<# OF PROCESSORS\>_

Presently, no analysis methods allow for multithreading, so this option does nothing. Eventually, this will permit AMDAT to run some analyses on \<# OF PROCESSORS\> processing cores. Defaults to 1 if not given..

_-c \<CONSTANT NAME\> \<CONSTANT VALUE\>_

This flag passes to AMDAT the value of a constant (aka, variable) with name by \<CONSTANT NAME\> with an initial value of \<CONSTANT VALUE\>. AMDAT will replace any text in the input file with the format "${\<CONSTANT NAME\>}" with "\<CONSTANT VALUE\>" during execution. This flag is not required. Constants alternatively may be specified in the input file using the _constant_ command.

_-h or –help_

This flag returns a help menu detailing the AMDAT verstion and date of release as well as the command line format detailed here.

-v

This flag returns the AMDAT version and date of release.

-w \<TIME\>

This flag causes the program to wait for \<TIME\> seconds before beginning.

\> \<LOGFILE\> 2\>&1

This is the standard linux command to redirect screen output to \<LOGFILE\>. This command must be used after all other flags.


<h2>Essential links</h2>

Look in the following links for information on how to:
* Write an [input file](input_file.md)
* Create and use [trajectory_lists](trajectory_lists.md) to select and analyze particle trajectories
* Create and use [trajectory_bin_lists](trajectory_bin_list.md) to efficiently perform spatially-resolved analysis of trajectories.
* Create and use [multibodies](multibodies.md) to create and analyze entities comprised of multiple particles
* Create and use [neighbor_lists](neighborlist.md) to analyze particle adjacencies
* Create and use [value_lists](value_list.md) to analyze additional per-particle quantities

Developer documentation is currently incomplete, but a stub can be found [here](developer_documentation.md).
