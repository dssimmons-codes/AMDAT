<h1>Input files</h1>

This page provides information on how to write an AMDAT input file. The input file is where one specifies the system to be analyzed and the analyses to be performed and is therefore the core of using AMDAT. AMDAT is packaged with example input files that provide the best starting point, along with this page and its links, for how to construct an input file.

<h2>System block</h2>

The input file always begins with a block of data specifying the trajectory file and metadata describing the trajectory. This block of data must be entered without blank lines. It is structured as follows

_\<system type\>_

_\<trajectory file type\>_

_\<filename(s)\>_

_\<time scheme\>_

_\<Additional lines of data as required by trajectory type, described in "Trajectory file types" section\>_

<h3>system types</h3>

AMDAT divides trajectories into two types, with the type denoted a the _\<system\_type\>_ keyword: constant volume, denoted by keyword _system_ or _system\_nv_, and non-constant volume, denoted by keyword _system\_np_. All AMDAT functionalities are available with constant volume systems, whereas the following limitations apply for non-constant volume systems.

1. The only trajectory file type presently supported for non-constant volume trajectories is _custom_.
2. In order to obtain unwrapped coordinates in non-constant volume systems, either unwrapped coordinates must be directly provided by the custom trajectory file or image index values must be provided (i\_x, i\_y, and i\_z in the LAMMPS custom dump file format). Otherwise, any analysis techniques relying on unwrapped coordinates will not work correctly
3. The "distinct van hove" analysis method does not presently work for non-constant volume systems

<h3>trajectory file types</h3>

LAMMPS can read multiple trajectory file types, each requiring a different set of addition lines of input metadata as noted in the structure above. Currently recognized file types and their intended usage are as follows. Follow the links in the table below for information on how to read in each file type.

| File type | Usage and comments |
|----------|----------|
| [xyz](read_xyz.md) | Reads in xyz files produced by LAMMPS. Generally reads in wrapped coordinates and attempts to infer unwrapped coordinates as best it can. Only works with fixed box size (system\_nv) systems. |
| [xyz_log](read_xyz_log.md | Reads in xyz files produced by LAMMPS, as with xyz. Also reads in a LAMMPS log file to obtain box size information. Only works with fixed box size (system\_nv) systems. |
| [custom](read_custom.md) | Reads in a LAMMPS custom trajectory file.|
| [custom_ manual](read_custom_manual.md) | Reads in a LAMMPS custom trajectory file. Provides more customizability regarding which columns are assocated with which coordinate data, and also enables read-in of additional columns as [value_list](value_list.md) objects for analysis.|
| [xtc](read_xtc.md) | Reads in xtc binary format produced by GROMACS. Developed on a much earlier version of GROMACS and requires updated testing. |
