<h1>Input files</h1>

This page provides information on how to write an AMDAT input file. The input file is where one specifies the system to be analyzed and the analyses to be performed and is therefore the core of using AMDAT. AMDAT is packaged with example input files that provide the best starting point, along with this page and its links, for how to construct an input file.

<h2>System block</h2>

The input file always begins with a block of data specifying the trajectory file and metadata describing the trajectory. This block of data must be entered without blank lines. It is structured as follows

```
<system type>
<trajectory file type>
<filename(s)>
<time scheme>
<Additional lines of data as required by trajectory type, described in "Trajectory file types" section below>
```

<h3>system types</h3>

AMDAT divides trajectories into two types, with the type denoted a the `<system_type>` keyword: constant volume, denoted by keyword `system` or `system_nv`, and non-constant volume, denoted by keyword `system_np`. All AMDAT functionalities are available with constant volume systems, whereas the following limitations apply for non-constant volume systems.

1. The only trajectory file type presently supported for non-constant volume trajectories is `custom`.
2. In order to obtain unwrapped coordinates in non-constant volume systems, either unwrapped coordinates must be directly provided by the custom trajectory file or image index values must be provided (`i_x`, `i_y`, and `i_z` in the LAMMPS custom dump file format). Otherwise, any analysis techniques relying on unwrapped coordinates will not work correctly.
3. The "distinct van hove",  "structure_factor", and "isf" analysis methods do not presently work for non-constant volume systems (additional constraints apply to some of these even for constant volume systems). This will hopefully be addressed in a future update.

<h3>trajectory file types</h3>

LAMMPS can read multiple trajectory file types, each requiring a different set of additional lines of input metadata as noted in the structure above. Currently recognized file types and their intended usage are as follows. Follow the links in the table below for information on how to read in each file type.

| File type | Usage and comments |
|----------|----------|
| [xyz](xyz.md) | Reads in xyz files produced by LAMMPS. Generally reads in wrapped coordinates and attempts to infer unwrapped coordinates as best it can. Only works with fixed box size (`system_nv`) systems.|
| [xyz_log](xyz_log.md) | Reads in xyz files produced by LAMMPS, as with xyz. Also reads in a LAMMPS log file to obtain box size information. Only works with fixed box size (`system_nv`) systems. |
| [custom](custom.md) | Reads in a LAMMPS custom trajectory file.|
| [custom_byid](custom_byid.md) | Reads in a LAMMPS custom trajectory file, making use of atom IDs in file.|
| [custom_ manual](custom_manual.md) | Reads in a LAMMPS custom trajectory file. Provides more customizability regarding which columns are assocated with which coordinate data, and also enables read-in of additional columns as [value_list](value_list.md) objects for analysis.|
| [xtc](xtc.md) | Reads in xtc binary format produced by GROMACS. Developed on a much earlier version of GROMACS and requires updated testing.|

<h3>time scheme</h3>

AMDAT can handle multiple types of time schema reflecting different algorithms defining the timesteps contained in the trajectory file and how they should be analyzed for dynamical (two-time) properties. All of these schema, and their syntax, are documented in the page on [time schemes](time_scheme.md). A summary of options and their purpose follows.

| Time scheme | Usage and comments |
|----------|----------|
| `snapshot` | Used to analyze the first frame in the trajectory file; logically does not allow analysis of dynamical properties.|
| `linear` | Used for trajectories that are linearly spaced in time; dynamic analyses employ all possible time gaps within the file. |
| `exponential` | Used for trajectories organized into exponential blocks to enable access to many orders of magnitude of dynamics without excessive file size. Analyzes only time gaps starting at the beginning of an exponential block. |

<h2>Analysis block</h2>

The remainder of the file provides a script specifying the analyses to be performed on the trajectory by AMDAT. Analysis tools are built around a set of data objects that store various types of data for analysis. The most important is the `trajectory_list`, which stores a set of particle trajectories for analysis. AMDAT can also define `trajectory_bin_list` objects, which store particle trajectories within a spatially resolved structure to allow for facile spatially resolved analysis. In addition, AMDAT allows definition of `multibodies`, which are data objects that define sets of multiple particles allowing multibody analysis. AMDAT also allows definition of `neighbor_lists` that track the neighbors of a set of central particles. Finally, it can store `value_lists`, which are versatile data objects that store a single value per time per particle in the list. This latter category also allows direct read-in from custom trajectory files, allowing ancillary data columns in the trajectory file to be associated with particles and analyzed accordingly.

<h3>General commands, logic structures, mathematics, and variables</h3>

AMDAT includes a modest set of control structures including loops, if structures and variables, summarized below. 

| Command | function |
|----------|----------|
| `#` | Comment |
| `constant` | Defines a constant with some fixed value. |
| `user_input` | Pauses script execution and allows the user to manually input commands. |
| `print` | Prints to screen. |
| `wait` | Pauses execution. |
| `exit` | Terminates execution. |
| `skip...endskip` | Skips execution of lines between these commands. |
| `limit` | Sets maximum number of frame pairings to use per time gap. |
| `processors` | Does nothing - placeholder command for when multithreading is implemented. |
| `for...end` | Loops over lines between for and end. |
| `break` | Breaks out of a loop. |
| `if...else` | Execute a block of lines conditionally. |
| `evaluate` | Evaluates a mathematical expression and saves the result to a constant. |
| `round` | Rounds a constant to the nearest integer. |
| `floor` | Rounds a constant down to the closest lower integer. |
| `ceil` | Rounds a constant up to the next higher integer. |


##### comment
Any line in the analysis block beginning with the `#` symbol is ignored.

##### constant

`constant <name> <value>`

Sets the constant with name `<name>` to value `<value>`

##### user_input

`user_input`

Pauses script execution and allows the user to manually input commands as though they were part of the input script. This allows the user to run additional analysis, view or change the value of a constant (using `print` or `constant`, respectively), or use any other command available. When prompted, enter the commands line-by-line, pressing "enter" when finished with a line. When completely finished entering commands, type `done` to execute the commands or `cancel` to cancel execution. When the commands are finished executing (or cancelled), the user will be prompted to enter more commands. This will repeat until the user enters the `done` command with no other commands entered. After exiting, execution will resume where it left off.

##### print

`print ARG1 ARG2...`

Prints everything after the command. Useful for printing the value of a constant: `print ${CONSTANT}`

##### wait

`wait <T>`

Pauses execution for `<T>` seconds

##### exit

`exit`

Terminates script execution

##### skip

```
skip

...

endskip
```

Skips input lines until command `endskip`

##### limit

Sets maximum number of frame pairings to use per time gap. Note that this method is 'dumb': it simply cuts off execution of the relevant algorithm at the specified number of data points and does not choose time separations that are statistically independent. We highly recommend not using this command – it should be reserved for situations in which analysis is taking unacceptably long.

##### processors

Currently a placeholder command which does nothing, as AMDAT does not currently support multithreading.

`processors <N>`

Sets the number of processing cores to be used in multithreaded analysis methods to `<N>`, where $0 < N \le \mathrm{MAXTHREADS}$, where `MAXTHREADS` is the number given by the –n flag when executing AMDAT (default 1).


##### for

```
for <constant_name> <initial_value> <final_value> <(optional)step>

...

end
```

Loops over the region between `for` and `end`, changing the value of `<constant_name>` by `<step>` (defaults to 1) from `<initial_value>` to `<final_value>` (exclusive, i.e., it never actually equals `final_value`, identical to `for (int const=<init>; const<final>; const=const+<step>)` in C++).

If `<initial_value>` is larger than `<final_value>` it will automatically decrement `<constant>` by `<step>` until it reaches `<final_value>`. Changing the value of `<constant_name>` will not affect the loop, `<constant_name>` will have its value changed to the next value as if you never changed it.

##### break

`break`

Breaks out of a loop, returning to the code immediately after the `end` command that corresponds to the loop. Only breaks out of one layer of loops. NOTE: old functionality of stopping AMDAT execution is now moved to `exit` command! (V0.410)

##### if...else

```
if <val1> <compare> <val2>`

...

else

...

end
```

Compares `<val1>` and `<val2>` using `<compare>`. `<compare>` can be: `==`, `!=`, `<`, `>`, `<=`, `>=`

If the comparison evaluates to true, the code between the `if` and `else` (if present) or `end` (if no else present) is executed. If it evaluates to false, then the code between `else` and `end` is executed, if present.

```
if ${const1} == 5

print It is five!

else

print It is not five!

end
```

##### evaluate/eval

`evaluate <constant_name> = <num> <operator> <num> <operator> ...`

Evaluates the mathematical expression to the right of the equals sign and saves its value into the constant `<constant\_name>`. Evaluation is done from left to right with parenthesis evaluated first. Currently supported operations are: `+` `-` `*` `/` `^` `%`

`eval C=((${A}^2)+(${B}^2))^(1/2)`

##### round

`round <constant_name>`

Rounds the value of constant `<constant_name>` to the closest integer.

##### floor

`floor <constant_name>`

Rounds the value of constant `<constant_name>` down to the closest integer.

##### ceiling/ceil

`ceiling <constant_name>`

Rounds the value of constant `<constant_name>` up to the closest integer. Can be shortened from `ceiling` to `ceil`.

<h3>Data objects and their use in trajectory analysis</h3>

AMDAT analysis can involve up to 5 types of data objects, as specified in the following table. The pages linked below on each data object type provide lists of commands to create, manipulate, and analyze these data objects. The creation and analysis of these data objects is the core of an AMDAT run.

| Data object | Purpose |
|----------|----------|
| [trajectory\_list](trajectory_lists.md) | Stores, at each time, a list of particles for analysis. |
| [trajectory\_bin\_list](trajectory_bin_list.md) | Stores, at each time, a list of particles, sorted by spatial location in the box, for analysis. |
| [multibody_list](multibodies.md) | Stores a list of multibodies (potentially different at each time), each of which is comprised of multiple particles. |
| [neighbor\_list](neighbor_list.md) | For a selected set of particles, stores a list of their neigboring particles. |
| [value\_list](value_list.md) | Stores a value for each of a list of particles. |

These types of data objects often interact. For example, there are analysis tools that perform an analysis on a `trajectory_list` and generate a `value_list` corresponding to a time-dependent value for each particle in the original list.

In addition to the analysis methods associated with these data objects, AMDAT contains several [miscellatious analysis tools](misc_tools.md) that do not directly operate on any data object but instead further process the results of a prior analysis.

AMDAT includes a partially implemented capability of storing completed analyses, but the details are lost in the mists of time and must be recovered with some investigation. It appears that use of an "s" flag at some point around the input command will trigger this functionality for some analysis tools, allowing use of the later commands `write_analysis` and `delete_analysis`. This will be clarified in a future manual update.
