<h1>Developer Documentation</h1>

This developer documentation is somewhat out of date, but it may still provide some help. This will be expanded at some point in the future.

Broadly speaking, there are two parts of the AMDAT code: 1) the data structure and 2) analysis classes.

The purpose of the data structure is to store trajectories. A `trajectory` is a sequence of coordinates in time, with some additional metadata. Common trajectories include trajectories of atoms (or particles), and potentially also molecular centers of mass. The core data structure includes the classes `System`, `Trajectory`, `Molecule`, `Atom_Trajectory`, and `Coordinate`. `Molecule` and `Atom_Trajectory` are both daughter classes of `Trajectory`, and any `Trajectory` contains an array of `Coordinates`.

The analysis classes comprise a set of classes which operate on this data structure to yield output data. There are generally two ways by which an analysis class may interact with the data structure:

1. It may use preexisting loops in the `System` class to loop over rational subsets of trajectories.
2. It may use an intermediary set of classes, `Trajectory_List` and its daughters, to loop over a subset of trajectories previously defined by the user.

Presently, all analysis classes should preferably be made compatible with both methods; however, method 2) should be given priority as method 1) may be phased out at a later time.

<!-- Classes shown in blue are those that inherit the functionality of parent class analysis. Classes shown in green are part of the -->

### 1. main

> (amdat.cpp)

### 2. class coordinate

> (coordinate.h, coordinate.cpp)

> Class to store a coordinate in three dimensions, with some associated methods and operator definitions.

### 3. class trajectory

> (trajectory.h, trajectory.cpp)

> Stores information about a particle, molecule, or point and its trajectory through time (stores an array of coordinates)

### 4. class atom\_trajectory

> (atom_trajectory.h, atom_trajectory.cpp)

> Stores information about an atom (particle) and its trajectory through time (stores an array of coordinates)

### 5. class molecule

> (molecule.h, molecule.cpp)

> Stores information about a molecule, its substituent atoms (stored in arrays) and in some cases the center-of-mass trajectory of the molecule

### 6. class system

> (system.h, system.cpp)

> This is the master class that reads trajectory files and stores information about their contents, including: molecules, atoms, time scheme, box size, and so on. Also contains methods to run loops over subsets of atoms and molecules as well as over times for use by analysis objects.

### 7. class analysis

> (analysis.h, analysis.cpp)

> Parent class for many analysis techniques, with methods for interfacing with the atom, molecule, and time loops provided by the system class.

> Basic public methods for running this

> `analyze(string)`

> `analyze(Trajectory_List*)`

### 8. class time\_correlation\_function

Parent class for several analysis classes that calculate a time correlation function.

(time\_correlation\_function.h, time\_correlation\_function.cpp)

class van\_hove\_self

Analysis class to calculate self van hove.

(van\_hove\_self.h, van\_hove\_self.cpp)

class progress

(progress.h, progress.cpp)

Class to track progress of a calculation.

### 9. class mean\_square\_displacement

(mean\_square\_displacement.h, mean\_square\_displacement.cpp)

Class to calculate mean square displacement as a function of time.

### 10. class van\_hove\_distinct

(van\_hove\_distinct.h, van\_hove\_distinct.cpp)

Class to calculate distinct van hove.

### 11. class spacial\_decomposition

(spacial\_decomposition.h, spacial\_decomposition.cpp)

Class to spatially decompose the system into smaller cells by sorting atoms at each time into these cells.

### 12. class control

(control.h, control.cpp)

Master class that reads input file and calls other objects accordingly.

### 13. class wave\_vectors

(wave\_vectors.h, wave\_vectors.cpp)

### 14. class spherical\_wave\_vectors

(spherical\_wave\_vectors.h, spherical\_wave\_vectors.cpp)

### 15. class wave\_density

(wave\_density.o, wave\_density.cpp)

### 16. class intermediate\_scattering\_function

(intermediate\_scattering\_function.h, intermediate\_scattering\_function.cpp)

### 17. class correlation\_2d

(correlation\_2d.h, correlation\_2d.cpp)

### 18. class incoherent\_scattering\_function

(incoherent\_scattering\_function.h, incoherent\_scattering\_function.cpp)

### 19. class debyewaller\_dist

(debyewaller\_dist.h, debyewaller\_dist.cpp)

### 20. class stiffness\_dist

(stiffness\_dist.h, stiffness\_dist.cpp)

### 21. class non\_gaussian\_parameter

(non\_gaussian\_parameter.h, non\_gaussian\_parameter.cpp)

### 22. class gaussian\_comparison

### 23. class particle\_list

### 24. class particle\_list\_exptime

### 25. class fast\_particles

### 26. class slow\_particles

### 27. class average\_particles

### 28. class radial\_debye\_waller

### 29. class tokenize

### 30. class mean\_square\_displacement\_2d

### 31. class velocity\_autocorrelation

### 32. class strings

### 33. class trajectory\_list

### 34. class static\_trajectory\_list

### 35. class exptime\_trajectory\_list
