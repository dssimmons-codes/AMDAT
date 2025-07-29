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

> ([amdat.cpp](../src/amdat.cpp))

### 2. class coordinate

> ([coordinate.h](../src/coordinate.h), [coordinate.cpp](../src/coordinate.cpp))

> Class to store a coordinate in three dimensions, with some associated methods and operator definitions.

### 3. class trajectory

> ([trajectory.h](../src/trajectory.h), [trajectory.cpp](../src/trajectory.cpp))

> Stores information about a particle, molecule, or point and its trajectory through time (stores an array of coordinates).

### 4. class atom_trajectory

> ([atom_trajectory.h](../src/atom_trajectory.h), [atom_trajectory.cpp](../src/atom_trajectory.cpp))

> Stores information about an atom (particle) and its trajectory through time (stores an array of coordinates).

### 5. class molecule

> ([molecule.h](../src/molecule.h), [molecule.cpp](../src/molecule.cpp))

> Stores information about a molecule, its substituent atoms (stored in arrays) and in some cases the center-of-mass trajectory of the molecule.

### 6. class system

> ([system.h](../src/system.h), [system.cpp](../src/system.cpp))

> This is the master class that reads trajectory files and stores information about their contents, including: molecules, atoms, time scheme, box size, and so on. Also contains methods to run loops over subsets of atoms and molecules as well as over times for use by analysis objects.

### 7. class analysis

> ([analysis.h](../src/analysis.h), [analysis.cpp](../src/analysis.cpp))

> Parent class for many analysis techniques, with methods for interfacing with the atom, molecule, and time loops provided by the system class.

> Basic public methods for running this:

> `analyze(string)`

> `analyze(Trajectory_List*)`

### 8. class time_correlation_function

> Parent class for several analysis classes that calculate a time correlation function.

> ([time_correlation_function.h](../src/time_correlation_function.h), [time_correlation_function.cpp](../src/time_correlation_function.cpp))

### 9. class van_hove_self

> ([van_hove_self.h](../src/van_hove_self.h), [van_hove_self.cpp](../src/van_hove_self.cpp))

> Analysis class to calculate self van hove.

### 10. class progress

> ([progress.h](../src/progress.h), [progress.cpp](../src/progress.cpp))

> Class to track progress of a calculation.

### 11. class mean_square_displacement

> ([mean_square_displacement.h](../src/mean_square_displacement.h), [mean_square_displacement.cpp](../src/mean_square_displacement.cpp))

> Class to calculate mean square displacement as a function of time.

### 12. class van_hove_distinct

> ([van_hove_distinct.h](../src/van_hove_distinct.h), [van_hove_distinct.cpp](../src/van_hove_distinct.cpp))

> Class to calculate distinct van hove.

### 13. class spacial_decomposition

> ([spacial_decomposition.h](../src/spacial_decomposition.h), [spacial_decomposition.cpp](../src/spacial_decomposition.cpp))

> Class to spatially decompose the system into smaller cells by sorting atoms at each time into these cells.

### 14. class control

> ([control.h](../src/control.h), [control.cpp](../src/control.cpp))

> Master class that reads input file and calls other objects accordingly.

### 15. class wave_vectors

> ([wave_vectors.h](../src/wave_vectors.h), [wave_vectors.cpp](../src/wave_vectors.cpp))

### 16. class spherical_wave_vectors

> ([spherical_wave_vectors.h](../src/spherical_wave_vectors.h), [spherical_wave_vectors.cpp](../src/spherical_wave_vectors.cpp))

### 17. class wave_density

> ([wave_density.h](../src/wave_density.h, [wave_density.cpp](../src/wave_density.cpp))

### 18. class intermediate_scattering_function

> ([intermediate_scattering_function.h](../src/intermediate_scattering_function.h), [intermediate_scattering_function.cpp](../src/intermediate_scattering_function.cpp))

### 19. class correlation_2d

> ([correlation_2d.h](../src/correlation_2d.h), [correlation_2d.cpp](../src/correlation_2d.cpp))

### 20. class incoherent_scattering_function

> ([incoherent_scattering_function.h](../src/incoherent_scattering_function.h), [incoherent_scattering_function.cpp](../src/incoherent_scattering_function.cpp))

### 21. class debyewaller_dist

> ([debyewaller_dist.h](../src/debyewaller_dist.h), [debyewaller_dist.cpp](../src/debyewaller_dist.cpp))

### 22. class stiffness_dist

> ([stiffness_dist.h](../src/stiffness_dist.h), [stiffness_dist.cpp](../src/stiffness_dist.cpp))

### 23. class non_gaussian_parameter

> ([non_gaussian_parameter.h](../src/non_gaussian_parameter.h), [non_gaussian_parameter.cpp](../src/non_gaussian_parameter.cpp))

### 24. class gaussian_comparison

### 25. class particle_list

### 26. class particle_list_exptime

### 27. class fast_particles

### 28. class slow_particles

### 29. class average_particles

### 30. class radial_debye_waller

### 31. class tokenize

### 32. class mean_square_displacement_2d

### 33. class velocity_autocorrelation

### 34. class strings

### 35. class trajectory_list

### 36. class static_trajectory_list

### 37. class exptime_trajectory_list
