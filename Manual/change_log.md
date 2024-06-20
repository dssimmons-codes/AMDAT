<h1>Change Log</h1>

AMDAT 0.500

Enabled custom definition of "multibodies" – sets of multiple particle trajectories within the same molecule on which analyes relevant to multibody correlations and so on and be performed
Coded "gyration radius" multibody_analysis method that calculates the mean gyration radius of a list of multibodies
Enabled custom defition of new trajectories corresponding to the center of mass or centroid of above trajectories
Eliminated all old functionality associated with molecule center of mass calculation and the old method for calculating molecule gyration radius (which worked only for a single molecule centered on the box center). This functionality has been subsumed by multibody methods.
Made extensive updates to storage of trajectory_lists, etc to employ hash tables for more efficient lookup
Eliminated capability for analysis methods to use system loops – now they must use trajectory lists only
Removed a number of deprecated functionalities such as particle lists and updated code to work without them
Disabled multithreading to to concerns over how data from separate processors was being merged
AMDAT 0.430

Added initial multithreading capability to the following analysis methods: msd, isfs, structure_factor
Added commands: processors, floor, ceiling.
Added support for parenthesis and modulus (%) in evaluate command. Removed support for the word 'equals' and requires using '='.
AMDAT version is now displayed in the screen output at the start of execution.
To update an analysis method to use multithreading you must remove all references to global variables that are changed in the analysis method and add "bool isThreadSafe(){return true;};" to the header file of the analysis method. In addition, you must actually cause the multithreading to happen, look at the system.cpp -> displacement_list method for an example. MSD/ISFS multithread in system.cpp, SF multithreads in trajectory_list.cpp
AMDAT 0.421

Added new analysis method mean_displacement, which calculated the mean displacement as a function of time. Its primary purpose is to look for momentum buildup in simulations.
AMDAT 0.420

Corrected time weighting out of exponential blocks – it was previously bugged . It is likely still bugged for frt = 1, but it should work properly for frt = 0.
Changed meaning of number of timesteps for linear time scheme. It is now number of frames minus 1.
Corrected bug in creation of trajectory list bins based on a plane, in the 'above' direction. It was previously not handling box wrapping correctly.
AMDAT 0.410

Nearly complete rewrite of input file handling. No longer creates temporary input file and instead uses a vector of strings with each line as a seperate element. No obvious difference to end user except that there are no temporary files, but allows for most of the new functions this version.
Initial implementation of commands: remove_bin_list, print, wait, constant, round, for, evaluate/eval, exit, if/else, user_input.
Constants are now able to have their value changed during execution. This is primarily done using the commands: constant, evaluate, for, and round.
Changed break command functionality to exit and repurposed break for breaking out of loops.
Added ability to remove bin lists after creation to free up memory. Uncommented several "delete" commands in trajectory_list_bins which has reduced the amount of memory needed to create additional bin lists.
Lines in the input script can be commented out using the "#" symbol (without quotes) in the same way that it is used in BASH. This includes using in the middle/end of lines (proper commenting of AMDAT scripts can now happen)
Comparison of box sizes between initial and current fatime now uses a float comparison that includes a small margin of error for rounding.
AMDAT 0.400

Critical fix - Unwrapping for isfs was not being handled correctly previously and has been fixed. Previous version worked correctly for cubic boxes if the box size size was used to determine wavenumbers. In other systems, typically very small errors were usually introduced by this problem.
Changes to strings
Strings analysis method now only accepts one timegap (multiple timegaps did not work properly before)
Strings analysis method now creates a trajectory list of the particles in strings; an additional argument for the name of the trajectory list has been added at the end.
AMDAT 0.392

Implemented trajectory_list_decay.
AMDAT 0.391

Fixed issues with certain directions in binning as by distance from a plane
AMDAT 0.390

updated binning as a function of distance to include number of bins and to include binning only above or below the plane.
qvector file paths are now defined in the makefile at compile time rather than in the input file at runtime
new qvector files have been freshly generated for 3d, 2d, and 1d
All analysis methods that use qvectors should now also work with 1d options. The list of options is now xyz, xy, xz, yz, x, y, and z.
AMDAT 0.381

Implemented binning as a function of distance from a plane or point.
Added –t flag
AMDAT 0.380

Implemented binning with structure factor, only works for symmetric.
Implemented binning as a function of distance from another trajectory list.
Corrected msd,composition,structure factor and isfs to return a value of 0 when no particles are in the trajectory list.
Added –w flag
Added error message if constants are not defined.
AMDAT 0.372

Fixed an error in new binning code that caused segmentation faults when particles were on the bin boundary.
AMDAT 0.371

New results have revealed that there is still a bug in center-of-mass determination. Center-of-mass analysis should not be trusted.
Fixed an error in new binning code that caused segmentation faults when particles were outside of binned region.
AMDAT 0.370

Modified and fixed analysis with bins. Bin analysis now defaults to using the non-persistence checked version. Adding a 4th argument to the target line of "1" now turns on persistence checking for dynamic analyses. Also made changes to improve memory management so that binning and binned analysis should not encounter serious memory problems.
AMDAT 0.360

Updated command line syntax for running AMDAT
Constants now taken from command line with –c flag. Constants can no longer be defined within the input file.
–i flag now required before input file
–h updated
–l flag no longer works. Use standard linux command to redirect screen file
–v flag updated to include release date
'Constants' parsing improved to prevent an eror that disallowed consecutive constants in some cases. For real this time.
Various manual updates.
AMDAT 0.351

Can no longer create list from a single bin, this is due to improved functionality for binning (higher density of bins now possible)
Changes to how bins are stored:
trajectory lists for bins are calculated on the fly rather than stored internally
For reference: on a systemwith 12GB RAM, no more than 100 bins are currently possible for a system of 4000 particles when calculating composition
Future changes will allow >1000 bins for all analysis techniques, stay tuned.
Fixed error in wave_vectors.cpp that caused incorrect decomposition of system into shells
Fixed elapsed time output for isfs.
'Constants' parsing improved to prevent an eror that disallowed consecutive constants in some cases.
Various manual updates for clarity.
AMDAT 0.350

Added bin_static_analysis and implemented w/ composition; added big 3 to composition
Added density composition to composition analysis method
Eliminated analysis_value_list (merged into value_list) and made extensive changes corresponding to this
Added autocorrelate_value_list
Added crosscorrelate_value_lists
Added vector_autocorrelation_function
Fixes to Strings (still needs a bit more additional verification, but looks like it probably is working right now)
now uses matrix of particle sizes to determine threshold
corrected algorithm error that was causing severe undercounting of strings
corrected algorithm error that was leading to strings of length zero
Incoherent scattering function (isfs) now works with bin trajectory lists
Altered incoherent scattering function to default to not using multiple parallel times for time spacings across blocks, but to instead just use the first time of each block. The user can now include an additional, optional, argument, of 0 or 1, with 1 specifying to use all times across blocks and 0 simply explicitly using the default setting.
Changed doubles in incoherent_scattering_function to floats for reduced memory usage, but retained double math in some of the core algorithm to prevent substantial rounding error.
ngp (non_gaussian_parameter) now works with trajectory_bin_lists and has big 3
AMDAT 0.340

Added clustered_list
Added invert_list
Fixed asymmetric structure factor
Added write_trajectory_list_full() method to write out xyz with dummy values for particles not in the trajectory list
AMDAT 0.330

Corrected output of thresholds from compare_gaussian
String method now writes image file to a file named by taking the filename given by user and appending "_img.xyz"
Created new Structure_Factor class that calculates structure factor with greatly reduced memory use. Modified control and input syntax accordingly.
Made corrections to molecule center of mass calculation. Eliminated valgrind memcheck errors. Some additional warnings regarding ambigious unwrapping are possible now.
AMDAT 0.321

Fixed Bug in system preventing the use of xyz files.
AMDAT 0.320

Changed functionality of structure_factor() in control and wave_density to reduce calculation time of wave densities.
Added copy constructor, default constructor, assignment operator, and destructors to
intermediate_scattering_function
wave_density
wave_vectors
Changes to remove unsupported dependencies for usage on TACC's Ranger SC cluster. Makefile and various methods changed.
Added ability to create a list from an individual bin from a previously created bin list.
Bug fixes
define_constant: multiple uses of a variable in a line
create_bin_list: type region now functions properly (added conditional statement to trajectory_list_bins::bin_boolean_calculation() ).
custom: boolean definitions in system updated.
AMDAT 0.311

Made corrections to xyz_log for reading box size
AMDAT 0.310

Made some corrections to read-in of custom dump files and template files
AMDAT 0.300

Added ability to read non-constant volume files.
Modified multiple analysis methods to work with new non-constant box dimensions.
Corrected problems with read-in of custom trajectory with template file
Improved algorithm for calculating molecule center of mass trajectories. The algorithm now searches for the first time at which center of mass calculation is unambiguous and then uses unwrapped particle displacements to determine the COM position at all other times. As a consequence, instead of requiring COM calculation to be unambiguous at all times, the algorithm only requires it to be unambiguous at one time. In the case that a molecule has no unambiguous time, the algorithm will print a single warning for that molecule COM, a considerable reduction from the large number of warnings frequently generated previously.
AMDAT 0.260

Can now read LAMMPS custom dump files (this feature is still in beta and should be tested further)
Added command "define_constant", which allows user to define a constant in input file which will be replaced any time it is found. Constant can be used with the syntax ${constant_name}.
Can now determine box bound information from a LAMMPS log file rather than from input file, if log file is provided
Corrected an error in the way trajectory type is determined for center of mass trajectories
Added ability to read in particle velocity information from a custom LAMMPs dump file and store it in trajectory objects; also added methods to class trajectory to interact with this velocity data.
AMDAT 0.250

Structure factor updated to work in 2D with plane options.
msd_2d compatible with binning.
Unused variables removed to eliminate warning messages on compilation.
Corrected center-of-mass calculations.
AMDAT 0.240

Binning scheme added with limited functionality, see documentation for capabilities.
New classes: trajectory_list_bins and bin_dynamics_analysis
Significant changes to control including new run_analysis template to be used with any method that uses binning
write_bins_xyz capability included, see documentation
Any analysis method must have the following changes to utilize binning (see msd for example)
Copy constructor
Destructor
Assignment operator
bin_hook method to allow bin_dyanmics_analysis's displacement kernel to call an analysis object's displacement kernel.
Must use run_analysis template in control
Assignment operator, copy constructor, and destructor added for trajectory_list, boolean_list, and mean_square_displacement.
bin_hook method added to msd.
hexatic_order_parameter replaced with n_fold_order_parameter
Added composition class, which allows computation of the time dependent and average compositions
added analysis_value_list, which stores an array of value_lists
update to control to be able to call stored analysis_value_lists
added skip command to control
AMDAT 0.230

Corrections and assignment operators for trajectory_list and Boolean_list
Made some OOP and const-correctness improvements to the core code. This may give some minor performance improvements and should help in error prevention.
Use of passing by const reference where appropriate to improve efficiency
Use of const class methods where appropriate
Returning constant pointers from class methods when appropriate to prevent unintended/unauthorized modification of internal data members
In trajectory list public 'set' command taking a Boolean_List pointer argument, switched from copying pointer over to doing a deep copy of the Boolean_List objects referred to by the pointer – this should prevent unintentional modification of the included trajectory list by external code
Began adding copy constructors and assignment operators to all classes where the defaults are inadequate due to use of pointers.
AMDAT 0.220

Eliminated spurious warning associated with molecule center of mass trajectory calculation
Updated value_list.h class with new write_pdb method that allows the class to write a pdb file that VMD will read and color according to the values stored in value_list.h
Corrected ngp calculation via trajectory list method – it was incorrect in prior versions
Slightly modified the way in which 'ambiguous unwrapping' warnings are generated in a way that should catch more ambiguous unwrappings
Added hexatic_order_parameter class, which allows computation of average and time dependent hexatic order parameter and also enables creation of hexatic order parameter maps. This class works only for 2-D and pseudo-2-D systems, and presently is not implemented for center of mass trajectories, only for single atom trajectories. It also presently excludes particles within a cutoff distance of the edge of the box.
Added displacement_map class, which generates a map (or maps) of particle displacements at user-specified times.
AMDAT 0.210

Added option for additional command line parameter specifying file to redirect log output to.
Corrected errors related to unwrapping and molecule COM calculation warning generation
AMDAT 0.200

Created new timetype snapshot, which loads a single frame (the first) from a trajectory
Internally reworked system loops to pass all trajectory objects by reference rather than passing their species, molecule, type, and atom indices
Corrected calculation of molecule center of mass trajectories and made the calculation automatic during system construction
Most analysis techniques should now work with molecule center of mass trajectories
Created new trajectory set (analysis target) commands moleculecom_all and moleculecom_species, which perform analysis on molecule centers of mass of all molecules and molecules of a given species, respectively.
Created new analysis.write(string) command which serves as a prototype for data writing commands for all analysis classes. All data writing by analysis classes must now be handled by this command.
Fixed structure factor and intermediate scattering function calculation from lists
Created new command write_list_trajectory, which writes an xyz trajectory file from a trajectory list.
AMDAT 0.130

Created class Boolean_List, which stores a 0 or 1 for every trajectory in the system. Objects of this class are used as an inclusion list by Trajectory_List class to allow for Boolean comparisons between trajectory lists.
Created templated class Value_List, which stores a value of variable type for every trajectory in the system, and also stores a Boolean_List, signifying which trajectories are included
Enhanced Trajectory_List classes to support Boolean operations between Trajectory_List objects and between times of a given object
Expanded rgtensor_stats command and class to output additional data about rgtensor distributions.
