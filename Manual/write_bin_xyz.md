<h1>write_bin_xyz</h1>

writes xyz file for each bin, prepending bin index to filename. Note that this may produce a large number of files. Moreover, if the number of particles in each bin changes with time (as is typically the case), this will write files containing a different number of particles at each time frame, which may not be tolerated by many visualization software packages. Note that, unlike most analysis commands, this command involves two lines. Syntax is as follows.

Line 1: _write\_bin\_xyz \<name of bin\_list\> \<stem of xyz filenames to write\>_
Line 2: _\<write type\> \<arguments\>_

_\<write type\>_ is either "all" or "single. If it is "all", then an xyz file will be written for every bin in the bin\_list. In this case, no arguments are required or allowed on the second line. If it is "single", then an xyz file will be written for only a single bin. In this case the bin must be specified with three arguments, providing the x, y, and z indices of the targeted bin (in that order).
