<h1>write_single_particle</h1>

Writes out an xyz trajectory file based on a single selected particle in the trajectory. In this xyz file, the coordinates of the selected particles at different times are written out as though they are the positions of multiple particles at one time, with earlier times corresponding to lower-indexed particles. This can facilitate visualization of the spatial trajectory of one particle as a single snapshot. Syntax is as follows.

_write\_single\_particle \<output filename\> \<index of particle to write out\>

The scheme for particle selection here is rather simplistic. It does not employ molecular structure in any way. It is simply the index of the particle as it is stored in AMDAT, beginning with index 0. This typically corresponds to the index of the particle in the trajectory file that was initially read in by AMDAT. Put another way, if one examines the initial frame in the original trajectory file and counts down particles from the top starting with 0, the index here will be the index of the particle in that count.
