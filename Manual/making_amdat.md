<h1>Making AMDAT</h1>

First, ensure that the packages in section A below are properly configured. Then, move to `src/voro++-0.4.6/` and run `make` to compile the voronoi package. To make AMDAT, `cd` to the main AMDAT directory and run `make`. Additionally, you must edit lines 9-11 of the makefile so that they point to the locations of the qvectors folders provided with AMDAT (see more on this in section IIIB below).

## A. Required packages

Fftw3 – must be installed

Xdrfile-1.1b – must be installed, provided with AMDAT distribution. To do so, cd to the xrdfile-1.1b directory included in the AMDAT release, run configure, and then run make.

tnt\_126 – directory must be specified in CPLUS\_INCLUDE\_PATH, provided with AMDAT distribution and available from NIST at http://math.nist.gov/tnt/download.html

jama125 – directory must be specified in CPLUS\_INCLUDE\_PATH, included with AMDAT distribution and available from NIST at http://math.nist.gov/tnt/download.html

## B. Makefile Options

Debugging flags can be turned on for use with valgrind by changing comment locations

The flag –fopenmp must be present for both the CFLAGS and CFLAGSmain variables in the makefile in order to build with openmp support (multithreading). Even if you do not plan on using multithreading, this flag is required (probably).

To compile for use on TACC's Ranger SC, set the variable SERVER in _makefile_ to TACC. This removes dependencies from the compilation that are not available on this system and disables several functions. This may also potentially help with compatibility on other large servers. See Analysis Block – Analysis Tools to see which functionality is removed.

Note that the default compilation scheme is to leave the SERVER variable blank (all functionality included).

It is necessary to specify paths to qvector files for inverse space calculations. The location of 3d, 2d, and 1d qvector files must be specified by the makefile variables WAVEVECTORS3D, WAVEVECTORS2D, and WAVEVECTORS1D, respectively. The paths must end in the stem of the files (_/qvector_ by default) and must be enclosed in _\"_ quotation characters. These data folders are supplied with the AMDAT distribution (qvectors.tar.gz) and may be located anywhere in the file system provided that their location is listed here.
