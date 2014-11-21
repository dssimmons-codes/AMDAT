## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

bin_PROGRAMS = trr2xtc

lib_LTLIBRARIES = libxdrfile@LIBSUFFIX@.la

LDADD = $(lib_LTLIBRARIES)

libxdrfile@LIBSUFFIX@_la_LIBADD         = 
libxdrfile@LIBSUFFIX@_la_DEPENDENCIES   = 
libxdrfile@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@

# Note: you don't have to list sources for "prog" if it is the single file prog.c

libxdrfile@LIBSUFFIX@_la_SOURCES = \
			xdrfile.c	\
			xdrfile_trr.c 	\
			xdrfile_xtc.c	

EXTRA_PROGRAMS = xdrfile_c_test

test:	xdrfile_c_test
	./xdrfile_c_test

CLEANFILES   = 	*~ \\\#*










