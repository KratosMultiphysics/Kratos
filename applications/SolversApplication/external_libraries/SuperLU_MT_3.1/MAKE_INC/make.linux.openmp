############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:        March 20, 2013  Version 2.1
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _OPENMP

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = libtmglib.a
SUPERLULIB   = libsuperlu_mt$(PLAT).a

BLASDEF   =
BLASLIB = ../lib/libblas$(PLAT).a
MATHLIB   = -lm
MPLIB     =

#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
 
# Definitions used by CPP:
PREDEFS         = -D__OPENMP
#

CC           = gcc -fopenmp
CFLAGS       = -O3 $(PREDEFS) -DPRNTlevel=0 -DDEBUGlevel=0
#
# There are a few files that should NOT be optimized, which appear in the
# end of SRC/Makefile, taking the flag $(NOOPTS)
# NOOPTS should be defined to explicitly turn off any optimization.
NOOPTS	 = -O0

############################################################################
FORTRAN      = gfortran -fopenmp
FFLAGS       = 
LOADER       = gcc -fopenmp
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DNoChange
