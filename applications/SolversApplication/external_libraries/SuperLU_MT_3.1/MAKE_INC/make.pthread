############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   August 15, 1997
#
#  Modified:	    September 1, 1999 version 1.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _PTHREAD

#
#  The name of the libraries to be created/linked to
#
TMGLIB       	= libtmglib$(PLAT).a
SUPERLULIB   	= libsuperlu_mt$(PLAT).a
BLASDEF	        = -DUSE_VENDOR_BLAS
BLASLIB		= -L/usr/lib -lf77blas -latlas
#BLASLIB = ../lib/libblas$(PLAT).a
MATHLIB		= -lm 
MPLIB		= -lpthread

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH          	= ar
ARCHFLAGS     	= cr
RANLIB       	= ranlib

############################################################################
#
# Definitions used by CPP:
PREDEFS      	= -D__PTHREAD  ## -DDEBUGlevel=0 -DPRNTlevel=0
#
CC           	= gcc
CFLAGS       	= $(PREDEFS) -g #-O3
# Uncomment the following to use 64-bit integer
# CFLAGS          += -D_LONGINT

#
# There are a few files that should NOT be optimized, which appear in the
# end of SRC/Makefile, taking the flag $(NOOPTS)
# NOOPTS should be defined to explicitly turn off any optimization.
NOOPTS	     	= -O0

BLASOPTS	=

############################################################################

FORTRAN		= gfortran
FFLAGS		= -O 
# Uncomment the following to use 64-bit integer
# FFLAGS         += -fdefault-integer-8

############################################################################

LOADER       	= gcc
LOADOPTS     	= -O #-pthread

############################################################################
#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        	= -DAdd_
