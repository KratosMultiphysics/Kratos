############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   August 15, 1997
#
#  Modified:   	    September 1, 1999 version 1.0
#                   March 20, 2013  version 2.1
#
############################################################################

include make.inc

all: install lib testing

lib: superlulib tmglib

clean: cleanlib cleantesting

install:
	( cd INSTALL; $(MAKE) )
#	( cd INSTALL; cp lsame.c ../SRC/; \
#	  cp dlamch.c ../SRC/; cp slamch.c ../SRC/ )

blaslib:
	( cd CBLAS; $(MAKE) )

superlulib:
	( cd SRC; $(MAKE) )

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

testing:
	( cd TESTING ; $(MAKE) )

cleanlib:
	( cd SRC; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )
	( cd lib; rm -f *.a )
	( rm -f *.a )

cleantesting:
	( cd INSTALL; $(MAKE) clean )
	( cd TESTING; $(MAKE) clean )
	( cd EXAMPLE; $(MAKE) clean )
