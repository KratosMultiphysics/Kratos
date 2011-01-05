include ../make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../SRC
LIBS	= $(SUPERLULIB) $(BLASLIB) -lm

F77EXM	= f77_main.o hbcode1.o c_fortran_dgssv.o

all:	f77exm

f77exm: $(F77EXM) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F77EXM) $(LIBS) -o $@

c_fortran_zgssv.o: c_fortran_zgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o f77exm
