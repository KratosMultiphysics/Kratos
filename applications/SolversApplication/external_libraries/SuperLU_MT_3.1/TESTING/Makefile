include ../make.inc

#######################################################################
#  This makefile creates the test programs for the linear equation
#  routines in SuperLU.  The test files are grouped as follows:
#
#       ALINTST -- Auxiliary test routines
#       SLINTST -- Single precision real test routines
#       DLINTST -- Double precision real test routines
#       CLINTST -- Double precision complex test routines
#       ZLINTST -- Double precision complex test routines
#
#  Test programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four test programs.
#  The executable files are called
#       pstest
#       pdtest
#       pctest
#       pztest
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

HEADER  = ../SRC

ALINTST = sp_ienv.o

SLINTST = psdrive.o sp_sconvert.o psgst01.o psgst02.o psgst04.o psgst07.o \
	  psgssv.o psgssvx.o

DLINTST = pddrive.o sp_dconvert.o pdgst01.o pdgst02.o pdgst04.o pdgst07.o \
	  pdgssv.o pdgssvx.o

CLINTST = pcdrive.o sp_cconvert.o pcgst01.o pcgst02.o pcgst04.o pcgst07.o \
	  pcgssv.o pcgssvx.o

ZLINTST = pzdrive.o sp_zconvert.o pzgst01.o pzgst02.o pzgst04.o pzgst07.o \
	  pzgssv.o pzgssvx.o

all: single double complex complex16

pstest: $(SLINTST) $(ALINTST) ../lib/$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(SLINTST) $(ALINTST) \
	$(TMGLIB) ../lib/$(SUPERLULIB) $(BLASLIB) -lm $(MPLIB) -o $@

pdtest: $(DLINTST) $(ALINTST) ../lib/$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(DLINTST) $(ALINTST) \
	$(TMGLIB) ../lib/$(SUPERLULIB) $(BLASLIB) -lm $(MPLIB) -o $@

pctest: $(CLINTST) $(ALINTST) ../lib/$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(CLINTST) $(ALINTST) \
	$(TMGLIB) ../lib/$(SUPERLULIB) $(BLASLIB) -lm $(MPLIB) -o $@

pztest: $(ZLINTST) $(ALINTST) ../lib/$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(ZLINTST) $(ALINTST) \
	$(TMGLIB) ../lib/$(SUPERLULIB) $(BLASLIB) -lm $(MPLIB) -o $@

single: pstest pstest.out

double: pdtest pdtest.out

complex: pctest pctest.out

complex16: pztest pztest.out

pstest.out: pstest pstest.csh
	@echo Testing SINGLE PRECISION linear equation routines 
	csh pstest.csh

pdtest.out: pdtest pdtest.csh
	@echo Testing DOUBLE PRECISION linear equation routines
	csh pdtest.csh

pctest.out: pctest pctest.csh
	@echo Testing SINGLE COMPLEX PRECISION linear equation routines
	csh pctest.csh

pztest.out: pztest pztest.csh
	@echo Testing DOUBLE COMPLEX PRECISION linear equation routines
	csh pztest.csh

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

clean:	
	rm -f *.o p*test *.out

