include ../make.inc

LUSRC   = ../SRC
HEADER	= -I$(LUSRC) -I$(MATLAB)/extern/include 
#
# For Matlab Version 4 compatibility: comment out -DV5, add -V4 in FFLAGS
FLAGS   = -O -DV5

all:            mexlusolve mexsuperlu

mexlusolve:	mexlusolve.c $(SUPERLULIB)
		${MATLAB}/bin/mex $(HEADER) ${FLAGS} mexlusolve.c \
		$(SUPERLULIB) $(BLASLIB) -lm

mexsuperlu:	mexsuperlu.c $(SUPERLULIB)
		${MATLAB}/bin/mex $(HEADER) ${FLAGS} mexsuperlu.c \
		$(SUPERLULIB) $(BLASLIB) -lm

clean:	
	rm -f *.o *.mex*
