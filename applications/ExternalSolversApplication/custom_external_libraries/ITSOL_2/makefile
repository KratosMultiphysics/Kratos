# this makefile is for LINUX machines only 
OBJS = $(addprefix OBJ/, fgmr.o iluk.o ilut.o arms2.o ilutpC.o ilutc.o \
	vbiluk.o vbilut.o auxill.o PQ.o piluNEW.o indsetC.o sets.o \
	MatOps.o tools.o systimer.o misc.o setblks.o svdInvC.o)
AR = ar

#
FC      =  gfortran
FCFLAGS =  -c -g -Wall -I./INC
CC      =  gcc
CCFLAGS =  -c -g -DLINUX -Wall -O3 -I./INC
LIB     = LIB/libitsol.a
#

# clear list of default suffixes, and declare default suffixes
.SUFFIXES: .f .c .o

# default rule to make .o files from .f files
OBJ/%.o : SRC/%.f
	$(FC) $(FCFLAGS) -o $@ -c $<

OBJ/%.o : SRC/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<

#
lib $(LIB): $(OBJS)
	$(AR) -rc $(LIB)  $(OBJS)
#
clean :
	rm -f ${OBJS} *~ core *.cache ${LIB} SRC/*~ SRC/*.cache

cleanall :
	rm -f ${OBJS} *~ core *.cache ${LIB} SRC/*~ SRC/*.cache \
        TESTS/*.o TESTS/*.ex 
