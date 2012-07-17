# this makefile is for LINUX machines only 
#
LINKS = -L../LIB -litsol \
        -L/project/scicom/scicom00/SOFT/lib/gfort64/ -llapack \
	-L/project/scicom/scicom00/SOFT/lib/gfort64/ -lblas 
	 

FC      =  gfortran
FCFLAGS =  -c -g -Wall -I../INC
CC      =  gcc
CCFLAGS =  -c -g -DLINUX -Wall -O3 -I../INC
LD      =  gfortran
LDFLAGS = 
#
# clear list of default suffixes, and declare default suffixes
.SUFFIXES: .f .c .o
# default rule to make .o files from .f files
.f.o  : ;       $(FC) $(FCFLAGS) $*.f -o $*.o
.c.o  : ;       $(CC) $(CCFLAGS) $*.c -o $*.o
#

default: arms.ex

all: arms.ex iluk.ex ilut.ex iluc.ex vbiluk.ex vbilut.ex

arms.ex: mainARMS.o 
	$(LD) $(LDFLAGS) mainARMS.o $(LINKS) -o arms.ex 

iluk.ex: mainILUK.o 
	$(LD) $(LDFLAGS) mainILUK.o  $(LINKS) -o iluk.ex

ilut.ex: mainILUT.o 
	$(LD) $(LDFLAGS) mainILUT.o $(LINKS) -o ilut.ex

iluc.ex: mainILUC.o   
	$(LD) $(LDFLAGS) mainILUC.o $(LINKS) -o iluc.ex

vbiluk.ex: mainVBILUK.o 
	$(LD) $(LDFLAGS) mainVBILUK.o $(LINKS) -o vbiluk.ex 

vbilut.ex: mainVBILUT.o 
	$(LD) $(LDFLAGS) mainVBILUT.o $(LINKS) -o vbilut.ex
#
clean :
	rm -f *.o *.ex *~ core *.cache OUT/*
