# makefile for TetGen
#
# Type "make" to compile TetGen into an executable program (tetgen).
# Type "make tetlib" to compile TetGen into a library (libtet.a).
# Type "make distclean" to delete all object (*.o) files.

# CC should be set to the name of your favorite C++ compiler.

CC = g++

# OPT is the level of optimiztion, default is -O. One should try -O2, -O3
#   ... to find the best optimization level.

OPT = -O3 -fPIC

# CFLAGS is a list of switches to compile TetGen.
#
# By default, TetGen uses double precision floating point numbers.  If you
#   prefer single precision, use the -DSINGLE switch. 
#
# The source code of TetGen includes a lot of assertions, which are mainly
#   used for catching bugs at that places.  These assertions somewhat slow
#   down the speed of TetGen.  They can be skipped by define the -DNDEBUG
#   switch.

#CFLAGS = -Wall -DSELF_CHECK -DNDEBUG
CFLAGS = -Wall -DNDEBUG

#CFLAGS = -funroll-all-loops -fomit-frame-pointer\
#         -Wall -Wabi -Wctor-dtor-privacy \
#         -Woverloaded-virtual -Wno-pmf-conversions -Wsign-promo \
#         -Wsynth  -Wchar-subscripts -Wconversion -Wsign-compare \
#         -Wcomment  -Wimplicit -Wmissing-braces -Wparentheses \
#         -Wreturn-type -Wswitch -Wswitch-default \
#         -Wswitch-enum -Wtrigraphs -W

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm

# The action starts here.

tetgen:	tetgen.cxx predicates.o
	$(CC) $(CFLAGS) $(OPT) -o tetgen tetgen.cxx predicates.o -lm

tetlib: tetgen.cxx predicates.o
	$(CC) $(CFLAGS) $(OPT) -DTETLIBRARY -c tetgen.cxx
	ar r libtet.a tetgen.o predicates.o

#predicates.o: predicates.cxx
#	$(CC) $(CFLAGS) -g -c predicates.cxx
predicates.o: predicates.cxx
	$(CC) $(CFLAGS) -fPIC -DLINUX -c predicates.cxx

distclean:
	$(RM) -rf $(SRC)*.o $(SRC)*.a




