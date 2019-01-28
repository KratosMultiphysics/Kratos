include ./make.inc

ARCH:=none

LIB:=none

help:
	@echo "Usage: make ARCH=<arch> LIB={feast,pfeast} {all, clean} "
	@echo "--- <arch> is your architecture type and the name of the directory in ../lib"
	@echo "     where the FEAST library will be created !!!! "
	@echo "--  LIB=feast  compile and create or clean FEAST     library > lfeast lfeast_dense lfeast_banded lfeast_sparse"
	@echo "--  LIB=pfeast compile and create or clean FEAST-MPI library > lpfeast lpfeast_dense lpfeast_banded lpfeast_sparse"
	@echo

ifeq ($(OPTION),2)
FUNCTIONS_WRAPPER:=c
else
FUNCTIONS_WRAPPER:=f90
endif

ifeq ($(LIB),none)
LIB:=feast
endif

ifeq ($(ARCH),none)
ARCH:=x64
endif

#==============================================
# Name of Libraries 
#==============================================

################ Sequential or Parallel
ifeq ($(LIB),feast)
FF=$(F90)
FFFLAGS=$(F90FLAGS)
NAME_dense = libfeast_dense.a
NAME_banded = libfeast_banded.a
NAME_sparse = libfeast_sparse.a
NAME = libfeast.a
else
FF=$(PF90)
FFFLAGS=$(PF90FLAGS)
NAME_dense = libpfeast_dense.a
NAME_banded = libpfeast_banded.a
NAME_sparse = libpfeast_sparse.a
NAME = libpfeast.a
endif

#sources Fortran90
CODE90_functions_wrapper = $(PWD)/kernel/f90_functions_wrapper.f90
CODE90 = $(PWD)/kernel/feast_tools.f90 $(PWD)/kernel/dzfeast.f90 $(PWD)/kernel/scfeast.f90 $(PWD)/kernel/feast_aux.f90 
CODE90_banded = $(PWD)/banded/dzfeast_banded.f90 $(PWD)/banded/scfeast_banded.f90 
CODE90_dense =  $(PWD)/dense/dzfeast_dense.f90 $(PWD)/dense/scfeast_dense.f90 
CODE90_sparse =  $(PWD)/sparse/dzfeast_sparse.f90 $(PWD)/sparse/scfeast_sparse.f90 

#sources C
CODEC_functions_wrapper = $(PWD)/kernel/c_functions_wrapper.c

ifeq ($(FUNCTIONS_WRAPPER),c)
	OBJECTS =  $(CODEC_functions_wrapper:.c=.o)  $(CODE90:.f90=.o) 
	OBJECTS_dense =  $(CODE90_dense:.f90=.o)  
	OBJECTS_banded = $(CODE90_banded:.f90=.o) 
	OBJECTS_sparse = $(CODE90_sparse:.f90=.o)   
else
	OBJECTS =  $(CODE90_functions_wrapper:.f90=.o)  $(CODE90:.f90=.o)  
	OBJECTS_dense =  $(CODE90_dense:.f90=.o) 
	OBJECTS_banded =  $(CODE90_banded:.f90=.o)
	OBJECTS_sparse =  $(CODE90_sparse:.f90=.o)
endif



####################################


# .SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: clean all

#.f90.o: $(CODE90)
%.o: %.f90
	$(FF) $(FFFLAGS) -c $< -o $@

#.f.o: $(CODEC)
%.o: %.f
	$(CC) $(CFLAGS)  -c $< -o $@


#all: $(NAME)

#$(NAME): $(OBJECTS) 

all:	$(OBJECTS) $(OBJECTS_dense) $(OBJECTS_banded) $(OBJECTS_sparse)  
	ar r $(NAME)  $(OBJECTS)
	ar r $(NAME_dense) $(OBJECTS_dense)
	ar r $(NAME_banded) $(OBJECTS_banded)
	ar r $(NAME_sparse) $(OBJECTS_sparse)
	@if test -d ../lib/$(ARCH); then \
	echo "directory $(ARCH) already exists"; \
	else echo "create directory $(ARCH)" && mkdir ../lib/$(ARCH); fi;
	mv *.a ../lib/$(ARCH)
	rm -f $(OBJECTS) $(OBJECTS_dense) $(OBJECTS_banded) $(OBJECTS_sparse) 


clean: 
	rm ../lib/$(ARCH)/$(NAME)
	rm ../lib/$(ARCH)/$(NAME_dense)
	rm ../lib/$(ARCH)/$(NAME_banded)
	rm ../lib/$(ARCH)/$(NAME_sparse)
	rm -r ../lib/$(ARCH)
touch:
	@touch $(CODE90) $(CODE90_dense) $(CODE90_banded) $(CODE90_sparse)













