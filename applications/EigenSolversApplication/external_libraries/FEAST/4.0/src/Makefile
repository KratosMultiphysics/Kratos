ARCH?=x64

F90?=ifort

MKL?=yes

CC?=gcc

MPI?=impi


help:
	@echo "---------------------------------------------------------------------------------------------------------------"
	@echo "Usage: make ARCH=<arch> F90=<f90> MPI=<mpi> MKL=<mkl> {feast, pfeast}"
	@echo ""
	@echo "Library: 'make feast' will create the FEAST library 'libfeast.a' in ../lib/<arch>"
	@echo "         'make pfeast' will create the PFEAST library 'libpfeast.a' in ../lib/<arch>"
	@echo ""
	@echo "Default: make {feast, pfeast}"
	@echo "         ==is equivalent to== "
	@echo "         make ARCH=x64 F90=ifort MPI=impi MKL=yes {feast, pfeast}"
	@echo ""
	@echo "Other possible options:"
	@echo "--- <arch> is the name of the directory in ../lib where libraries will be located. Default is x64"
	@echo "--- <f90>  is your own Fortran90 compiler (possible choices: ifort, gfortran, pgf90). Default is ifort"
	@echo "--- <mpi>  (mandatory for pfeast only). it is your MPI library (possible choices: impi, mpich, openmpi)."
	@echo "           Default is impi (intel MPI)"
	@echo "--- <mkl>  enables Intel-MKL math library instructions (possible choices: yes,no). Defaults is yes"
	@echo "          if MKL=yes, Your own code must be linked with Intel MKL"
	@echo "          if MKL=no,  Your own code can be linked with any BLAS/LAPACK"
	@echo "                      (i) MKL-PARDISO and cluster-PARDISO will not be available for the FEAST sparse interfaces"
	@echo "                     (ii) in-built sparse mat-vec routines (for IFEAST sparse interfaces) will be slower"
	@echo "---------------------------------------------------------------------------------------------------------------"

############################################################
######### Fortran Compiler set-up
############################################################


# intel fortran
ifeq ($(F90),ifort)
F90FLAGS= -O3 -qopenmp -cpp -fPIC
ifeq ($(MKL),yes)
F90FLAGS += -DMKL
endif
#F90FLAGS += -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces  #-heap-arrays 1
endif

# gnu fortran
ifeq ($(F90),gfortran)
F90FLAGS= -O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none -cpp #-fPIC
ifeq ($(MKL),yes)
F90FLAGS += -DMKL
endif
endif

# portland group fortran compiler
ifeq ($(F90),pgf90)
F90FLAGS= -O3 -mp -cpp -fPIC
ifeq ($(MKL),yes)
F90FLAGS += -DMKL
endif
endif


##############################################################
######### include MPI directives
##############################################################


### intel mpi
ifeq ($(MPI),impi)
PF90=mpiifort -f90=$(F90)
PF90FLAGS= $(F90FLAGS) -DMPI
endif


### mpich
ifeq ($(MPI),mpich)
PF90= mpif90.mpich -fc=$(F90)
PF90FLAGS= $(F90FLAGS) -DMPI
endif

### openmpi ...requires shell environment variable "OMPI_FC=$(F90)"
ifeq ($(MPI),openmpi)
export OMPI_FC=$(F90) # for BASH shell
PF90= mpif90.openmpi 
PF90FLAGS= $(F90FLAGS) -DMPI
endif



##############################################################

#==============================================
# Name of Libraries 
#==============================================

NAME = libfeast.a
PNAME = libpfeast.a


###################################
#### SOURCE CODE FILES
###################################


CODE90 = $(PWD)/kernel/dzfeast.f90 $(PWD)/kernel/feast_aux.f90  $(PWD)/kernel/feast_tools.f90  $(PWD)/kernel/libnum.f90 #$(PWD)/kernel/scfeast.f90 
CODE90_dense = $(PWD)/dense/dzfeast_dense.f90 $(PWD)/dense/dzfeast_pev_dense.f90
CODE90_banded = $(PWD)/banded/dzfeast_banded.f90
CODE90_sparse = $(PWD)/sparse/dzfeast_sparse.f90 $(PWD)/sparse/dzifeast_sparse.f90 $(PWD)/sparse/dzlsprim.f90 $(PWD)/sparse/sclsprim.f90  
CODE90_sparse += $(PWD)/sparse/dzfeast_pev_sparse.f90 $(PWD)/sparse/dzifeast_pev_sparse.f90
CODE90_psparse =  $(PWD)/sparse/pdzfeast_sparse.f90 $(PWD)/sparse/pdzfeast_pev_sparse.f90 $(PWD)/sparse/pdzifeast_sparse.f90 $(PWD)/sparse/pdzifeast_pev_sparse.f90

################################################
############ OBJECTS############################
################################################

OBJECTS =  $(CODE90:.f90=.o)
OBJECTS_dense =  $(CODE90_dense:.f90=.o)
OBJECTS_banded =  $(CODE90_banded:.f90=.o)
OBJECTS_sparse =  $(CODE90_sparse:.f90=.o)

OBJECTS_p =  $(CODE90:.f90=p.o)
OBJECTS_pdense =  $(CODE90_dense:.f90=p.o)
OBJECTS_pbanded =  $(CODE90_banded:.f90=p.o)
OBJECTS_psparse =  $(CODE90_sparse:.f90=p.o) $(CODE90_psparse:.f90=p.o)



####################################

# .SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: clean all

#.f90.o: $(CODE90)
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%p.o: %.f90
	$(PF90) $(PF90FLAGS) -c $< -o $@

#.c.o: $(CODEC)
%.o: %.c
	$(CC) $(CFLAGS)  -c $< -o $@


###################################################
########### COMMANDS
###################################################

feast:	$(OBJECTS)  $(OBJECTS_sparse) $(OBJECTS_banded) $(OBJECTS_dense)
	ar cr $(NAME)  $(OBJECTS_sparse) $(OBJECTS_banded) $(OBJECTS_dense) $(OBJECTS) 
#ar r $(NAME_dense) $(OBJECTS_dense)
#ar r $(NAME_banded) $(OBJECTS_banded)
#ar r $(NAME_sparse) $(OBJECTS_sparse)
	@if test -d ../lib/$(ARCH); then \
	echo "directory $(ARCH) already exists"; \
	else echo "create directory $(ARCH)" && mkdir ../lib/$(ARCH); fi;
	mv *.a ../lib/$(ARCH)
	rm -f $(OBJECTS) $(OBJECTS_sparse) $(OBJECTS_banded) $(OBJECTS_dense)

pfeast: $(OBJECTS_p) $(OBJECTS_psparse) $(OBJECTS_pbanded) $(OBJECTS_pdense)
	ar cr $(PNAME) $(OBJECTS_psparse) $(OBJECTS_pbanded) $(OBJECTS_pdense) $(OBJECTS_p) 
	@if test -d ../lib/$(ARCH); then \
	echo "directory $(ARCH) already exists"; \
	else echo "create directory $(ARCH)" && mkdir ../lib/$(ARCH); fi;
	mv *.a ../lib/$(ARCH)
	rm -f $(OBJECTS_p) $(OBJECTS_psparse) $(OBJECTS_pbanded) $(OBJECTS_pdense)  


clean:	
	rm ../lib/$(ARCH)/$(NAME)
	rm ../lib/$(ARCH)/$(PNAME)
#touch:
#	@touch $(CODE90) $(CODE90_dense) $(CODE90_banded) $(CODE90_sparse)
