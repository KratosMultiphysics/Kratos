############################################################################
#  PATHs to Compilers and libraries
############################################################################ 

############################################################################
########## Select C compiler and flags (needed for C and/or MPI-C examples)
#############################################################################
### Default
CC     = gcc     ##### FEAST library compiled with option 2 (noruntime dependency)
CFLAGS = -O3  

### other possibilities (examples)
#CC     = gcc -lifcoremt ### FEAST compiled with option 1 and ifort
#CC     = gcc -lgfortran ### FEAST compiled with option 1 and gfortran 

#############################################################################
########## Select F90 compiler and flags (needed for Fortran and/or MPI examples)
############################################################################
### Default
FORT  =  ifort #-lgfortran    
FFLAGS = -O3 -openmp 

### other possibility (example)
#FORT  =  gfortan #-lifcoremt   
#FFLAGS = -O3 -ffree-line-length-none -ffixed-line-length-none

#########################################################################
########## Select MPI-C compiler and flags (for C-MPI examples)
########################################################################

####### Default (If FEAST-MPI compiled with MPICH2)
PCC  =  mpicc.mpich2 -cc=$(CC) 
PCFLAGS = $(CFLAGS)


####### other possibilities (examples)
### If FEAST-MPI compiled with Intel-MPI
#PCC  =  mpiicc -cc=$(CC) 
## If FEAST-MPI compiled with OpenMPI...requires shell environment variable "OMPI_CC=$(CC)"
#export OMPI_CC=$(CC) # for BASH shell
#PCC  =  mpicc.openmpi 


#########################################################################
########## Select MPI-F90 compiler and flags (for Fortran-MPI examples)
#########################################################################

####### Default (If FEAST-MPI compiled with MPICH2)
PFORT  =  mpif90.mpich2 -f90=$(FORT)
PFFLAGS = $(FFLAGS)


######### other possibilities (examples)
### If FEAST-MPI compiled with Intel-MPI
#PFORT  =  mpiifort -f90=$(FORT)
### If FEAST-MPI compiled with OpenMPI ...requires shell environment variable "OMPI_FC=$(FORT)"
#export OMPI_FC=$(FORT) # for BASH shell
#PFORT  =  mpif90.openmpi 


########################################################################################################
####### Select PATH to FEAST Library 
#######################################################################################################

#==============================================
# Specify PATH to FEAST (main directory)- path can be defined as shell variable
#==============================================
FEASTROOT=$(PWD)/../../../../

LOCLIBS = -L$(FEASTROOT)/lib/x64 

###### NAMES for FEAST-SMP (for F90 or C examples)
FEAST= -lfeast 
FEAST_DENSE= -lfeast_dense
FEAST_BANDED= -lfeast_banded 
FEAST_SPARSE= -lfeast_sparse 

####### NAMES for FEAST-MPI (For MPI-F90 or MPI-C examples) 
####### *ATTENTION* compilated FEAST library should be compatible with MPI definition in PFORT or PCC above 

## example for intel-mpi
#PFEAST= -lpfeast_impi
#PFEAST_DENSE= -lpfeast_dense_impi
#PFEAST_BANDED= -lpfeast_banded_impi
#PFEAST_SPARSE= -lpfeast_sparse_impi

## example for MPICH2 
PFEAST= -lpfeast#_mpich2
PFEAST_DENSE= -lpfeast_dense#_mpich2
PFEAST_BANDED= -lpfeast_banded#_mpich2
PFEAST_SPARSE= -lpfeast_sparse#_mpich2

## example for OpenMPI
#PFEAST= -lpfeast_openmpi
#PFEAST_DENSE= -lpfeast_dense_openmpi
#PFEAST_BANDED= -lpfeast_banded_openmpi
#PFEAST_SPARSE= -lpfeast_sparse_openmpi

###############################################################  
######## Select PATH to INCLUDE FEAST files (for C or MPI-C examples)
##############################################################
INCL=-I/$(FEASTROOT)/include



#################################################################################################
$######## Select libraries 
#########                    MKL-PARDISO (mandatory only if sparse examples are compiled)
#########                    LAPACK (mandatory) 
#########                    BLAS (mandatory)
######### Remark: * Intel MKL includes Pardiso, Lapack and Blas
#########         * FEAST already includes the SPIKE-SMP solver for banded interfaces
#################################################################################################

########### Default (using Intel- MKL for all)

### Fortran examples
# with intel fortran
FLIBS= -mkl 
# with gfortran
#FLIBS=-Wl,--start-group -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core  -Wl,--end-group -liomp5 -lpthread

### C examples
# with gcc  
CLIBS=-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -Wl,--end-group -limf -lsvml -lirc -liomp5 -lpthread -lm
