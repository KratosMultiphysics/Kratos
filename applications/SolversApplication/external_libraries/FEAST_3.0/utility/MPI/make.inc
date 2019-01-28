###########################################
# USER SPECIFICATIONS #####################
###########################################

#==============================================
# Specify Compiler
#==============================================

## Fortran compiler
FC = ifort # -lgfortran  #### add flag if FEAST library compiled with Option-1 (pure gfortran) 
#FC=gfortan # -lifcoremt #### add flag if FEAST library compiled with Option-1 (pure ifort)

## MPI implementation -- Many possibilities to choose from (examples)
#PFC  =  mpiifort -f90=$(FC)
PFC  =  mpif90.mpich2 -f90=$(FC)
#PFC  =  mpif90.mpich -f90=$(FC)
#PFC  =  mpif90.openmpi # with shell variable "OMPI_FC=$(FC)" 

#==============================================
# Specify flags
#==============================================
FLAGS =-O3 -openmp  

#==============================================
# Specify PATH to FEAST (main directory)- can be defined as shell variable
#==============================================
FEASTROOT=$(PWD)/../../

#==============================================
# Specify libraries path (FEAST+PARDISO+LAPACK+BLAS) 
#==============================================
LOCLIBS = -L/$(FEASTROOT)/lib/x64   

#==============================================================
# Specify all the libraries feast -  mkl==pardiso -  lapack -blas  
#==============================================================
LIB=$(LOCLIBS) -lpfeast_sparse -lpfeast -mkl

###########################################
# END USER SPECIFICATIONS #################
###########################################
