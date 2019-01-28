###########################################
# USER SPECIFICATIONS #####################
###########################################

#==============================================
# Specify Compiler
#==============================================

## Fortran compiler
FC = ifort # -lgfortran  #### add flag if FEAST library compiled with Option-1 (pure gfortran) 
#FC=gfortan # -lifcoremt #### add flag if FEAST library compiled with Option-1 (pure ifort)

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
LIB=$(LOCLIBS) -lfeast_sparse -lfeast -mkl

###########################################
# END USER SPECIFICATIONS #################
###########################################
