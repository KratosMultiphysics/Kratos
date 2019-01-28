include ./make.inc

help:
	@echo
	@echo "Usage: make {all allC allF pallF pallC rall rallC rallF clean cleanC cleanF pcleanF pcleanC}"
	@echo
	@echo "************* WARNING ************************** "
	@echo "1- Before use: Correct accordingly compilers and libraries paths" 
	@echo "    in file make.inc !!!!"
	@echo "2- This Makefile provides a (general) complete testing "
	@echo "    DENSE, BANDED, or SPARSE examples can also be compiled separately"
	@echo "    from sub-directories "
	@echo "3- MPI examples can only be run from sub-directories"
	@echo
	@echo "************* Compilation ************************ "
	@echo "make all    -- compile all C and F90 examples "
	@echo "make allC   -- compile all C examples "
	@echo "make allF   -- compile all F90 examples "
	@echo "make pallF  -- compile all MPI-F90 examples "
	@echo "make pallC  -- compile all MPI-C examples "
	@echo "************* Both Compilation and Run ********** "
	@echo "make rall   -- compile and run all C and F90 examples "
	@echo "make rallC  -- compile and run all C examples "
	@echo "make rallF  -- compile and run all F90 examples "
	@echo "************* Cleaning ************************** "
	@echo "make clean  -- clean all C and F90 examples"
	@echo "make cleanC -- clean all C examples"
	@echo "make cleanF -- clean all F90 examples"
	@echo "make pcleanF-- clean all MPI-F90 examples"
	@echo "make pcleanC-- clean all MPI-C examples"
	@echo ""

#########################

allC:
	cd ./Hermitian/C/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all; cd ../../../Non-Hermitian/C/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\

allF:
	cd ./Hermitian/Fortran/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all; cd ../../../Non-Hermitian/Fortran/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\


pallF:
	cd ./Hermitian/Fortran-MPI/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all; cd ../../../Non-Hermitian/Fortran-MPI/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\


pallC:
	cd ./Hermitian/C-MPI/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all; cd ../../../Non-Hermitian/C-MPI/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\



all: allF allC pallF pallC


cleanC:
	cd ./Hermitian/C/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean; cd ../../../Non-Hermitian/C/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\

cleanF:
	cd ./Hermitian/Fortran/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean; cd ../../../Non-Hermitian/Fortran/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\

pcleanF:
	cd ./Hermitian/Fortran-MPI/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean; cd ../../../Non-Hermitian/Fortran-MPI/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\

pcleanC:
	cd ./Hermitian/C-MPI/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean; cd ../../../Non-Hermitian/C-MPI/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\


clean: cleanF cleanC pcleanF pcleanC




rallC:
	cd ./Hermitian/C/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall; cd ../../../Non-Hermitian/C/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall;\


rallF:
	cd ./Hermitian/Fortran/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall; cd ../../../Non-Hermitian/Fortran/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall;\



rall: rallF rallC
