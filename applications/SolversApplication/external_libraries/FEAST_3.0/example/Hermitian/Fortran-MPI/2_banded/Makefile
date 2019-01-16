include ../../../make.inc


help:
	@echo "Usage: make {all rall clean} "
	@echo "make all    -- compile/run all F90 examples for feast_banded interface "
	@echo "make clean  -- clean all F90 examples for feast_banded interface "
	@echo "!!!!Please correct accordingly compiler and libraries paths, change compiler options " 
	@echo " in file ../../make.inc !!!!"
	@echo

#==============================================================
# Include the LIB (feast and  lapack - blas)  
#==============================================================
LIB = $(LOCLIBS)  $(PFEAST_BANDED) $(PFEAST) $(FLIBS)
#==============================================================
# List of codes to be compiled 
#==============================================================
EXAMPLES = pdriver_sfeast_sbgv pdriver_dfeast_sbgv pdriver_cfeast_hbev pdriver_zfeast_hbev 3pdriver_zfeast_hbev
#==============================================================
# Compile Link Execute
#==============================================================
all: examples 


examples: 
	@echo $(EXAMPLES)
	@for file in $(EXAMPLES); \
	do \
		echo $(PFORT)  $(PFFLAGS) $(INCL) -c $$file.f90;\
		$(PFORT)  $(PFFLAGS) $(INCL) -c $$file.f90 ;\
		echo $(PFORT)   -o $$file $$file.o $(LIB) ;\
		$(PFORT)   -o $$file $$file.o  $(LIB) ;\
	done


#==========================================
# Clean up directory: delete object files 
#==========================================
clean: 
	-@rm  $(EXAMPLES) *.o
