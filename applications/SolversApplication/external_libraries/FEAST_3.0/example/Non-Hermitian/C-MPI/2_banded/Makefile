include ../../../make.inc

help:
	@echo "Usage: make {all rall clean} "
	@echo "make all    -- compile all C examples for feast_dense interface "
	@echo "make rall   -- compile and run all C examples for feast_dense interface "
	@echo "make clean  -- clean all C examples for feast_dense interface "
	@echo "!!!!Please correct accordingly compiler and libraries paths, change compiler options " 
	@echo " in file ../../make.inc !!!!"
	@echo

#==============================================================
# Include the LIB (feast and  lapack - blas)  
#==============================================================
LIB = $(LOCLIBS) $(PFEAST_BANDED) $(PFEAST) $(CLIBS)
#==============================================================
# List of codes to be compiled 
#==============================================================
EXAMPLES = pdriver_zfeast_sbev pdriver_zfeast_sbevx pdriver_dfeast_gbgv  
#==============================================================
# Compile Link Execute
#==============================================================
all: examples 

examples: 
	@echo $(EXAMPLES)
	@for file in $(EXAMPLES); \
	do \
		echo $(PCC)  $(PCFLAGS) $(INCL) -c $$file.c;\
		$(PCC)  $(PCFLAGS) $(INCL) -c $$file.c ;\
		echo $(PCC)   -o $$file $$file.o $(LIB) ;\
		$(PCC)   -o $$file $$file.o  $(LIB) ;\
	done

rall: 	examples
	@for file in $(EXAMPLES); \
	do \
                ./$$file; \
	done


#==========================================
# Clean up directory: delete object files 
#==========================================
clean: 
	-@rm  $(EXAMPLES) *.o
