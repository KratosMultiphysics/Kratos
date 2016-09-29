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
LIB = $(LOCLIBS) $(FEAST_DENSE) $(FEAST) $(CLIBS)
#==============================================================
# List of codes to be compiled 
#==============================================================
EXAMPLES = helloworld driver_sfeast_sygv driver_dfeast_sygv driver_cfeast_heev driver_zfeast_heev
#==============================================================
# Compile Link Execute
#==============================================================
all: examples 

examples: 
	@echo $(EXAMPLES)
	@for file in $(EXAMPLES); \
	do \
		echo $(CC)  $(CFLAGS) $(INCL) -c $$file.c;\
		$(CC)  $(CFLAGS) $(INCL) -c $$file.c ;\
		echo $(CC)   -o $$file $$file.o $(LIB) ;\
		$(CC)   -o $$file $$file.o  $(LIB) ;\
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
