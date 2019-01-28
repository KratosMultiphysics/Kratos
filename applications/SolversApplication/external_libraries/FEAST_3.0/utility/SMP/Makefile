include ./make.inc


help:
	@echo "Usage: make {all clean} "
	@echo "===> Please edit file *make.inc*  " 
	@echo "===> then change accordingly compiler, options, libraries paths, etc."
	@echo

#============================================================
# COMPILE and LINK
#============================================================
EXAMPLES = driver_feast_sparse 


all: examples 


examples: 
	@for file in $(EXAMPLES); \
	do \
		echo $(FC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		$(FC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		echo $(FC)   -o $$file $$file.o $(LIB) ;\
		$(FC)   -o $$file $$file.o $(LIB) ;\
	done


#==========================================================
# Clean up directory: delete object files and modules
#==========================================================
clean: 
	-@rm  $(EXAMPLES) *.o
