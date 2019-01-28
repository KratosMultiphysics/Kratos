include ./make.inc


help:
	@echo "Usage: make {all clean} "
	@echo "===> Please edit file *make.inc*  " 
	@echo "===> then change accordingly compiler, options, libraries paths, etc."
	@echo

#============================================================
# COMPILE and LINK
#============================================================
EXAMPLES = pdriver_feast_sparse 


all: examples 


examples: 
	@for file in $(EXAMPLES); \
	do \
		echo $(PFC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		$(PFC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		echo $(PFC)   -o $$file $$file.o $(LIB) ;\
		$(PFC)   -o $$file $$file.o $(LIB) ;\
	done


#==========================================================
# Clean up directory: delete object files and modules
#==========================================================
clean: 
	-@rm  $(EXAMPLES) *.o
