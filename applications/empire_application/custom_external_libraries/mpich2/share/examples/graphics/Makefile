
##### User configurable options #####

MAKE                 = make --no-print-directory

MPI_CC               = /home/stefan/tmp/empire_package/mpich2/bin/mpicc
MPI_F77              = /home/stefan/tmp/empire_package/mpich2/bin/mpif77
MPI_CLINKER          = $(MPI_CC)
MPI_FLINKER          = $(MPI_F77)
MPE_CLINKER          = /home/stefan/tmp/empire_package/mpich2/bin/mpicc -mpe=
MPE_FLINKER          = /home/stefan/tmp/empire_package/mpich2/bin/mpif77 -mpe=
#
DEBUG_LIBS           = 
MPI_CFLAGS           = 
MPI_FFLAGS           = 
MPI_LIBS             = 
CFLAGS               =    -O2
FFLAGS               =   -O2
MPE_INCDIR           = /home/stefan/tmp/empire_package/mpich2/include
MPE_LIBDIR           = /home/stefan/tmp/empire_package/mpich2/lib
MPE_SBINDIR          = /home/stefan/tmp/empire_package/mpich2/sbin
ANIM_LIBS            = -L$(MPE_LIBDIR) -lampe -lmpe -lX11 $(MPI_LIBS)
C_XGRAPH_LIBS        = -L$(MPE_LIBDIR) -lmpe -lX11 $(MPI_LIBS)
F_XGRAPH_LIBS        = -L$(MPE_LIBDIR) -lmpe -lX11 $(MPI_LIBS)
#
C_INCS               = -I$(MPE_INCDIR)
F_INCS               = -I$(MPE_INCDIR)
LDFLAGS              =  
MPE_BUILD_FORTRAN2C  = yes

### End User configurable options ###

SHELL         = /bin/sh
RM            = rm 

VPATH=.:$(srcdir)

srcdir        = .

# MPE include path, x_INCS, goes in front of MPI include path, MPI_xFLAGS,
# in case MPI include path is in MPI_xFLAGS containing old MPE include path.
MPICC_FLAGS   = $(C_INCS) $(MPI_CFLAGS) $(MPE_COPTS)
MPIFC_FLAGS   = $(F_INCS) $(MPI_FFLAGS) $(MPE_FOPTS)
MPICC_LDFLAGS = $(LDFLAGS) $(MPI_CFLAGS) $(MPE_LDOPTS)
MPIFC_LDFLAGS = $(LDFLAGS) $(MPI_FFLAGS) $(MPE_LDOPTS)
EXECS         = cpi_anim cxgraphics

.SUFFIXES: .c .f .o 

.c.o:
	$(MPI_CC) $(MPICC_FLAGS) -c $<
.f.o:
	$(MPI_F77) $(MPIFC_FLAGS) -c $<

default: $(EXECS)
	-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	     $(MAKE) fxgraphics ; \
	 fi

cpi_anim: cpi.o
	$(MPE_CLINKER)mpianim $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm 
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(ANIM_LIBS) $(DEBUG_LIBS) -lm 

cxgraphics: cxgraphics.o
	$(MPE_CLINKER)graphics $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(C_XGRAPH_LIBS) $(DEBUG_LIBS) -lm

fxgraphics: fxgraphics.o
	$(MPE_FLINKER)graphics $(MPIFC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(F_XGRAPH_LIBS) $(DEBUG_LIBS) -lm

linktest: linktest_C

linktest_C: 
	@-$(MPE_SBINDIR)/mpetestlink cpi_anim \
	  "*** Link C program with the MPI animation library ........................"
	@-$(MPE_SBINDIR)/mpetestlink cxgraphics \
	  "*** Link C program with the X11 graphics library ........................."
	@-echo ; $(MAKE) clean

linktest_F77: 
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	      $(MPE_SBINDIR)/mpetestlink fxgraphics \
	      "*** Link F77 program with the X11 graphics library ......................." ; \
	      $(MAKE) clean ; \
	      echo ; \
	  fi

runtest_C: 
	@-$(MPE_SBINDIR)/mpetestexec cpi_anim \
	  "*** Test C program with the MPI animation library ........................"
	@-$(MPE_SBINDIR)/mpetestexec cxgraphics \
	  "*** Test C program with the X11 graphics library ........................."
	@-echo ; $(MAKE) clean

runtest_F77: 
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	      $(MPE_SBINDIR)/mpetestexec fxgraphics \
	      "*** Test F77 program with the X11 graphics library ......................." ; \
	      $(MAKE) clean ; \
	      echo ; \
	  fi

clean:
	@-$(RM) -f work.pc work.pcl
	@-$(RM) -f *.o *~ PI* $(EXECS)
	@-$(RM) -f cpi_anim cxgraphics fxgraphics

distclean: clean
	@-$(RM) -f Makefile
