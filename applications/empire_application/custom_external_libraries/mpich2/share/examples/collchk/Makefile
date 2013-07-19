
##### User configurable options #####

MAKE                 = make --no-print-directory

MPI_CC               = /home/stefan/tmp/empire_package/mpich2/bin/mpicc
MPI_F77              = /home/stefan/tmp/empire_package/mpich2/bin/mpif77
MPI_CLINKER          = $(MPI_CC)
MPI_FLINKER          = $(MPI_F77)
MPE_CLINKER          = /home/stefan/tmp/empire_package/mpich2/bin/mpicc -mpe=
MPE_FLINKER          = /home/stefan/tmp/empire_package/mpich2/bin/mpif77 -mpe=

DEBUG_LIBS           = 
MPI_CFLAGS           = 
MPI_FFLAGS           = 
MPI_LIBS             = 
CFLAGS               =    -O2
FFLAGS               =   -O2
MPE_INCDIR           = /home/stefan/tmp/empire_package/mpich2/include
MPE_LIBDIR           = /home/stefan/tmp/empire_package/mpich2/lib
MPE_SBINDIR          = /home/stefan/tmp/empire_package/mpich2/sbin
F2CMPI_LIBS          = -lmpe_f2cmpi
COLLCHK_LIBS         = -L$(MPE_LIBDIR) -lmpe_collchk -lmpe $(MPI_LIBS)
COLLCHK_FLIBS        = -L$(MPE_LIBDIR) \
                       $(F2CMPI_LIBS) -lmpe_collchk -lmpe \
                       $(MPI_LIBS)

C_INCS               = -I$(MPE_INCDIR)
F_INCS               = -I$(MPE_INCDIR)
LDFLAGS              =  
MPE_BUILD_FORTRAN2C  = yes

#### End User configurable options ###

SHELL         = /bin/sh
RM            = rm

VPATH=.:$(srcdir)

srcdir        = .

CC_FLAGS      = $(C_INCS) $(CFLAGS) $(MPE_COPTS)
FC_FLAGS      = $(F_INCS) $(FFLAGS) $(MPE_FOPTS)
CC_LDFLAGS    = $(LDFLAGS) $(CFLAGS) $(MPE_LDOPTS)
FC_LDFLAGS    = $(LDFLAGS) $(FFLAGS) $(MPE_LDOPTS)
# MPE include path, x_INCS, goes in front of MPI include path, MPI_xFLAGS,
# in case MPI include path is in MPI_xFLAGS containing old MPE include path.
MPICC_FLAGS   = $(C_INCS) $(MPI_CFLAGS) $(MPE_COPTS)
MPIFC_FLAGS   = $(F_INCS) $(MPI_FFLAGS) $(MPE_FOPTS)
MPICC_LDFLAGS = $(LDFLAGS) $(MPI_CFLAGS) $(MPE_LDOPTS)
MPIFC_LDFLAGS = $(LDFLAGS) $(MPI_FFLAGS) $(MPE_LDOPTS)

.SUFFIXES: .c .f .o 

C_SRC         = wrong_composite.c wrong_scatterv.c \
                wrong_int_byte.c wrong_root.c \
                composite.c scatterv.c \
                time_bcast.c time_allreduce.c \
                time_alltoallv.c time_alltoallv0.c
F_SRC         = wrong_reals.f
C_OBJS        = $(C_SRC:.c=.o)
F_OBJS        = $(F_SRC:.f=.o)
C_EXECS       = $(C_OBJS:.o=)

OPT_C_EXECS   = time_bcast_nochk time_allreduce_nochk \
                time_alltoallv0_nochk time_alltoallv_nochk

OPT_F_EXECS   = wrong_reals wrong_reals_nochk

.c.o:
	$(MPI_CC) $(MPICC_FLAGS) -c $<

.f.o:
	$(MPI_F77) $(MPIFC_FLAGS) -c $<

.o:
	$(MPE_CLINKER)mpicheck $(MPICC_LDFLAGS) -o $@ $<
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $< $(COLLCHK_LIBS)

time_bcast_nochk: time_bcast.o
	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(MPI_LIBS)

time_allreduce_nochk: time_allreduce.o
	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(MPI_LIBS)

time_alltoallv_nochk: time_alltoallv.o
	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(MPI_LIBS)

time_alltoallv0_nochk: time_alltoallv0.o
	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(MPI_LIBS)

wrong_int_byte: wrong_int_byte.o
	$(MPE_CLINKER)mpicheck $(MPICC_LDFLAGS) -o $@ $?
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(COLLCHK_LIBS)

wrong_reals: wrong_reals.o
	$(MPE_FLINKER)mpicheck $(MPIFC_LDFLAGS) -o $@ $?
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(COLLCHK_FLIBS)

wrong_reals_nochk: wrong_reals.o
	$(MPE_FLINKER)mpicheck $(MPIFC_LDFLAGS) -o $@ $?
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(MPI_LIBS)

default: $(C_OBJS) $(C_EXECS) $(OPT_C_EXECS)

linktest: linktest_C

linktest_C:
	@-$(MPE_SBINDIR)/mpetestlink wrong_int_byte \
	  "*** Link C program with the MPI collective/datatype checking library ....."
	@-echo ; $(MAKE) clean

linktest_F77:
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	    $(MPE_SBINDIR)/mpetestlink wrong_reals \
	    "*** Link F77 program with the MPI collective/datatype checking library ..." ; \
	    $(MAKE) clean ; \
	    echo ; \
	  fi


runtest_C:
	@-$(MPE_SBINDIR)/mpetestexec wrong_int_byte \
	  "*** Test C program with the MPI collective/datatype checking library ....."
	@-echo ; $(MAKE) clean

runtest_F77:
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	    $(MPE_SBINDIR)/mpetestexec wrong_reals \
	    "*** Test F77 program with the MPI collective/datatype checking library ..." ; \
	    $(MAKE) clean ; \
	    echo ; \
	  fi

clean:
	@-$(RM) -f $(C_OBJS) $(F_OBJS) *.o
	@-$(RM) -f $(C_EXECS) $(OPT_C_EXECS) $(OPT_F_EXECS)

distclean: clean
	@-$(RM) -f Makefile
