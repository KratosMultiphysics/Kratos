
##### User configurable options #####

MAKE                 = make --no-print-directory

CPP                  = gcc -E
CC                   = gcc
CLINKER              = gcc
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
PROF_LIBS            = -L$(MPE_LIBDIR) -lmpe $(MPI_LIBS)
LOG_LIBS             = -L$(MPE_LIBDIR) -llmpe -lmpe $(MPI_LIBS)
TRACE_LIBS           = -L$(MPE_LIBDIR) -ltmpe $(MPI_LIBS)
F2CMPI_LIBS          = -lmpe_f2cmpi
PROF_FLIBS           = -L$(MPE_LIBDIR) \
                       $(F2CMPI_LIBS) -lmpe \
                       $(MPI_LIBS)
LOG_FLIBS            = -L$(MPE_LIBDIR) \
                       $(F2CMPI_LIBS) -llmpe -lmpe \
                       $(MPI_LIBS)
#
C_INCS               = -I$(MPE_INCDIR)
F_INCS               = -I$(MPE_INCDIR)
LDFLAGS              =  
MPE_BUILD_FORTRAN2C  = yes

### End User configurable options ###

SHELL         = /bin/sh
RM            = rm
MV            = mv
CP            = cp

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

C_EXECS       = cpilog cpilog_pack srtest comm1_isr
EXTRA_C_EXECS = cpi_trace cpi_log ncpilog ncpilog_pack iotest comm1_isr_loop log_cost
F_EXECS       = fpilog fpilog_pack fcomm1_isr
EXTRA_F_EXECS = nfpilog_pack

.SUFFIXES: .c .f .no .o 

.c.o:
	$(MPI_CC) $(MPICC_FLAGS) -c $<

.f.o:
	$(MPI_F77) $(MPIFC_FLAGS) -c $<

.c.no:
	$(MPI_CC) -DNO_MPI_LOGGING $(MPICC_FLAGS) -c $<
	$(MV) -f $*.o $*.no

.F.o:
	$(RM) -f mpe2tmp.f mpe2tmp.o
	$(CPP) $< > mpe2tmp.f && \
	$(MPI_F77) $(MPIFC_FLAGS) -c mpe2tmp.f
	$(MV) -f mpe2tmp.o $*.o

.F.no:
	$(RM) -f mpe2tmp.f mpe2tmp.o
	$(CPP) -DNO_MPI_LOGGING $< > mpe2tmp.f && \
	$(MPI_F77) $(MPIFC_FLAGS) -c mpe2tmp.f
	$(MV) -f mpe2tmp.o $*.no

default: $(C_EXECS)
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	      for pgm in $(F_EXECS) ; do \
	          $(MAKE) $$pgm ; \
	      done ; \
	  fi

all: default $(EXTRA_C_EXECS)
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	      for pgm in $(EXTRA_F_EXECS) ; do \
	          $(MAKE) $$pgm ; \
	      done ; \
	  fi

cpi_trace: cpi.o
	$(MPE_CLINKER)mpitrace $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(TRACE_LIBS) $(DEBUG_LIBS) -lm

cpi_log: cpi.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

cpilog: cpilog.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm 
#	$(MPI_CLINKER) $(MPICC_LDFLAGS)  -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm 

cpilog_pack: cpilog_pack.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm 
#	$(MPI_CLINKER) $(MPICC_LDFLAGS)  -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm 

ncpilog: cpilog.no
	$(CP) -f $? mpe2tmp.o
	$(MPE_CLINKER)log $(MPICC_LDFLAGS) -o $@ mpe2tmp.o $(DEBUG_LIBS) -lm 
	@-$(RM) -f mpe2tmp.o
#	$(MPI_CLINKER) $(MPICC_LDFLAGS)  -o $@ $? $(PROF_LIBS) $(DEBUG_LIBS) -lm 

ncpilog_pack: cpilog_pack.no
	$(CP) -f $? mpe2tmp.o
	$(MPE_CLINKER)log $(MPICC_LDFLAGS)  -o $@ mpe2tmp.o $(DEBUG_LIBS) -lm 
	@-$(RM) -f mpe2tmp.o
#	$(MPI_CLINKER) $(MPICC_LDFLAGS)  -o $@ $? $(PROF_LIBS) $(DEBUG_LIBS) -lm 

fpilog: fpilog.o
	$(MPE_FLINKER)mpilog $(MPIFC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(LOG_FLIBS) $(DEBUG_LIBS) -lm

fpilog_pack: fpilog_pack.o
	$(MPE_FLINKER)mpilog $(MPIFC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(LOG_FLIBS) $(DEBUG_LIBS) -lm

nfpilog_pack: fpilog_pack.no
	$(CP) -f $? mpe2tmp.o
	$(MPE_FLINKER)log $(MPIFC_LDFLAGS) -o $@ mpe2tmp.o $(DEBUG_LIBS) -lm
	@-$(RM) -f mpe2tmp.o
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(PROF_FLIBS) $(DEBUG_LIBS) -lm

srtest: srtest.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

pthread_sendrecv: pthread_sendrecv.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

pthread_allreduce: pthread_allreduce.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

iotest: iotest.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

openmp_mpilog: openmp_mpilog.o
	$(MPE_FLINKER)mpilog $(MPIFC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(LOG_FLIBS) $(DEBUG_LIBS) -lm

comm1_isr: comm1_isr.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

comm1_isr_loop: comm1_isr_loop.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS)  -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

fcomm1_isr: fcomm1_isr.o
	$(MPE_FLINKER)mpilog $(MPIFC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_FLINKER) $(MPIFC_LDFLAGS) -o $@ $? $(LOG_FLIBS) $(DEBUG_LIBS) -lm

comm2_connect_accept: comm2_connect_accept.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

comm2_spawn_parent: comm2_spawn_parent.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

comm2_spawn_child: comm2_spawn_child.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(MPICC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

comm2_join: comm2_join.o
	$(MPE_CLINKER)mpilog $(MPICC_LDFLAGS) -o $@ $? $(DEBUG_LIBS) -lm
#	$(MPI_CLINKER) $(CC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

log_cost: log_cost.c
	$(CC) $(CC_FLAGS) -DCLOG_NOMPI -c $?
	$(CLINKER) $(CC_LDFLAGS) -DCLOG_NOMPI -o $@ log_cost.o -L$(MPE_LIBDIR) -lmpe_nompi

linktest: linktest_C

linktest_C: 
	@-$(MPE_SBINDIR)/mpetestlink cpi_trace \
	  "*** Link C program with the MPI tracing library .........................."
	@-$(MPE_SBINDIR)/mpetestlink cpi_log \
	  "*** Link C program with the MPI logging library .........................."
	@-$(MPE_SBINDIR)/mpetestlink cpilog \
	  "*** Link C program with the MPI and manual logging libraries ............."
	@-echo ; $(MAKE) clean

linktest_F77: 
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	    $(MPE_SBINDIR)/mpetestlink fpilog \
	    "*** Link F77 program with the MPI and manual logging libraries ..........." ; \
	    $(MAKE) clean ; \
	    echo ; \
	  fi


runtest_C: 
	@-$(MPE_SBINDIR)/mpetestexec cpi_trace \
	  "*** Test C program with the MPI tracing library .........................."
	@-$(MPE_SBINDIR)/mpetestexeclog cpi_log \
	  "*** Test C program with the MPI logging library .........................."
	@-$(MPE_SBINDIR)/mpetestexeclog cpilog \
	  "*** Test C program with the MPI and manual logging libraries ............."
	@-echo ; $(MAKE) clean

runtest_F77: 
	@-if [ "$(MPE_BUILD_FORTRAN2C)" = "yes" ] ; then \
	      $(MPE_SBINDIR)/mpetestexeclog fpilog \
	      "*** Test F77 program with the MPI and manual logging libraries ..........." ; \
	      $(MAKE) clean ; \
	      echo ; \
	  fi

clean:
	@-$(RM) -f work.pc work.pcl
	@-$(RM) -f *.o *.no *~ PI*
	@-$(RM) -f mpe2tmp.f
	@-$(RM) -f $(C_EXECS) $(F_EXECS)
	@-$(RM) -f $(EXTRA_C_EXECS) $(EXTRA_F_EXECS)
	@-$(RM) -f *.clog2 *.slog2

distclean: clean
	@-$(RM) -f Makefile
