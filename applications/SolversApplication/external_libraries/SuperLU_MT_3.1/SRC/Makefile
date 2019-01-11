# makefile for sparse supernodal LU, implemented in ANSI C, targeted for 
# shared memory machines (SMPs).
#
include ../make.inc

#######################################################################
#  This is the makefile to create a library for SuperLU.
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions
#       SLUSRC -- Single precision real SuperLU routines
#       DLUSRC -- Double precision real SuperLU routines
#       CLUSRC -- Single precision complex SuperLU routines
#       ZLUSRC -- Double precision complex SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

ALLAUX = superlu_timer.o dclock.o sp_ienv.o lsame.o xerbla.o \
	util.o pmemory.o qrnzcnt.o cholnzcnt.o await.o \
	get_perm_c.o mmd.o colamd.o sp_coletree.o \
	pxgstrf_scheduler.o sp_colorder.o \
	pxgstrf_mark_busy_descends.o \
	pxgstrf_pruneL.o pxgstrf_super_bnd_dfs.o \
	pxgstrf_relax_snode.o heap_relax_snode.o \
	pxgstrf_synch.o pxgstrf_finalize.o

SZAUX  = slamch.o
DZAUX  = dlamch.o

SLUSRC = sreadhb.o sreadrb.o smatgen.o \
	psgssv.o psgssvx.o sgstrs.o sgsrfs.o \
	sgscon.o slacon.o slangs.o sgsequ.o slaqgs.o spivotgrowth.o \
	psmemory.o psutil.o smyblas2.o \
	psgstrf.o psgstrf_init.o \
	psgstrf_thread.o psgstrf_thread_init.o psgstrf_thread_finalize.o \
	psgstrf_factor_snode.o psgstrf_snode_dfs.o psgstrf_snode_bmod.o \
	psgstrf_panel_dfs.o psgstrf_panel_bmod.o \
	psgstrf_copy_to_ucol.o psgstrf_pivotL.o \
	psgstrf_column_dfs.o psgstrf_column_bmod.o \
	psgstrf_bmod1D.o psgstrf_bmod2D.o \
	psgstrf_bmod1D_mv2.o psgstrf_bmod2D_mv2.o \
	ssp_blas2.o ssp_blas3.o

DLUSRC = dreadhb.o dreadrb.o dmatgen.o \
	pdgssv.o pdgssvx.o dgstrs.o dgsrfs.o \
	dgscon.o dlacon.o dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o \
	pdmemory.o pdutil.o dmyblas2.o \
	pdgstrf.o pdgstrf_init.o \
	pdgstrf_thread.o pdgstrf_thread_init.o pdgstrf_thread_finalize.o \
	pdgstrf_factor_snode.o pdgstrf_snode_dfs.o pdgstrf_snode_bmod.o \
	pdgstrf_panel_dfs.o pdgstrf_panel_bmod.o \
	pdgstrf_copy_to_ucol.o pdgstrf_pivotL.o \
	pdgstrf_column_dfs.o pdgstrf_column_bmod.o \
	pdgstrf_bmod1D.o pdgstrf_bmod2D.o \
	pdgstrf_bmod1D_mv2.o pdgstrf_bmod2D_mv2.o \
	dsp_blas2.o dsp_blas3.o 

CLUSRC = scomplex.o creadhb.o creadrb.o cmatgen.o scsum1.o icmax1.o \
	pcgssv.o pcgssvx.o cgstrs.o cgsrfs.o \
	cgscon.o clacon.o clangs.o cgsequ.o claqgs.o cpivotgrowth.o \
	pcmemory.o pcutil.o cmyblas2.o \
	pcgstrf.o pcgstrf_init.o \
	pcgstrf_thread.o pcgstrf_thread_init.o pcgstrf_thread_finalize.o \
	pcgstrf_factor_snode.o pcgstrf_snode_dfs.o pcgstrf_snode_bmod.o \
	pcgstrf_panel_dfs.o pcgstrf_panel_bmod.o \
	pcgstrf_copy_to_ucol.o pcgstrf_pivotL.o \
	pcgstrf_column_dfs.o pcgstrf_column_bmod.o \
	pcgstrf_bmod1D.o pcgstrf_bmod2D.o \
	pcgstrf_bmod1D_mv2.o pcgstrf_bmod2D_mv2.o \
	csp_blas2.o csp_blas3.o

ZLUSRC = dcomplex.o zreadhb.o zreadrb.o zmatgen.o dzsum1.o izmax1.o \
	pzgssv.o pzgssvx.o zgstrs.o zgsrfs.o \
	zgscon.o zlacon.o zlangs.o zgsequ.o zlaqgs.o zpivotgrowth.o \
	pzmemory.o pzutil.o zmyblas2.o \
	pzgstrf.o pzgstrf_init.o \
	pzgstrf_thread.o pzgstrf_thread_init.o pzgstrf_thread_finalize.o \
	pzgstrf_factor_snode.o pzgstrf_snode_dfs.o pzgstrf_snode_bmod.o \
	pzgstrf_panel_dfs.o pzgstrf_panel_bmod.o \
	pzgstrf_copy_to_ucol.o pzgstrf_pivotL.o \
	pzgstrf_column_dfs.o pzgstrf_column_bmod.o \
	pzgstrf_bmod1D.o pzgstrf_bmod2D.o \
	pzgstrf_bmod1D_mv2.o pzgstrf_bmod2D_mv2.o \
	zsp_blas2.o zsp_blas3.o 

all:    single double complex complex16

single: $(SLUSRC) $(ALLAUX) $(SZAUX)
	$(ARCH) $(ARCHFLAGS) ../lib/$(SUPERLULIB) $(SLUSRC) $(ALLAUX) $(SZAUX)
	$(RANLIB) ../lib/$(SUPERLULIB)

double: $(DLUSRC) $(ALLAUX) $(DZAUX)
	$(ARCH) $(ARCHFLAGS) ../lib/$(SUPERLULIB) $(DLUSRC) $(ALLAUX) $(DZAUX)
	$(RANLIB) ../lib/$(SUPERLULIB)

complex: $(CLUSRC) $(ALLAUX) $(SZAUX)
	$(ARCH) $(ARCHFLAGS) ../lib/$(SUPERLULIB) $(CLUSRC) $(ALLAUX) $(SZAUX)
	$(RANLIB) ../lib/$(SUPERLULIB)

complex16: $(ZLUSRC) $(ALLAUX) $(DZAUX)
	$(ARCH) $(ARCHFLAGS) ../lib/$(SUPERLULIB) $(ZLUSRC) $(ALLAUX) $(DZAUX)
	$(RANLIB) ../lib/$(SUPERLULIB)

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $(NOOPTS) $<
slamch.o: slamch.c ; $(CC) -c $(NOOPTS) $<
superlu_timer.o:  superlu_timer.c ; $(CC) -c $(PREDEFS) $(NOOPTS) $<
dclock.o:  dclock.c ; $(CC) -c $(PREDEFS) $(NOOPTS) $<

# Do not optimize this on Cray; related to 'volatile' variable.
await.o: await.c
	$(CC) -c $(NOOPTS) $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

clean:	
	rm -f *.o core ../lib/$(SUPERLULIB)
