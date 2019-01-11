/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 * -- SuperLU MT routine (version 2.2) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Last modified: August 18, 2014
 *
 */
#ifdef unix
#include <unistd.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "slu_mt_ddefs.h"


void superlu_abort_and_exit(char* msg)
{
    fprintf(stderr,"%s\n", msg);
    exit (-1);
}

void *superlu_malloc(size_t size)
{
    void *buf;
    buf = (void *) malloc(size);
    return (buf);
}

void superlu_free(void *addr)
{
    free (addr);
}

/* Deallocate the structure pointing to the actual storage of the matrix. */
void
Destroy_SuperMatrix_Store(SuperMatrix *A)
{
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperMatrix_Store ...\n");
#endif
}

/* A is of type Stype==NC */
void
Destroy_CompCol_Matrix(SuperMatrix *A)
{
    SUPERLU_FREE( ((NCformat *)A->Store)->rowind );
    SUPERLU_FREE( ((NCformat *)A->Store)->colptr );
    SUPERLU_FREE( ((NCformat *)A->Store)->nzval );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_Matrix ...\n");
#endif
}

/* A is of type Stype==NCP */
void
Destroy_CompCol_Permuted(SuperMatrix *A)
{
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colbeg );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colend );
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_Permuted ...\n");
#endif
}

/* A is of type Stype==NCP */
void
Destroy_CompCol_NCP(SuperMatrix *A)
{
    SUPERLU_FREE ( ((NCPformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colbeg );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colend );
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_NCP ...\n");
#endif
}

/* A is of type Stype==SC */
void
Destroy_SuperNode_Matrix(SuperMatrix *A)
{
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCformat *)A->Store)->sup_to_col );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperNode_Matrix ...\n");
#endif
}

/* A is of type Stype==SCP */
void
Destroy_SuperNode_SCP(SuperMatrix *A)
{
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colend );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperNode_SCP ...\n");
#endif
}

/*
 * Reset repfnz[*] for the current column 
 */
void
pxgstrf_resetrep_col(const int_t nseg, const int_t *segrep, int_t *repfnz)
{
    register int_t i, irep;
    
    for (i = 0; i < nseg; ++i) {
	irep = segrep[i];
	repfnz[irep] = EMPTY;
    }
}


/*
 * Count the total number of nonzeros in factors L and U,  and in the 
 * symmetrically reduced L. 
 */
void
countnz(const int_t n, int_t *xprune, int_t *nnzL, int_t *nnzU, GlobalLU_t *Glu)
{
    register int_t nsuper, fsupc, i, j, nnzL0, jlen, irep;
    register int_t nnzsup = 0;
    register int_t *xsup, *xsup_end, *xlsub, *xlsub_end, *supno;
	
    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    supno     = Glu->supno;
    *nnzU     = Glu->nextu;
    nnzL0     = 0;
    *nnzL     = 0;
    nsuper    = supno[n];

    if ( n <= 0 ) return;

    /* 
     * For each supernode ...
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jlen = xlsub_end[fsupc] - xlsub[fsupc];
	nnzsup += jlen * (xsup_end[i] - fsupc);
			  
	for (j = fsupc; j < xsup_end[i]; j++) {
	    *nnzL += jlen;
	    *nnzU += j - fsupc + 1;
	    jlen--;
	}
	irep = SUPER_REP(i);
	if ( SINGLETON(supno[irep]) )
	    nnzL0 += xprune[irep] - xlsub_end[irep];
	else 
	    nnzL0 += xprune[irep] - xlsub[irep];
    }

#if ( PRNTlevel==1 )
    printf(".. # supernodes = " IFMT "\n", nsuper+1);
    printf(".. # edges in symm-reduced L = " IFMT "\n", nnzL0);
    if ( Glu->dynamic_snode_bound )
      printf(".. # NZ in LUSUP " IFMT ", dynamic bound " IFMT ", utilization %.2f\n",
	     nnzsup, Glu->nextlu, (float)nnzsup/Glu->nextlu);
    else
      printf(".. # NNZ in LUSUP " IFMT ", static bound " IFMT ", utilization %.2f\n",
	     nnzsup, Glu->nzlumax, (float)nnzsup/Glu->nzlumax);
#endif
}



/*
 * Fix up the data storage lsub for L-subscripts. It reclaims the
 * storage for the adjancency lists of the pruned graph, and applies
 * row permuation to the row subscripts of matrix $L$.
 */
void
fixupL(const int_t n, const int_t *perm_r, GlobalLU_t *Glu)
{
    register int_t nsuper, fsupc, nextl, i, j, jstrt;
    register int_t *xsup, *xsup_end, *lsub, *xlsub, *xlsub_end;

    if ( n <= 1 ) return;

    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    lsub      = Glu->lsub;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    nsuper    = Glu->supno[n];
    nextl     = 0;
    
    /* 
     * For each supernode ...
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jstrt = xlsub[fsupc];
	xlsub[fsupc] = nextl;
	for (j = jstrt; j < xlsub_end[fsupc]; j++) {
	    lsub[nextl] = perm_r[lsub[j]]; /* Now indexed into P*A */
	    nextl++;
  	}
	xlsub_end[fsupc] = nextl;
    }
    xlsub[n] = nextl;

#if ( PRNTlevel==1 )
    printf(".. # edges in supernodal graph of L = " IFMT "\n", nextl);
    fflush(stdout);
#endif
}

/*
 * Print all definitions to be used by CPP.
 */
int cpp_defs()
{
    printf("CPP Defs:\n");
#ifdef PRNTlevel
    printf("\tPRNTlevel=\n", PRNTlevel);
#endif 
#ifdef DEBUGlevel
    printf("\tDEBUGlevel=%d\n", DEBUGlevel);
#endif 
#ifdef PROFILE
    printf("\tPROFILE\n");
#endif
#ifdef PREDICT_OPT
    printf("\tPREDICT_OPT\n");
#endif
#ifdef USE_VENDOR_BLAS
    printf("\tUSE_VENDOR_BLAS\n");
#endif
#ifdef GEMV2
    printf("\tGEMV2\n");
#endif
#ifdef SCATTER_FOUND
    printf("\tSCATTER_FOUND\n");
#endif

    return 0;
}

/*
 * Compress the data storage LUSUP[*] for supernodes. It removes the
 * memory holes due to untightness of the upper bounds by A'A-supernode.
 */
void
compressSUP(const int_t n, GlobalLU_t *Glu)
{
    register int_t nextlu, i, j, jstrt;
    int_t *xlusup, *xlusup_end;
    double *lusup;

    if ( n <= 1 ) return;

    lusup     = Glu->lusup;
    xlusup    = Glu->xlusup;
    xlusup_end= Glu->xlusup_end;
    nextlu     = 0;
    
    for (j = 0; j < n; ++j) {
	jstrt = xlusup[j];
	xlusup[j] = nextlu;
	for (i = jstrt; i < xlusup_end[j]; ++i, ++nextlu)
	    lusup[nextlu] = lusup[i];
	xlusup_end[j] = nextlu;
    }
    xlusup[n] = nextlu;
    printf("\tcompressSUP() nextlu" IFMT "\n", nextlu);
}

int_t check_mem_leak(char *where)
{
#ifdef unix
    void *addr;
    addr = (void *)sbrk(0);
    printf("\tsbrk(0) %s: addr = %ld\n", where, (size_t) addr);
#endif
    return 0;
}

/*
 * Diagnostic print of segment info after pdgstrf_panel_dfs().
 */
void print_panel_seg(int_t n, int_t w, int_t jcol, int_t nseg, 
		     int_t *segrep, int_t *repfnz)
{
    int_t j, k;
    
    for (j = jcol; j < jcol+w; j++) {
	printf("\tcol" IFMT ":\n", j);
	for (k = 0; k < nseg; k++)
	    printf("\t\tseg"IFMT ", segrep"IFMT ", repfnz"IFMT "\n", k, 
			segrep[k], repfnz[(j-jcol)*n + segrep[k]]);
    }
}

/*
 * Allocate storage for various statistics.
 */
void
StatAlloc(const int_t n, const int_t nprocs, const int_t panel_size, 
	  const int_t relax, Gstat_t *Gstat)
{
    register int_t w;

    w = SUPERLU_MAX( panel_size, relax ) + 1;
    Gstat->panel_histo = intCalloc(w);
    Gstat->utime = (double *) SUPERLU_MALLOC(NPHASES * sizeof(double));
    Gstat->ops   = (flops_t *) SUPERLU_MALLOC(NPHASES * sizeof(flops_t));
    
    if ( !(Gstat->procstat =
	   (procstat_t *) SUPERLU_MALLOC(nprocs*sizeof(procstat_t))) )
	SUPERLU_ABORT( "SUPERLU_MALLOC failed for procstat[]" );

#if (PRNTlevel==1)
    printf(".. StatAlloc(): n " IFMT ", nprocs " IFMT ", panel_size " IFMT ", relax " IFMT "\n",
		n, nprocs, panel_size, relax);
#endif
#ifdef PROFILE    
    if ( !(Gstat->panstat =
	   (panstat_t*) SUPERLU_MALLOC(n * sizeof(panstat_t))) )
	SUPERLU_ABORT( "SUPERLU_MALLOC failed for panstat[]" );
    Gstat->panhows = intCalloc(3);
    Gstat->height = intCalloc(n+1);
    if ( !(Gstat->flops_by_height =
	   (float *) SUPERLU_MALLOC(n * sizeof(float))) )
	SUPERLU_ABORT("SUPERLU_MALLOC failed for flops_by_height[]");
    
#endif
    
#ifdef PREDICT_OPT
    if ( !(cp_panel = (cp_panel_t *) SUPERLU_MALLOC(n * sizeof(cp_panel_t))) )
	SUPERLU_ABORT( "SUPERLU_MALLOC failed for cp_panel[]" );
    if ( !(desc_eft = (desc_eft_t *) SUPERLU_MALLOC(n * sizeof(desc_eft_t))) )
	SUPERLU_ABORT( "SUPERLU_MALLOC failed for desc_eft[]" );
    cp_firstkid = intMalloc(n+1);
    cp_nextkid = intMalloc(n+1);
#endif
    
}

/*
 * Initialize various statistics variables.
 */
void
StatInit(const int_t n, const int_t nprocs, Gstat_t *Gstat)
{
    register int i;
    
    for (i = 0; i < NPHASES; ++i) {
	Gstat->utime[i] = 0;
	Gstat->ops[i] = 0;
    }
    
    for (i = 0; i < nprocs; ++i) {
	Gstat->procstat[i].panels = 0;
	Gstat->procstat[i].fcops = 0.0;
	Gstat->procstat[i].skedwaits = 0;
	Gstat->procstat[i].skedtime = 0.0;
	Gstat->procstat[i].cs_time = 0.0;
	Gstat->procstat[i].spintime = 0.0;
	Gstat->procstat[i].pruned = 0;
	Gstat->procstat[i].unpruned = 0;
    }

#ifdef PROFILE    
    for (i = 0; i < n; ++i) {
	Gstat->panstat[i].fctime = 0.0;
	Gstat->panstat[i].flopcnt = 0.0;
	Gstat->panstat[i].pipewaits = 0;
	Gstat->panstat[i].spintime = 0.0;
	Gstat->flops_by_height[i] = 0.0;
    }
    for (i = 0; i < 3; ++i) Gstat->panhows[i] = 0;
    Gstat->dom_flopcnt = 0.;
    Gstat->flops_last_P_panels = 0;
#endif
    
#ifdef PREDICT_OPT
    for (i = 0; i < n; ++i)
	cp_panel[i].est = cp_panel[i].pdiv = 0;
#endif
    
#if ( PRNTlevel==1 )
    printf(".. StatInit(): n " IFMT ", nprocs " IFMT "\n", n, nprocs);
#endif
}


/* Print timings used in factorization and solve. */
void
PrintStat(Gstat_t *Gstat)
{
    double         *utime;
    flops_t        *ops;

    utime = Gstat->utime;
    ops   = Gstat->ops;
    printf("Factor time  = %8.2f\n", utime[FACT]);
    if ( utime[FACT] != 0.0 )
      printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	     ops[FACT]*1e-6/utime[FACT]);

    printf("Solve time   = %8.2f\n", utime[SOLVE]);
    if ( utime[SOLVE] != 0.0 )
      printf("Solve flops = %e\tMflops = %8.2f\n", ops[SOLVE],
	     ops[SOLVE]*1e-6/utime[SOLVE]);

}

void
StatFree(Gstat_t *Gstat)
{
    SUPERLU_FREE (Gstat->panel_histo);
    SUPERLU_FREE (Gstat->utime);
    SUPERLU_FREE (Gstat->ops);
    SUPERLU_FREE (Gstat->procstat);
    
#ifdef PROFILE
    SUPERLU_FREE (Gstat->panstat);
    SUPERLU_FREE (Gstat->panhows);
    SUPERLU_FREE (Gstat->height);
    SUPERLU_FREE (Gstat->flops_by_height);
#endif

#ifdef PREDICT_OPT
    SUPERLU_FREE (Gstat->cp_panel);
    SUPERLU_FREE (Gstat->desc_eft);
    SUPERLU_FREE (Gstat->cp_firstkid);
    SUPERLU_FREE (Gstat->cp_nextkid);
#endif    

#if (PRNTlevel==1)
    printf(".. StatFree(): Free Stat variables.\n");
#endif
}

flops_t
LUFactFlops(Gstat_t *Gstat)
{
    return (Gstat->ops[FACT]);
}

flops_t
LUSolveFlops(Gstat_t *Gstat)
{
    return (Gstat->ops[SOLVE]);
}



/* 
 * Fills an integer array with a given value.
 */
void ifill(int_t *a, int_t alen, int_t ival)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = ival;
}



/* 
 * Get the statistics of the supernodes 
 */
#define NBUCKS 10
static 	int_t	max_sup_size;

void super_stats(int_t nsuper, int_t *xsup, int_t *xsup_end)
{
    register int_t nsup1 = 0;
    int_t          i, isize, whichb, bl, bh;
    int_t          bucket[NBUCKS];

    max_sup_size = 0;

    /* Histogram of the supernode sizes */
    ifill (bucket, NBUCKS, 0);

    for (i = 0; i <= nsuper; i++) {
        isize = xsup_end[i] - xsup[i];
	if ( isize == 1 ) nsup1++;
	if ( max_sup_size < isize ) max_sup_size = isize;	
        whichb = (float) isize / max_sup_size * NBUCKS;
        if (whichb >= NBUCKS) whichb = NBUCKS - 1;
        bucket[whichb]++;
    }
    
    printf("** Supernode statistics:\n\tno of supernodes = " IFMT "\n", nsuper+1);
    printf("\tmax supernode size = " IFMT "\n", max_sup_size);
    printf("\tno of size 1 supernodes = " IFMT "\n", nsup1);

    printf("\tHistogram of supernode size:\n");
    for (i = 0; i < NBUCKS; i++) {
        bl = (float) i * max_sup_size / NBUCKS;
        bh = (float) (i+1) * max_sup_size / NBUCKS;
        printf("\t" IFMT "-" IFMT "\t\t" IFMT "\n", bl+1, bh, bucket[i]);
    }

}

void panel_stats(int_t n, int_t max_w, int_t* in_domain, Gstat_t *Gstat)
{
    register int_t i, w;
    float *histo_flops, total;

    histo_flops = (float *) SUPERLU_MALLOC( max_w * sizeof(float) );

    for (i = 0; i < max_w; ++i) histo_flops[i] = 0;
    total = 0;
    for (i = 0; i < n; i += w) {
	w = Gstat->panstat[i].size;
	if ( in_domain[i] != TREE_DOMAIN ) {
	    histo_flops[w - 1] += Gstat->panstat[i].flopcnt;
	    total += Gstat->panstat[i].flopcnt;
	}
    }

    if ( total != 0.0 ) {
	printf("** Panel & flops distribution: nondomain flopcnt %e\n", total);
	for (i = 1; i <= max_w; i++)
	    printf("\t" IFMT "\t" IFMT "\t%e (%.2f)\n", i, Gstat->panel_histo[i],
		   histo_flops[i-1], histo_flops[i-1]/total);
    }
    SUPERLU_FREE (histo_flops);
}



float SpaSize(int_t n, int_t np, float sum_npw)
{
    return (sum_npw*8 + np*8 + n*4)/1024.;
}

float DenseSize(int_t n, float sum_nw)
{
    return (sum_nw*8 + n*8)/1024.;;
}


/*
 * Check whether repfnz[] == EMPTY after reset.
 */
void check_repfnz(int_t n, int_t w, int_t jcol, int_t *repfnz)
{
    int_t jj, k;

    for (jj = jcol; jj < jcol+w; jj++) 
	for (k = 0; k < n; k++)
	    if ( repfnz[(jj-jcol)*n + k] != EMPTY ) {
		fprintf(stderr, "col " IFMT ", repfnz_col[" IFMT "] = " IFMT "\n", jj,
			k, repfnz[(jj-jcol)*n + k]);
		SUPERLU_ABORT("repfnz[] not empty.");
	    }
}


int_t PrintInt10(char *name, int_t len, int_t *x)
{
    register int_t i;
    
    printf("(len=" IFMT ") %s:", len, name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) printf("\n[" IFMT "-" IFMT "]", i, i+9);
	printf(IFMT, x[i]);
    }
    printf("\n");
    return 0;
}

/* Print a summary of the testing results. */
void
PrintSumm(char *type, int_t nfail, int_t nrun, int_t nerrs)
{
    if ( nfail > 0 )
	printf("%3s driver: " IFMT " out of " IFMT " tests failed to pass the threshold\n",
	       type, nfail, nrun);
    else
	printf("All tests for %3s driver passed the threshold (" IFMT " tests run)\n", type, nrun);

    if ( nerrs > 0 )
	printf(IFMT " error messages recorded\n", nerrs);
}


/* Print the adjacency list for graph of L, including the pruned graph,
   graph of U, and L supernodes partition */
int_t PrintGLGU(int_t n, int_t *xprune, GlobalLU_t *Glu)
{
    register int_t nsuper = Glu->nsuper;
    PrintInt10("LSUB", Glu->xlsub_end[n-1], Glu->lsub);
    PrintInt10("XLSUB", n, Glu->xlsub);
    PrintInt10("XLSUB_END", n, Glu->xlsub_end);
    PrintInt10("XPRUNE", n, xprune);
    PrintInt10("USUB", Glu->xusub_end[n-1], Glu->usub);
    PrintInt10("XUSUB", n, Glu->xusub);
    PrintInt10("XUSUB_END", n, Glu->xusub_end);
    PrintInt10("SUPNO", n, Glu->supno);
    PrintInt10("XSUP", nsuper+1, Glu->xsup);
    PrintInt10("XSUP_END", nsuper+1, Glu->xsup_end);
    return 0;
}

#if 0
/*
 * Print the statistics of the relaxed snodes for matlab process
 */
void relax_stats(int_t start, int_t end, int_t step)
{
    FILE *fp;
    int i;

    fp = fopen("relax.m", "w");
    
    fprintf(fp,"relax = [\n");
    for (i = start; i <= end; i += step) fprintf(fp, "%d ", i);
    fprintf(fp, "];\n");

    fprintf(fp, "fctime = [\n");
    for (i = start; i <= end; i += step) 
	fprintf(fp, "%15.8e\n ", stat_relax[i].fctime);
    fprintf(fp, "];\n");

    fprintf(fp, "mflops = [\n");
    for (i = start; i <= end; i += step)
	fprintf(fp, "%15.8e\n ", (float)stat_relax[i].flops / 1e6);
    fprintf(fp, "];\n");

    fprintf(fp, "mnzs = [\n");
    for (i = start; i <= end; i += step)
 	fprintf(fp, "%15.8e\n ", stat_relax[i].nzs / 1e6);
    fprintf(fp, "];\n");

    fclose(fp);
}

/*
 * Obtain the distribution of time/flops/nzs on the snode size.
 */
void snode_profile(int_t nsuper, int_t *xsup)
{
    FILE *fp;
    int_t i, j;
    int_t ssize;

    if ( !(stat_snode = (stat_snode_t *) SUPERLU_MALLOC((max_sup_size+1) *
	sizeof(stat_snode_t))) ) ABORT("SUPERLU_MALLOC fails for stat_snode[].");

    for (i = 0; i <= max_sup_size; i++) {
	stat_snode[i].ncols = 0;
	stat_snode[i].flops = 0;
	stat_snode[i].nzs = 0;
	stat_snode[i].fctime = 0.0;
    }	

    for (i = 0; i <= nsuper; i++) {

	ssize = xsup[i+1] - xsup[i];   
	stat_snode[ssize].ncols += ssize;

        for (j=xsup[i]; j<xsup[i+1]; j++) { 
	    stat_snode[ssize].flops += stat_col[j].flops;	    
	    stat_snode[ssize].nzs += stat_col[j].nzs;	    
	    stat_snode[ssize].fctime += stat_col[j].fctime;	    
	}

    }

    fp = fopen("snode.m", "w");
    
    fprintf(fp, "max_sup_size = " IFMT ";\n", max_sup_size);

    fprintf(fp,"ncols = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, IFMT "  ", stat_snode[i].ncols);
    fprintf(fp, "];\n");

    fprintf(fp, "fctime = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", stat_snode[i].fctime);
    fprintf(fp, "];\n");

    fprintf(fp, "mflops = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", (float) stat_snode[i].flops / 1e6);
    fprintf(fp, "];\n");

    fprintf(fp, "mnzs = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", (float) stat_snode[i].nzs / 1e6);
    fprintf(fp, "];\n");

    fclose(fp);

    SUPERLU_FREE (stat_snode);

}
#endif

int print_int_vec(char *what, int_t n, int_t *vec)
{
    int_t i;
    printf("%s\n", what);
    for (i = 0; i < n; ++i) printf(IFMT "\t" IFMT "\n", i, vec[i]);
    return 0;
}


/*
 * Print the parallel execution statistics.
 */
int_t ParallelProfile(const int_t n, const int_t supers, const int_t panels, 
		const int_t procs, Gstat_t *Gstat)
{
    register int_t i, imax, pruned, unpruned, waits, itemp, cs_numbers;
    register float loadmax, loadtot, temp, thresh, loadprint;
    register float waittime, cs_time;
    double    *utime = Gstat->utime;
    procstat_t *pstat;
    panstat_t *pan;
    void print_flops_by_height(int_t, panstat_t *, int_t *, float *);
    
    printf("\n---- Parallel Profile Per Processor ----\n");
    printf("%4s%16s%8s%10s%10s%10s%10s%8s\n", "proc", "factops",
	   "seconds", "skedwaits", "skedtime", "CS-time",
	   "spin-time", "[%tot]");
    for (i = 0; i < procs; ++i) {
	pstat = &(Gstat->procstat[i]);
	if ( pstat->fctime != 0 ) {
	    temp = pstat->spintime/pstat->fctime*100.;
	    printf( IFMT "%16e%8.2f" IFMT "%10.3f%10.3f%10.3f%8.1f\n", 
		   i, pstat->fcops, pstat->fctime, pstat->skedwaits,
		   pstat->skedtime, pstat->cs_time, pstat->spintime, temp);
	}
    }

    printf("%4s%8s%12s%14s\n",
	   "proc", "#panels", "dfs_pruned","dfs_unpruned");
    pruned = unpruned = 0;
    cs_time = 0.0;
    for (i = 0; i < procs; ++i) {
	pstat = &(Gstat->procstat[i]);
	printf(IFMT IFMT IFMT IFMT "\n", i, pstat->panels, pstat->pruned, pstat->unpruned);
	pruned += Gstat->procstat[i].pruned;
	unpruned += Gstat->procstat[i].unpruned;
	cs_time += Gstat->procstat[i].cs_time;
    }
    temp = pruned + unpruned;
    if ( temp != 0 ) {
    	printf("%12s%26s\n", "", "--------------------");
    	printf("%12s" IFMT IFMT "%14.0f\n", "total", pruned, unpruned, temp);
    	printf("%12s%12.2f%14.2f\n", "frac.", pruned/temp, unpruned/temp);
    }

    printf("%16s" IFMT "\n", "piped-panels", Gstat->panhows[PIPE]);
    printf("%16s" IFMT "\n", "nonpiped-DADs", Gstat->panhows[DADPAN]);
    printf("%16s" IFMT "\n", "nonpiped-panels", Gstat->panhows[NOPIPE]);

    /* work load distribution */
    loadmax = loadtot = Gstat->procstat[0].fcops;
    imax = 0;
    for (i = 1; i < procs; ++i) {
	temp = Gstat->procstat[i].fcops;
	loadtot += temp;
	if ( temp > loadmax ) {
	    loadmax = temp;
	    imax = i;
	}
    }
    printf("%25s%8.2f\n", "Load balance [mean/max]", loadtot/loadmax/procs);

    /* Delays due to pipelining. */
    waits = waittime = 0;
    for (i = 0; i < n; i += Gstat->panstat[i].size) { /* add up all panels */
	waits += Gstat->panstat[i].pipewaits;
	waittime += Gstat->panstat[i].spintime;
    }
    printf("%25s" IFMT ",\tper-panel %.1f\n", "total #delays in pipeline",
	    waits, (float)waits/panels);
    temp = waittime / procs;
    printf("%25s%8.2f\t [%.1f]\n", "mean spin time per-proc", 
	   temp, temp/utime[FACT]*100);
    
    /* Delays due to scheduling. */
    waits = waittime = 0;
    for (i = 0; i < procs; ++i) {
	waits += Gstat->procstat[i].skedwaits;
	waittime += Gstat->procstat[i].skedtime;
    }
    printf("%25s" IFMT "\n", "total #delays in schedule", waits);
    temp = waittime / procs;
    printf("%25s%8.2f\t [%.1f]\n", "mean sched. time per-proc", 
	   temp, temp/utime[FACT]*100);

    /* estimated overhead in spin-locks */
#if ( MACH==CRAY_PVP )    /* measured for mutex lock/unlock on 4 cpus */
#define TMUTEX          4.42e-6
#define FLOPS_PER_LOCK  221
#elif ( MACH==SUN )
#define TMUTEX          4.36e-6
#define FLOPS_PER_LOCK  109
#elif ( MACH==SGI || MACH==ORIGIN )
#define TMUTEX          2.02e-6
#define FLOPS_PER_LOCK  364
#elif ( MACH==DEC || PTHREAD )
#define TMUTEX          2.71e-6
#define FLOPS_PER_LOCK  407
#else
#define TMUTEX          2.00e-6
#define FLOPS_PER_LOCK  500
#endif
    cs_numbers = n + 3*supers + panels + procs; 
    itemp = cs_numbers * FLOPS_PER_LOCK;     /* translated to flops */
    temp = cs_numbers * TMUTEX;
    printf("mutex-lock overhead (est.) %8.2f, #locks " IFMT ", equiv. flops %e\n", 
	   temp, cs_numbers, (float) itemp);
    printf("time in critical section   %8.2f\t [%.1f]\n",
	   cs_time/procs, cs_time/procs/utime[FACT]*100);

    printf("\n---- Parallel Profile Per Panel ----\n");
    printf("%8s%8s%16s%8s%8s%12s%8s\n", "panel", "height",
	    "factops", "[tot%]", "msec", "spin(msec)", "Mflops");
    thresh = 0.005 * loadtot;
    loadprint = 0;
    itemp = 0;
    for (i = 0; i < n; i += Gstat->panstat[i].size) {
	pan = &(Gstat->panstat[i]);
	if ( pan->flopcnt > thresh ) {
	    loadprint += pan->flopcnt;
	    ++itemp;
	    if ( pan->fctime != 0 ) temp = pan->flopcnt/pan->fctime*1e-6;
	    printf(IFMT IFMT IFMT "%16e%8.1f%8.2f%12.2f%8.2f\n", i, pan->size,
		    Gstat->height[i], pan->flopcnt, pan->flopcnt/loadtot*100.,
		    pan->fctime*1e3, pan->spintime*1e3, temp);
	}
    }
    printf("Total panels " IFMT ",  height(T) " IFMT ", height(T)/n= %.4f\n", 
	   panels, Gstat->height[n], (float)Gstat->height[n]/n);
    printf("Printed flops %e [%.1f], printed panels " IFMT " [%.1f]\n",
	    loadprint, loadprint/loadtot*100.,
	    itemp, (float)itemp/panels);

/*    print_flops_by_height(n, panstat, height, flops_by_height);*/
	
    printf("---- End ParallelProfile().\n\n");
    fflush(stdout);
    return 0;
}

/*
 * Print the distribution of flops by the height of etree.
 */
void
print_flops_by_height(int_t n, panstat_t *panstat,
		      int_t *height, float *flops_by_height)
{
    register int_t i, w, ht;
    register float flops;

    for (i = 0; i < n; i += w) {
	w = panstat[i].size;
	ht = height[i];
	flops_by_height[ht] += panstat[i].flopcnt;
    }

    printf("\n%8s\t%8s\n", "height", "flops");
    ht = height[n-1]; /* root */
    for (i = 0; i <= ht; ++i) {
	flops = flops_by_height[i];
	if ( flops != 0.0 ) printf(IFMT "\t%e\n", i, flops);
    }
}

   
/*
 * Print the analysis of the optimal runtime.
 */
int_t
CPprofile(const int_t n, cp_panel_t *cp_panel, pxgstrf_shared_t *pxgstrf_shared)
{
    Gstat_t *Gstat = pxgstrf_shared->Gstat;
    register int_t maxpan, i, j, treecnt;
    register float eft, maxeft; /* earliest (possible) finish time */
    flops_t  *ops;

    /* Find the longest (weighted) path in the elimination forest. */
    treecnt = 0;
    maxeft = 0;
    for (i = Gstat->cp_firstkid[n]; i != EMPTY; i = Gstat->cp_nextkid[i]) {
/*	printf("Root " IFMT ", height " IFMT "\n", i, height[i]);*/
	j = (pxgstrf_shared->pan_status[i].size > 0) ? 
	  i : (i + pxgstrf_shared->pan_status[i].size);
	eft   = cp_panel[j].est + cp_panel[j].pdiv;
	if ( eft > maxeft ) {
	    maxeft = eft;
	    maxpan = j;
	}
	++treecnt;
    }
    
    ops   = Gstat->ops;
    printf("\n** Runtime prediction model: #trees " IFMT "\n", treecnt);
    printf("Last panel " IFMT ", seq-time %e, EFT %e, ideal-speedup %.2f\n",
	   maxpan, ops[FACT], maxeft, ops[FACT]/maxeft);

#if ( DEBUGlevel>=2 )
    printf("last panel " IFMT "\n", maxpan);
    for (i = 0; i < n; i += pxgstrf_shared->pan_status[i].size)
	printf("%6d %8s%e\t%8s%8.0f\n", i, "est  ", cp_panel[i].est,
	       "pdiv  ", cp_panel[i].pdiv);
#endif    
    return 0;
}


/***************************************************************
 * Utilities to print the supermatrix.
 ***************************************************************/

#define PERLINE  10
#define FMT      "%7.4f "

/* A is of type Stype==SCP */
void
Print_SuperNode_SCP(SuperMatrix *A)
{
    int_t i, j, c;
    int_t n = A->ncol;
    SCPformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int_t *colbeg = Astore->nzval_colbeg, *colend = Astore->nzval_colend;
    printf("SuperNode_SCP: nnz " IFMT ", nsuper " IFMT "\n", Astore->nnz, Astore->nsuper);
    printf("valL=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colbeg[j]; i < colend[j]; ++i) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	    ++c;
	}
    }
    printf("];\n");
    fflush(stdout);
    /*    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colend );*/
}

/* A is of type NCP */
void
Print_CompCol_NC(SuperMatrix *A)
{
    int_t i, j, c;
    int_t n = A->ncol;
    NCformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int_t *colptr = Astore->colptr;
    printf("CompCol_NC: nnz " IFMT "\n", Astore->nnz);
    printf("valA=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colptr[j]; i < colptr[j+1]; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	}
    }
    printf("];\n");
    fflush(stdout);
}

/* A is of type NCP */
void
Print_CompCol_NCP(SuperMatrix *A)
{
    int_t i, j, c;
    int_t n = A->ncol;
    NCPformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int_t *colbeg = Astore->colbeg, *colend = Astore->colend;
    printf("SuperNode_NCP: nnz " IFMT "\n", Astore->nnz);
    printf("nzval[U]\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colbeg[j]; i < colend[j]; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	}
    }
    printf("\n");
    fflush(stdout);
}

/* A is of type DN */
void
Print_Dense(SuperMatrix *A)
{
    int_t i, j, c;
    int_t m = A->nrow, n = A->ncol;
    DNformat *Astore = A->Store;
    int_t lda = Astore->lda;
    double *nzval = Astore->nzval;
    printf("Dense: lda " IFMT "\n", lda);
    printf("val=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i + j*lda]);
      }
    }
    printf("];\n");
    fflush(stdout);
}

