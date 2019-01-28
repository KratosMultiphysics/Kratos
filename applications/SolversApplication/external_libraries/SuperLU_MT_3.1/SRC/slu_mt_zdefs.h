/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 * April 20, 2015
 *
 * Sparse matrix types and function prototypes.
 *
 */

#ifndef __SLU_MT_ZDEFS /* allow multiple inclusions */
#define __SLU_MT_ZDEFS

/*
 * File name:           slu_mt_zdefs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

/****************************
  Include thread header file
  ***************************/
#if defined ( _SOLARIS )
#include <thread.h>
#include <sched.h>
#elif defined( _DEC )
#include <pthread.h>
#include <unistd.h>
#include <sys/mman.h>
#elif defined ( _OPENMP )
#include <omp.h>
#elif defined ( _PTHREAD )
#include <pthread.h>
#elif defined ( _CRAY )
#include <fortran.h>
#include <string.h>
#endif

/* Define my integer type int_t */
#ifdef _LONGINT
typedef long long int int_t;
#define IFMT "%lld"
#else
typedef int int_t; /* default */
#define IFMT "%8d"
#endif

#include "slu_mt_machines.h"
#include "slu_mt_Cnames.h"
#include "supermatrix.h"
#include "slu_mt_util.h"
#include "pxgstrf_synch.h"

#include "slu_dcomplex.h"

/*
 * *************************************************
 *  Global data structures used in LU factorization
 * *************************************************
 * 
 *   nsuper: number of supernodes = nsuper+1, numbered between 0 and nsuper.
 *
 *   (supno, xsup, xsup_end):
 *      supno[i] is the supernode number to which column i belongs;
 *	xsup[s] points to the first column of supernode s;
 *      xsup_end[s] points to one past the last column of supernode s.
 *	Example: supno  0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	          xsup  0 1 2 4 7
 *            xsup_end  1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (lsub, xlsub, xlsub_end):
 *      lsub[*] contains the compressed subscripts of the supernodes;
 *      xlsub[j] points to the starting location of the j-th column in
 *               lsub[*]; 
 *      xlsub_end[j] points to one past the ending location of the j-th
 *               column in lsub[*].
 *	Storage: original row subscripts in A.
 *
 *      During the course of sparse LU factorization, we also use
 *	(lsub, xlsub, xlsub_end, xprune) to represent symmetrically 
 *      pruned graph. Contention will occur when one processor is
 *      performing DFS on supernode S, while another processor is pruning
 *      supernode S. We use the following data structure to deal with
 *      this problem. Suppose each supernode contains columns {s,s+1,...,t},
 *      with first column s and last column t.
 *
 *      (1) if t > s, only the subscript sets for column s and column t
 *          are stored. Column t represents pruned adjacency structure.
 *
 *                  --------------------------------------------
 *          lsub[*]    ... |   col s    |   col t   | ...
 *                  --------------------------------------------
 *                          ^            ^           ^
 *                       xlsub[s]    xlsub_end[s]  xlsub_end[s+1]
 *                                   xlsub[s+1]      :
 *                                       :           :
 *                                       :         xlsub_end[t]
 *                                   xlsub[t]      xprune[t] 
 *                                   xprune[s]    
 *
 *      (2) if t == s, i.e., a singleton supernode, the subscript set
 *          is stored twice:
 *
 *                  --------------------------------------
 *          lsub[*]    ... |      s     |     s     | ...
 *                  --------------------------------------
 *                          ^            ^           ^
 *                       xlsub[s]   xlsub_end[s]  xprune[s]
 *
 *      There are two subscript sets for each supernode, the last column
 *      structures (for pruning) will be removed after the numerical LU
 *      factorization phase:
 *        o lsub[j], j = xlsub[s], ..., xlsub_end[s]-1
 *          is the structure of column s (i.e. structure of this supernode).
 *          It is used for the storage of numerical values.
 *	  o lsub[j], j = xlsub[t], ..., xlsub_end[t]-1
 *	    is the structure of the last column t of this supernode.
 *	    It is for the purpose of symmetric pruning. Therefore, the
 *	    structural subscripts can be rearranged without making physical
 *	    interchanges among the numerical values.
 *
 *       DFS will traverse the first subscript set if the supernode
 *       has not been pruned; otherwise it will traverse the second
 *       subscript list, i.e., the part of the pruned graph.
 *
 *   (lusup, xlusup, xlusup_end):
 *      lusup[*] contains the numerical values of the supernodes;
 *      xlusup[j] points to the starting location of the j-th column in
 *                storage vector lusup[*]; 
 *      xlusup_end[j] points to one past the ending location of the j-th 
 *                column in lusup[*].
 *	Each supernode is stored in column-major, consistent with Fortran
 *      two-dimensional array storage.
 *
 *   (ucol, usub, xusub, xusub_end):
 *      ucol[*] stores the numerical values of the U-columns above the
 *              supernodes. 
 *      usub[k] stores the row subscripts of nonzeros ucol[k];
 *      xusub[j] points to the starting location of column j in ucol/usub[]; 
 *      xusub_end[j] points to one past the ending location column j in
 *                   ucol/usub[].
 *	Storage: new row subscripts; that is indexed intp PA.
 *
 */
typedef struct {
    int_t     *xsup;    /* supernode and column mapping */
    int_t     *xsup_end;
    int_t     *supno;   
    int_t     *lsub;    /* compressed L subscripts */
    int_t	    *xlsub;
    int_t     *xlsub_end;
    doublecomplex  *lusup;   /* L supernodes */
    int_t     *xlusup;
    int_t     *xlusup_end;
    doublecomplex  *ucol;    /* U columns */
    int_t     *usub;
    int_t	    *xusub;
    int_t     *xusub_end;
    int_t     nsuper;   /* current supernode number */
    int_t     nextl;    /* next position in lsub[] */
    int_t     nextu;    /* next position in usub[]/ucol[] */
    int_t     nextlu;   /* next position in lusup[] */
    int_t     nzlmax;   /* current max size of lsub[] */
    int_t     nzumax;   /*    "    "    "      ucol[] */
    int_t     nzlumax;  /*    "    "    "     lusup[] */
    /* ---------------------------------------------------------------
     *  Memory managemant for L supernodes 
     */
    int_t  *map_in_sup;  /* size n+1 - the address offset of each column
                        * in lusup[*], which is divided into regions 
			* by the supernodes of Householder matrix H.
			* If column k starts a supernode in H,
			* map_in_sup[k] is the next open position in
			* lusup[*]; otherwise map_in_sup[k] gives the
			* offset (negative) to the leading column
			* of the supernode in H.
			*/
    int_t  dynamic_snode_bound;
    /* --------------------------------------------------------------- */
} GlobalLU_t;


/* 
 * *********************************************************************
 * The pxgstrf_shared_t structure contains the shared task queue and
 * the synchronization variables to facilitate parallel factorization. 
 * It also contains the shared L and U data structures.
 * *********************************************************************
 */
typedef struct {
    /* ----------------------------------------------------------------
     * Global variables introduced in parallel code for synchronization.
     */
    volatile int_t tasks_remain; /* number of untaken panels */
    int_t          num_splits;   /* number of panels split at the top */
    queue_t      taskq;        /* size ncol - shared work queue */
    mutex_t      *lu_locks;    /* 5 named mutual exclusive locks */
    volatile int_t *spin_locks;  /* size ncol - mark every busy column */
    pan_status_t *pan_status;  /* size ncol - panel status */
    int_t          *fb_cols;     /* size ncol - mark farthest busy column */
    /* ---------------------------------------------------------------- */
    int_t        *inv_perm_c;
    int_t        *inv_perm_r;
    int_t        *xprune;
    int_t        *ispruned;
    SuperMatrix *A;
    GlobalLU_t *Glu;
    Gstat_t    *Gstat;
    int_t        *info;
} pxgstrf_shared_t;

/* Arguments passed to each thread. */
typedef struct {
    int_t  pnum; /* process number */
    int_t  info; /* error code returned from each thread */       
    superlumt_options_t *superlumt_options;
    pxgstrf_shared_t  *pxgstrf_shared; /* shared for LU factorization */
} pzgstrf_threadarg_t;


/* *********************
   Function prototypes
   *********************/

#ifdef __cplusplus
extern "C" {
#endif


/* ----------------
   Driver routines 
   ---------------*/
extern void
pzgssv(int_t, SuperMatrix *, int_t *, int_t *, SuperMatrix *, SuperMatrix *, 
       SuperMatrix *, int_t *);
extern void
pzgssvx(int_t, superlumt_options_t *, SuperMatrix *, int_t *, int_t *,  
	equed_t *, double *, double *, SuperMatrix *, SuperMatrix *,
	SuperMatrix *, SuperMatrix *, 
	double *, double *, double *, double *, superlu_memusage_t *, 
	int_t *);

/* ---------------
   Driver related 
   ---------------*/
extern void zgsequ (SuperMatrix *, double *, double *, double *,
                    double *, double *, int_t *);
extern void zlaqgs (SuperMatrix *, double *, double *, double, 
		    double, double, equed_t *);
extern void zgscon (char *, SuperMatrix *, SuperMatrix *,
		    double, double *, int_t *);
extern double zPivotGrowth(int_t, SuperMatrix *, int_t *,
			   SuperMatrix *, SuperMatrix *);
extern void zgsrfs (trans_t, SuperMatrix *, SuperMatrix *, SuperMatrix *,
		    int_t *, int_t *, equed_t, double *, double *, SuperMatrix *,
		    SuperMatrix *, double *, double *, Gstat_t *, int_t *);
extern int_t  sp_ztrsv (char *, char *, char *, SuperMatrix *, SuperMatrix *,
		      doublecomplex *, int_t *);
extern int_t  sp_zgemv (char *, doublecomplex, SuperMatrix *, doublecomplex *,
		      int_t, doublecomplex, doublecomplex *, int_t);
extern int_t  sp_zgemm (char *, int_t, int_t, int_t, doublecomplex, SuperMatrix *, 
		      doublecomplex *, int_t, doublecomplex, doublecomplex *, int_t);

/* ----------------------
   Factorization related
   ----------------------*/
extern void pxgstrf_scheduler (const int_t, const int_t, const int_t *,
			       int_t *, int_t *, pxgstrf_shared_t *);
extern int_t  zParallelInit (int_t, pxgstrf_relax_t *, superlumt_options_t *,
			  pxgstrf_shared_t *);
extern int_t  ParallelFinalize ();
extern void pzgstrf_StackFree ();
extern int_t  queue_init (queue_t *, int_t);
extern int_t  queue_destroy (queue_t *);
extern int_t  EnqueueRelaxSnode (queue_t *, int_t, pxgstrf_relax_t *,
			      pxgstrf_shared_t *);
extern int_t  EnqueueDomains(queue_t *, struct Branch *, pxgstrf_shared_t *);
extern int_t  Enqueue (queue_t *, qitem_t);
extern int_t  Dequeue (queue_t *, qitem_t *);
extern int_t  NewNsuper (const int_t, pxgstrf_shared_t *, int_t *);
extern int_t  lockon(int_t *);
extern void PartDomains(const int_t, const float, SuperMatrix *, int_t *, int_t *);

extern void
zCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, doublecomplex *,
		      int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
void
zCreate_CompCol_Permuted(SuperMatrix *, int_t, int_t, int_t, doublecomplex *, int_t *,
			 int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, doublecomplex *, int_t,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, int_t, int_t, int_t, doublecomplex *, int_t *, int_t *,
			int_t *, int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Permuted(SuperMatrix *, int_t, int_t, int_t, doublecomplex *, 
			   int_t *, int_t *, int_t *, int_t *, int_t *, int_t *, 
			   int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_Dense_Matrix(int_t, int_t, doublecomplex *, int_t, doublecomplex *, int_t);

extern void Destroy_SuperMatrix_Store(SuperMatrix *);
extern void Destroy_CompCol_Matrix(SuperMatrix *);
extern void Destroy_CompCol_Permuted(SuperMatrix *);
extern void Destroy_CompCol_NCP(SuperMatrix *);
extern void Destroy_SuperNode_Matrix(SuperMatrix *);
extern void Destroy_SuperNode_SCP(SuperMatrix *);

extern void zallocateA (int_t, int_t, doublecomplex **, int_t **, int_t **);
extern void StatAlloc (const int_t, const int_t, const int_t, const int_t, Gstat_t*);
extern void StatInit  (const int_t, const int_t, Gstat_t*);
extern void StatFree  (Gstat_t*);
extern void get_perm_c(int_t, SuperMatrix *, int_t *);
extern void zsp_colorder (SuperMatrix *, int_t *, superlumt_options_t *,
			 SuperMatrix *);
extern int_t  sp_coletree (int_t *, int_t *, int_t *, int_t, int_t, int_t *);
extern int_t  zPresetMap (const int_t, SuperMatrix *, pxgstrf_relax_t *, 
		       superlumt_options_t *, GlobalLU_t *);
extern int_t  qrnzcnt (int_t, int_t, int_t *, int_t *, int_t *, int_t *, int_t *, int_t *,
		     int_t *, int_t *, int_t *, int_t *);
extern int_t  DynamicSetMap(const int_t, const int_t, const int_t, pxgstrf_shared_t*);
extern void pzgstrf (superlumt_options_t *, SuperMatrix *, int_t *, 
		     SuperMatrix *, SuperMatrix *, Gstat_t *, int_t *);
extern void pzgstrf_init (int_t, fact_t, trans_t, yes_no_t, int_t, int_t, double, yes_no_t, double,
			  int_t *, int_t *, void *, int_t, SuperMatrix *,
			  SuperMatrix *, superlumt_options_t *, Gstat_t *);
extern pzgstrf_threadarg_t*
pzgstrf_thread_init (SuperMatrix *, SuperMatrix *, SuperMatrix *,
		     superlumt_options_t*, pxgstrf_shared_t*, Gstat_t*, int_t*);
extern void
pzgstrf_thread_finalize (pzgstrf_threadarg_t *, pxgstrf_shared_t *,
			 SuperMatrix *, int_t *, SuperMatrix *, SuperMatrix *);
extern void pzgstrf_finalize(superlumt_options_t *, SuperMatrix *);
extern void pxgstrf_finalize(superlumt_options_t *, SuperMatrix *);
extern void pzgstrf_relax_snode (const int_t, superlumt_options_t *,
				 pxgstrf_relax_t *);
extern int_t
pzgstrf_factor_snode (const int_t, const int_t, SuperMatrix *, const double,
		      yes_no_t *, int_t *, int_t *, int_t*, int_t*, int_t*, int_t*,
		      doublecomplex *, doublecomplex *, pxgstrf_shared_t *, int_t *);
extern void
pxgstrf_mark_busy_descends (int_t, int_t, int_t *, pxgstrf_shared_t *, int_t *, int_t *);
extern int_t  pzgstrf_snode_dfs (const int_t, const int_t, const int_t, const int_t *,
			       const int_t *, const int_t *, int_t*, int_t *, int_t *,
			       pxgstrf_shared_t *);
extern int_t  pzgstrf_snode_bmod (const int_t, const int_t, const int_t, const int_t,
				doublecomplex *, doublecomplex *, GlobalLU_t*, Gstat_t*);
extern void pzgstrf_panel_dfs (const int_t, const int_t, const int_t, const int_t,
			       SuperMatrix *, int_t*, int_t*, int_t*, int_t*, int_t*, 
			       int_t*, int_t*, int_t*, int_t*, int_t*, int_t*, int_t*, int_t*,
			       doublecomplex*, GlobalLU_t *);
extern void pzgstrf_panel_bmod (const int_t, const int_t, const int_t, const int_t,
				const int_t, int_t*, int_t*, int_t*, int_t*, int_t*, int_t*,
				int_t*, int_t*, doublecomplex*, doublecomplex*, 
				pxgstrf_shared_t *);
extern void pzgstrf_bmod1D (const int_t, const int_t, const int_t, const int_t, 
			    const int_t, const int_t, const int_t, int_t, int_t,
			    int_t *, int_t *, int_t *, int_t *, doublecomplex *, doublecomplex *, 
			    GlobalLU_t *, Gstat_t *);
extern void pzgstrf_bmod2D (const int_t, const int_t, const int_t, const int_t,
			    const int_t, const int_t, const int_t, int_t, int_t,
			    int_t *, int_t *, int_t *, int_t *, doublecomplex *, doublecomplex *,
			    GlobalLU_t *, Gstat_t *);
extern void pzgstrf_bmod1D_mv2 (const int_t, const int_t, const int_t, const int_t, 
				const int_t, const int_t, const int_t, int_t, int_t,
				int_t *, int_t *, int_t *, int_t *, doublecomplex *, 
				doublecomplex *, GlobalLU_t *, Gstat_t *);
extern void pzgstrf_bmod2D_mv2 (const int_t, const int_t, const int_t, const int_t,
				const int_t, const int_t, const int_t, int_t, int_t,
				int_t *, int_t *, int_t *, int_t *, doublecomplex *, doublecomplex *,
				GlobalLU_t *, Gstat_t *);
extern void pxgstrf_super_bnd_dfs (const int_t, const int_t, const int_t, 
				   const int_t, const int_t, SuperMatrix*,
				   int_t*, int_t*, int_t*, int_t *, int_t *, int_t *,
				   int_t *, pxgstrf_shared_t *);
extern int_t  pzgstrf_column_dfs(const int_t, const int_t, const int_t, const int_t,
			       int_t*, int_t*, int_t*, int_t, int_t*, int_t*, int_t*, int_t*,
			       int_t *, int_t *, int_t *, int_t *, pxgstrf_shared_t *);
extern int_t  pzgstrf_column_bmod(const int_t, const int_t, const int_t, const int_t, 
				int_t*, int_t*, doublecomplex*, doublecomplex*,
				pxgstrf_shared_t *, Gstat_t *);
extern int_t  pzgstrf_pivotL (const int_t, const int_t, const double, yes_no_t*,
			    int_t*, int_t*, int_t*, int_t*, GlobalLU_t*, Gstat_t*);
extern int_t  pzgstrf_copy_to_ucol (const int_t, const int_t, const int_t, const int_t *,
				  const int_t *, const int_t *, doublecomplex*,
				  pxgstrf_shared_t*);
extern void pxgstrf_pruneL (const int_t, const int_t *, const int_t, const int_t,
			    const int_t *, const int_t *, int_t*, int_t *,
			    GlobalLU_t *);
extern void pxgstrf_resetrep_col (const int_t, const int_t *, int_t *);
extern void countnz (const int_t, int_t*, int_t *, int_t *, GlobalLU_t *);
extern void fixupL (const int_t, const int_t *, GlobalLU_t *);
extern void compressSUP (const int_t, GlobalLU_t *);
extern int_t  spcoletree (int_t *, int_t *, int_t *, int_t, int_t, int_t *);
extern int_t  *TreePostorder (int_t, int_t *);
extern void zreadmt (int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void zreadhb (int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void zGenXtrue (int_t, int_t, doublecomplex *, int_t);
extern void zFillRHS (trans_t, int_t, doublecomplex *, int_t, 
		      SuperMatrix *, SuperMatrix *);
extern void zgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
		    int_t*, int_t*, SuperMatrix*, Gstat_t *, int_t *);
extern void zlsolve (int_t, int_t, doublecomplex *, doublecomplex *);
extern void zusolve (int_t, int_t, doublecomplex *, doublecomplex *);
extern void zmatvec (int_t, int_t, int_t, doublecomplex *, doublecomplex *, doublecomplex *);


/* ---------------
   BLAS 
   ---------------*/
extern int zgemm_(char*, char*, int*, int*, int*, doublecomplex*,
                  doublecomplex*, int*, doublecomplex*, int*, doublecomplex*,
                  doublecomplex*, int*);
extern int ztrsm_(char*, char*, char*, char*, int*, int*, doublecomplex*,
                  doublecomplex*, int*, doublecomplex*, int*);
extern int ztrsv_(char*, char*, char*, int*, doublecomplex*, int*,
                  doublecomplex*, int*);
extern int zgemv_(char*, int*, int*, doublecomplex*, doublecomplex*, 
		   int*, doublecomplex*, int*, doublecomplex*, doublecomplex*, int*);

/* ---------------
   Memory related 
   ---------------*/
extern float pzgstrf_MemInit (int_t, int_t, superlumt_options_t *,
			SuperMatrix *, SuperMatrix *, GlobalLU_t *);
extern float pzgstrf_memory_use(const int_t, const int_t, const int_t);
extern int_t  pzgstrf_WorkInit (int_t, int_t, int_t **, doublecomplex **);
extern void pxgstrf_SetIWork (int_t, int_t, int_t *, int_t **, int_t **, int_t **,
		      int_t **, int_t **, int_t **, int_t **);
extern void pzgstrf_SetRWork (int_t, int_t, doublecomplex *, doublecomplex **, doublecomplex **);
extern void pzgstrf_WorkFree (int_t *, doublecomplex *, GlobalLU_t *);
extern int_t  pzgstrf_MemXpand (int_t, int_t, MemType, int_t *, GlobalLU_t *);

extern int_t  *intMalloc (int_t);
extern int_t  *intCalloc (int_t);
extern doublecomplex *doublecomplexMalloc(int_t);
extern doublecomplex *doublecomplexCalloc(int_t);
extern int_t  memory_usage ();
extern int_t  superlu_zQuerySpace (int_t, SuperMatrix *, SuperMatrix *, int_t, 
				 superlu_memusage_t *);
extern int_t  Glu_alloc (const int_t, const int_t, const int_t, const MemType,
		       int_t *, pxgstrf_shared_t *);

/* -------------------
   Auxiliary routines
   -------------------*/
extern double  SuperLU_timer_();
extern int_t     sp_ienv(int_t);
extern double  dlamch_();
extern int     lsame_(char *, char *);
extern int     xerbla_(char *, int *);
extern void    superlu_abort_and_exit(char *);
extern void    ifill(int_t *, int_t, int_t);
extern void    zfill(doublecomplex *, int_t, doublecomplex);
extern void    zinf_norm_error(int_t, SuperMatrix *, doublecomplex *);
extern void    dstat_allocate(int_t);
extern void    snode_profile(int_t, int_t *);
extern void    super_stats(int_t, int_t *, int_t *);
extern void    panel_stats(int_t, int_t, int_t *, Gstat_t *);
extern void    PrintSumm(char *, int_t, int_t, int_t);
extern void    zPrintPerf(SuperMatrix *, SuperMatrix *, superlu_memusage_t *,
			 double, double, double *, double *, char *,
			 Gstat_t *);
extern void    zCompRow_to_CompCol(int_t m, int_t n, int_t nnz, 
                           doublecomplex *a, int_t *colind, int_t *rowptr,
                           doublecomplex **at, int_t **rowind, int_t **colptr);


/* -----------------------
   Routines for debugging
   -----------------------*/
extern void    print_lu_col(int_t, char *, int_t, int_t, int_t, int_t *, GlobalLU_t *);
extern void    print_panel_seg(int_t, int_t, int_t, int_t, int_t *, int_t *);
extern void    zcheck_zero_vec(int_t, char *, int_t, doublecomplex *);
extern void    check_repfnz(int_t, int_t, int_t, int_t *);

#ifdef __cplusplus
	   }
#endif


#endif /* __SLU_MT_ZDEFS */

