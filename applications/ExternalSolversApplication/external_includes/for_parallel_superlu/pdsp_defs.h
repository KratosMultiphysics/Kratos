
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Sparse matrix types and function prototypes.
 *
 */

#ifndef __SUPERLU_dSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_dSP_DEFS

/*
 * File name:           pdsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

/****************************
  Include thread header file
  ***************************/
// #if defined ( _SOLARIS )
// #include <thread.h>
// #include <sched.h>
// #elif defined( _DEC )
// #include <pthread.h>
// #include <unistd.h>
// #include <sys/mman.h>
// #elif defined ( _OPENMP )
#include"/opt/intel/cce/10.1.015/include/omp.h"
//#include <omp.h>
// #elif defined ( _PTHREAD )
// #include <pthread.h>
// #elif defined ( _CRAY )
// #include <fortran.h>
// #include <string.h>
// #endif

#if defined ( MAX )
#undef	MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#if defined ( MIN )
#undef	MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#if defined ( ABS )
#undef	ABS
#define ABS(a)	   (((a) < 0) ? -(a) : (a))
#endif

/* Define my integer type int_t */
typedef int int_t; /* default */

#include "external_includes/superlu/slu_mt_machines.h"
#include "external_includes/superlu/slu_mt_Cnames.h"
#include "external_includes/superlu/supermatrix.h"
#include "external_includes/superlu/slu_mt_util.h"
#include "external_includes/superlu/pxgstrf_synch.h"


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
    int     *xsup;    /* supernode and column mapping */
    int     *xsup_end;
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    int     *xlsub_end;
    double  *lusup;   /* L supernodes */
    int     *xlusup;
    int     *xlusup_end;
    double  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     *xusub_end;
    int     nsuper;   /* current supernode number */
    int     nextl;    /* next position in lsub[] */
    int     nextu;    /* next position in usub[]/ucol[] */
    int     nextlu;   /* next position in lusup[] */
    int     nzlmax;   /* current max size of lsub[] */
    int     nzumax;   /*    "    "    "      ucol[] */
    int     nzlumax;  /*    "    "    "     lusup[] */
    /* ---------------------------------------------------------------
     *  Memory managemant for L supernodes 
     */
    int  *map_in_sup;  /* size n+1 - the address offset of each column
                        * in lusup[*], which is divided into regions 
			* by the supernodes of Householder matrix H.
			* If column k starts a supernode in H,
			* map_in_sup[k] is the next open position in
			* lusup[*]; otherwise map_in_sup[k] gives the
			* offset (negative) to the leading column
			* of the supernode in H.
			*/
    int  dynamic_snode_bound;
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
    volatile int tasks_remain; /* number of untaken panels */
    int          num_splits;   /* number of panels split at the top */
    queue_t      taskq;        /* size ncol - shared work queue */
    mutex_t      *lu_locks;    /* 5 named mutual exclusive locks */
    volatile int *spin_locks;  /* size ncol - mark every busy column */
    pan_status_t *pan_status;  /* size ncol - panel status */
    int          *fb_cols;     /* size ncol - mark farthest busy column */
    /* ---------------------------------------------------------------- */
    int        *inv_perm_c;
    int        *inv_perm_r;
    int        *xprune;
    int        *ispruned;
    SuperMatrix *A;
    GlobalLU_t *Glu;
    Gstat_t    *Gstat;
    int        *info;
} pxgstrf_shared_t;

/* Arguments passed to each thread. */
typedef struct {
    int  pnum; /* process number */
    int  info; /* error code returned from each thread */       
    superlumt_options_t *superlumt_options;
    pxgstrf_shared_t  *pxgstrf_shared; /* shared for LU factorization */
} pdgstrf_threadarg_t;

/*-----------------*/
#if 0 

/* The structure to record a relaxed supernode. */
typedef struct {
    int fcol;
    int size;
} pxgstrf_relax_t;

/* Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    int size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

/* The structure to keep track of memory usage. */
typedef struct {
    float for_lu;
    float total_needed;
    int   expansions;
} superlu_memusage_t;


/* *******
   Macros
   *******/

#define SUPER_REP(s)    ( xsup_end[s]-1 )
#define SUPER_FSUPC(s)  ( xsup[s] )
#define SINGLETON(s)    ( (xsup_end[s] - xsup[s]) == 1 )
#define ISPRUNED(j)     ( ispruned[j] )
#define STATE(j)        ( pxgstrf_shared->pan_status[j].state )
#define DADPANEL(j)     ( etree[j + pxgstrf_shared->pan_status[j].size-1] )

#ifdef PROFILE
#define TIC(t)          t = SuperLU_timer_()
#define TOC(t2, t1)     t2 = SuperLU_timer_() - t1
#else
#define TIC(t)
#define TOC(t2, t1)
#endif

/*-----------------*/
#endif


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
pdgssv(int, SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
       SuperMatrix *, int *);
extern void
pdgssvx(int, superlumt_options_t *, SuperMatrix *, int *, int *,  
	equed_t *, double *, double *, SuperMatrix *, SuperMatrix *,
	SuperMatrix *, SuperMatrix *, 
	double *, double *, double *, double *, superlu_memusage_t *, 
	int *);

/* ---------------
   Driver related 
   ---------------*/
extern void dgsequ (SuperMatrix *, double *, double *, double *,
                    double *, double *, int *);
extern void dlaqgs (SuperMatrix *, double *, double *, double, 
		    double, double, equed_t *);
extern void dgscon (char *, SuperMatrix *, SuperMatrix *,
		    double, double *, int *);
extern double dPivotGrowth(int, SuperMatrix *, int *,
			   SuperMatrix *, SuperMatrix *);
extern void dgsrfs (trans_t, SuperMatrix *, SuperMatrix *, SuperMatrix *,
		    int *, int *, equed_t, double *, double *, SuperMatrix *,
		    SuperMatrix *, double *, double *, Gstat_t *, int *);
extern int  sp_dtrsv (char *, char *, char *, SuperMatrix *, SuperMatrix *,
		      double *, int *);
extern int  sp_dgemv (char *, double, SuperMatrix *, double *,
		      int, double, double *, int);
extern int  sp_dgemm (char *, int, int, int, double, SuperMatrix *, 
		      double *, int, double, double *, int);

/* ----------------------
   Factorization related
   ----------------------*/
extern void pxgstrf_scheduler (const int, const int, const int *,
			       int *, int *, pxgstrf_shared_t *);
extern int  dParallelInit (int, pxgstrf_relax_t *, superlumt_options_t *,
			  pxgstrf_shared_t *);
extern int  ParallelFinalize ();
extern int  queue_init (queue_t *, int);
extern int  queue_destroy (queue_t *);
extern int  EnqueueRelaxSnode (queue_t *, int, pxgstrf_relax_t *,
			      pxgstrf_shared_t *);
extern int  EnqueueDomains(queue_t *, struct Branch *, pxgstrf_shared_t *);
extern int  Enqueue (queue_t *, qitem_t);
extern int  Dequeue (queue_t *, qitem_t *);
extern int  NewNsuper (const int, mutex_t *, int *);
extern int  lockon(int *);
extern void PartDomains(const int, const float, SuperMatrix *, int *, int *);

extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		      int *, int *, Stype_t, Dtype_t, Mtype_t);
void
dCreate_CompCol_Permuted(SuperMatrix *, int, int, int, double *, int *,
			 int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, double *, int *, int *,
			int *, int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Permuted(SuperMatrix *, int, int, int, double *, 
			   int *, int *, int *, int *, int *, int *, 
			   int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix(int, int, double *, int, double *, int);

extern void Destroy_SuperMatrix_Store(SuperMatrix *);
extern void Destroy_CompCol_Matrix(SuperMatrix *);
extern void Destroy_CompCol_Permuted(SuperMatrix *);
extern void Destroy_CompCol_NCP(SuperMatrix *);
extern void Destroy_SuperNode_Matrix(SuperMatrix *);
extern void Destroy_SuperNode_SCP(SuperMatrix *);

extern void dallocateA (int, int, double **, int **, int **);
extern void StatAlloc (const int, const int, const int, const int, Gstat_t*);
extern void StatInit  (const int, const int, Gstat_t*);
extern void StatFree  (Gstat_t*);
extern void get_perm_c(int, SuperMatrix *, int *);
extern void dsp_colorder (SuperMatrix *, int *, superlumt_options_t *,
			 SuperMatrix *);
extern int  sp_coletree (int *, int *, int *, int, int, int *);
extern int  dPresetMap (const int, SuperMatrix *, pxgstrf_relax_t *, 
		       superlumt_options_t *, GlobalLU_t *);
extern int  qrnzcnt (int, int, int *, int *, int *, int *, int *, int *,
		     int *, int *, int *, int *);
extern int  DynamicSetMap(const int, const int, const int, pxgstrf_shared_t*);
extern void pdgstrf (superlumt_options_t *, SuperMatrix *, int *, 
		     SuperMatrix *, SuperMatrix *, Gstat_t *, int *);
extern void pdgstrf_init (int, fact_t, trans_t, yes_no_t, int, int, double, yes_no_t, double,
			  int *, int *, void *, int, SuperMatrix *,
			  SuperMatrix *, superlumt_options_t *, Gstat_t *);
extern pdgstrf_threadarg_t*
pdgstrf_thread_init (SuperMatrix *, SuperMatrix *, SuperMatrix *,
		     superlumt_options_t*, pxgstrf_shared_t*, Gstat_t*, int*);
extern void
pdgstrf_thread_finalize (pdgstrf_threadarg_t *, pxgstrf_shared_t *,
			 SuperMatrix *, int *, SuperMatrix *, SuperMatrix *);
extern void pdgstrf_finalize(superlumt_options_t *, SuperMatrix *);
extern void pdgstrf_relax_snode (const int, superlumt_options_t *,
				 pxgstrf_relax_t *);
extern int
pdgstrf_factor_snode (const int, const int, SuperMatrix *, const double,
		      yes_no_t *, int *, int *, int*, int*, int*, int*,
		      double *, double *, pxgstrf_shared_t *, int *);
extern void
pxgstrf_mark_busy_descends (int, int, int *, pxgstrf_shared_t *, int *, int *);
extern int  pdgstrf_snode_dfs (const int, const int, const int, const int *,
			       const int *, const int *, int*, int *, int *,
			       pxgstrf_shared_t *);
extern int  pdgstrf_snode_bmod (const int, const int, const int, const int,
				double *, double *, GlobalLU_t*, Gstat_t*);
extern void pdgstrf_panel_dfs (const int, const int, const int, const int,
			       SuperMatrix *, int*, int*, int*, int*, int*, 
			       int*, int*, int*, int*, int*, int*, int*, int*,
			       double*, GlobalLU_t *);
extern void pdgstrf_panel_bmod (const int, const int, const int, const int,
				const int, int*, int*, int*, int*, int*, int*,
				int*, int*, double*, double*, 
				pxgstrf_shared_t *);
extern void pdgstrf_bmod1D (const int, const int, const int, const int, 
			    const int, const int, const int, int, int,
			    int *, int *, int *, int *, double *, double *, 
			    GlobalLU_t *, Gstat_t *);
extern void pdgstrf_bmod2D (const int, const int, const int, const int,
			    const int, const int, const int, int, int,
			    int *, int *, int *, int *, double *, double *,
			    GlobalLU_t *, Gstat_t *);
extern void pdgstrf_bmod1D_mv2 (const int, const int, const int, const int, 
				const int, const int, const int, int, int,
				int *, int *, int *, int *, double *, 
				double *, GlobalLU_t *, Gstat_t *);
extern void pdgstrf_bmod2D_mv2 (const int, const int, const int, const int,
				const int, const int, const int, int, int,
				int *, int *, int *, int *, double *, double *,
				GlobalLU_t *, Gstat_t *);
extern void pxgstrf_super_bnd_dfs (const int, const int, const int, 
				   const int, const int, SuperMatrix*,
				   int*, int*, int*, int *, int *, int *,
				   int *, pxgstrf_shared_t *);
extern int  pdgstrf_column_dfs(const int, const int, const int, const int,
			       int*, int*, int*, int, int*, int*, int*, int*,
			       int *, int *, int *, int *, pxgstrf_shared_t *);
extern int  pdgstrf_column_bmod(const int, const int, const int, const int, 
				int*, int*, double*, double*,
				pxgstrf_shared_t *, Gstat_t *);
extern int  pdgstrf_pivotL (const int, const int, const double, yes_no_t*,
			    int*, int*, int*, int*, GlobalLU_t*, Gstat_t*);
extern int  pdgstrf_copy_to_ucol (const int, const int, const int, const int *,
				  const int *, const int *, double*,
				  pxgstrf_shared_t*);
extern void pxgstrf_pruneL (const int, const int *, const int, const int,
			    const int *, const int *, int*, int *,
			    GlobalLU_t *);
extern void pxgstrf_resetrep_col (const int, const int *, int *);
extern void countnz (const int, int*, int *, int *, GlobalLU_t *);
extern void fixupL (const int, const int *, GlobalLU_t *);
extern void compressSUP (const int, GlobalLU_t *);
extern int  spcoletree (int *, int *, int *, int, int, int *);
extern int  *TreePostorder (int, int *);
extern void dreadmt (int *, int *, int *, double **, int **, int **);
extern void dreadhb (int *, int *, int *, double **, int **, int **);
extern void dGenXtrue (int, int, double *, int);
extern void dFillRHS (trans_t, int, double *, int, 
		      SuperMatrix *, SuperMatrix *);
extern void dgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
		    int*, int*, SuperMatrix*, Gstat_t *, int *);
extern void dlsolve (int, int, double *, double *);
extern void dusolve (int, int, double *, double *);
extern void dmatvec (int, int, int, double *, double *, double *);


/* ---------------
   BLAS 
   ---------------*/
extern int dgemm_(char*, char*, int*, int*, int*, double*,
                  double*, int*, double*, int*, double*,
                  double*, int*);
extern int dtrsm_(char*, char*, char*, char*, int*, int*, double*,
                  double*, int*, double*, int*);
extern int dtrsv_(char*, char*, char*, int*, double*, int*,
                  double*, int*);
extern int dgemv_(char*, int*, int*, double*, double*, 
		   int*, double*, int*, double*, double*, int*);

/* ---------------
   Memory related 
   ---------------*/
extern int  pdgstrf_MemInit (int, int, superlumt_options_t *,
			SuperMatrix *, SuperMatrix *, GlobalLU_t *);
extern int  pdgstrf_memory_use(const int, const int, const int);
extern int  pdgstrf_WorkInit (int, int, int **, double **);
extern void pxgstrf_SetIWork (int, int, int *, int **, int **, int **,
		      int **, int **, int **, int **);
extern void pdgstrf_SetRWork (int, int, double *, double **, double **);
extern void pdgstrf_WorkFree (int *, double *, GlobalLU_t *);
extern int  pdgstrf_MemXpand (int, int, MemType, int *, GlobalLU_t *);

extern int  *intMalloc (int);
extern int  *intCalloc (int);
extern double *doubleMalloc(int);
extern double *doubleCalloc(int);
extern int  memory_usage ();
extern int  superlu_dQuerySpace (int, SuperMatrix *, SuperMatrix *, int, 
				 superlu_memusage_t *);
extern int  Glu_alloc (const int, const int, const int, const MemType,
		       int *, pxgstrf_shared_t *);

/* -------------------
   Auxiliary routines
   -------------------*/
extern double  SuperLU_timer_();
extern int     sp_ienv(int);
extern double  dlamch_();
extern int     lsame_(char *, char *);
extern int     xerbla_(char *, int *);
extern void    superlu_abort_and_exit(char *);
extern void    ifill(int *, int, int);
extern void    dfill(double *, int, double);
extern void    inf_norm_error(int, SuperMatrix *, double *);
extern void    dstat_allocate(int);
extern void    snode_profile(int, int *);
extern void    super_stats(int, int *, int *);
extern void    panel_stats(int, int, int *, Gstat_t *);
extern void    PrintSumm(char *, int, int, int);
extern void    dPrintPerf(SuperMatrix *, SuperMatrix *, superlu_memusage_t *,
			 double, double, double *, double *, char *,
			 Gstat_t *);
extern void    dCompRow_to_CompCol(int m, int n, int nnz, 
                           double *a, int *colind, int *rowptr,
                           double **at, int **rowind, int **colptr);


/* -----------------------
   Routines for debugging
   -----------------------*/
extern void    print_lu_col(int, char *, int, int, int, int *, GlobalLU_t *);
extern void    print_panel_seg(int, int, int, int, int *, int *);
extern void    dcheck_zero_vec(int, char *, int, double *);
extern void    check_repfnz(int, int, int, int *);

#ifdef __cplusplus
	   }
#endif


#endif /* __SUPERLU_DSP_DEFS */

