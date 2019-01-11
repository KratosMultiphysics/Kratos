/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */

#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "slu_mt_machines.h"

/* Macros */
#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define SUPERLU_ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#define USER_MALLOC(size) superlu_malloc(size)
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)

#define MAX(x, y) 	   ( (x) > (y) ? (x) : (y) )
#define MIN(x, y) 	   ( (x) < (y) ? (x) : (y) )
#define SUPERLU_MAX(x, y)  ( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y)  ( (x) < (y) ? (x) : (y) )

/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colbeg[col] )
#define L_SUB_END(col)       ( Lstore->rowind_colend[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colbeg[col] )
#define L_NZ_END(col)        ( Lstore->nzval_colend[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_colbeg[superno] )
#define L_LAST_SUPC(superno) ( Lstore->sup_to_colend[superno] )
#define U_NZ_START(col)      ( Ustore->colbeg[col] )
#define U_NZ_END(col)        ( Ustore->colend[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )

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

/* 
 * Constants 
 */
#define EMPTY	(-1)
#define FALSE	0
#define TRUE	1

/**********************
  Enumerated constants
  *********************/
typedef enum {NO, YES}                       yes_no_t;
typedef enum {NOTRANS, TRANS, CONJ}          trans_t;
typedef enum {DOFACT, EQUILIBRATE, FACTORED} fact_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
	      METIS_AT_PLUS_A, PARMETIS, MY_PERMC} colperm_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}       equed_t;
typedef enum {LUSUP, UCOL, LSUB, USUB}       MemType;

/* Number of marker arrays used in the symbolic factorization, 
   each of size nrow. */
#define NO_MARKER     3

#define LOCOL    70
#define HICOL    78
#define BADROW   44
#define BADCOL   35
#define BADPAN   BADCOL
#define BADREP   35

/*
 * Type definitions
 */
typedef float    flops_t;
typedef unsigned char Logical;

#if ( MACH==DEC || MACH==PTHREAD )
#include <pthread.h>
typedef pthread_mutex_t mutex_t;
#elif ( MACH==SGI || MACH==ORIGIN )
typedef int mutex_t;
#elif ( MACH==CRAY_PVP || MACH==OPENMP )
typedef int mutex_t;
#endif


typedef enum {
    RELAX,
    COLPERM,
    ETREE,
    EQUIL,
    FINDDOMAIN,
    FACT,
    DFS,
    FLOAT,
    TRSV,
    GEMV,
    RCOND,
    TRISOLVE,
    SOLVE,
    REFINE,
    FERR,
    NPHASES
} PhaseType;

/* 
 * *********************************************************************
 * The superlumt_options_t structure contains the shared variables used 
 * for factorization, which are passed to each thread.
 * *********************************************************************
 * 
 * nprocs (int)
 *        Number of processes (or threads) to be spawned and used to perform
 *        the LU factorization by pdgstrf().
 *
 * fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, whether the matrix A should
 *        be equilibrated before it is factored.
 *        = DOFACT: The matrix A will be factored, and the factors will be
 *              stored in L and U.
 *        = EQUILIBRATE: The matrix A will be equilibrated if necessary, then
 *              factored into L and U.
 *
 * trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * refact (yes_no_t)
 *        Specifies whether this is first time or subsequent factorization.
 *        = NO:  this factorization is treated as the first one;
 *        = YES: it means that a factorization was performed prior to this
 *               one. Therefore, this factorization will re-use some
 *               existing data structures, such as L and U storage, column
 *               elimination tree, and the symbolic information of the
 *               Householder matrix.
 *
 * panel_size (int)
 *        A panel consists of at most panel_size consecutive columns.
 *
 * relax  (int)
 *        To control degree of relaxing supernodes. If the number
 *        of nodes (columns) in a subtree of the elimination tree is less
 *        than relax, this subtree is considered as one supernode,
 *        regardless of the row structures of those columns.
 *
 * diag_pivot_thresh (double)
 *        Diagonal pivoting threshold. At step j of the Gaussian elimination,
 *        if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *        use A_jj as pivot, else use A_ij with maximum magnitude. 
 *        0 <= diag_pivot_thresh <= 1. The default value is 1, 
 *        corresponding to partial pivoting.
 *
 * drop_tol (double) (NOT IMPLEMENTED)
 *	  Drop tolerance parameter. At step j of the Gaussian elimination,
 *        if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
 *        0 <= drop_tol <= 1. The default value of drop_tol is 0,
 *        corresponding to not dropping any entry.
 *
 * usepr  (yes_no_t)
 *        Whether the pivoting will use perm_r specified by the user.
 *        = YES: use perm_r; perm_r is input, unchanged on exit.
 *        = NO:  perm_r is determined by partial pivoting, and is output.
 *
 * SymmetricMode (yest_no_t)
 *        Specifies whether to use symmetric mode.
 *
 * PrintStat (yes_no_t)
 *        Specifies whether to print solver's statistics.
 *
 * perm_c (int*) dimension A->ncol
 *	  Column permutation vector, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *        When search for diagonal, perm_c[*] is applied to the
 *        row subscripts of A, so that diagonal threshold pivoting
 *        can find the diagonal of A, instead of that of A*Pc.
 *
 * perm_r (int*) dimension A->nrow
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *        If usepr = NO, perm_r is output argument;
 *        If usepr = YES, the pivoting routine will try to use the input
 *           perm_r, unless a certain threshold criterion is violated.
 *           In that case, perm_r is overwritten by a new permutation
 *           determined by partial pivoting or diagonal threshold pivoting.
 *
 * work   (void*) of size lwork
 *        User-supplied work space and space for the output data structures.
 *        Not referenced if lwork = 0;
 *
 * lwork  (int)
 *        Specifies the length of work array.
 *        = 0:  allocate space internally by system malloc;
 *        > 0:  use user-supplied work array of length lwork in bytes,
 *              returns error if space runs out.
 *        = -1: the routine guesses the amount of space needed without
 *              performing the factorization, and returns it in
 *              superlu_memusage->total_needed; no other side effects.
 *
 * etree  (int*)
 *        Elimination tree of A'*A, dimension A->ncol.
 *        Note: etree is a vector of parent pointers for a forest whose
 *        vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *        On input, the columns of A should be permutated so that the
 *        etree is in a certain postorder.
 *
 * colcnt_h (int*)
 *        Column colunts of the Householder matrix.
 *
 * part_super_h (int*)
 *        Partition of the supernodes in the Householder matrix.
 *	  part_super_h[k] = size of the supernode beginning at column k;
 * 	                  = 0, elsewhere.
 *
 *
 */
typedef struct {
    int_t        nprocs;
    fact_t     fact;
    trans_t    trans;
    yes_no_t   refact;
    int_t        panel_size;
    int_t        relax;
    double     diag_pivot_thresh;
    double     drop_tol;
    colperm_t  ColPerm;
    yes_no_t   usepr;
    yes_no_t   SymmetricMode;
    yes_no_t   PrintStat;

    /* The following arrays are persistent during repeated factorizations. */
    int_t  *perm_c;
    int_t  *perm_r;
    void *work;
    int_t  lwork;

    /* The following structural arrays are computed internally by 
       sp_colorder(). The user needs to allocate space on input.
       These 3 arrays are computed in the first factorization, and are 
       re-used in the subsequent factors of the matrices with the same
       nonzero structure. */
    int_t  *etree;
    int_t  *colcnt_h;
    int_t  *part_super_h;
} superlumt_options_t;

/* ----------------------------------------------
    The definitions below are used for profiling.
   ---------------------------------------------- */

/* The statistics to be kept by each processor. */
typedef struct {
    int_t	    panels;    /* number of panels taken */
    float   fcops;     /* factor floating-point operations */
    double  fctime;    /* factor time */
    int_t     skedwaits; /* how many times the processor fails to get a task */
    double  skedtime;  /* time spent in the scheduler */
    double  cs_time;   /* time spent in the critical sections */
    double  spintime;  /* spin-wait time */
    int_t     pruned;
    int_t     unpruned;
} procstat_t;


/* Statistics about each panel. */

typedef struct {
    int_t    size;      /* size of the panel */
    int_t    pnum;      /* which processor grabs this panel */
    double starttime; /* at waht time this panel is assigned to a proc */
    double fctime;    /* factorization time */
    float  flopcnt;   /* floating-point operations */
    int_t    pipewaits; /* how many times the panel waited during pipelining */
    double spintime;  /* spin waiting time during pipelining */
} panstat_t;

/* How was a panel selected by the scheduler */
typedef enum {NOPIPE, DADPAN, PIPE} how_selected_t;

/* Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    int_t size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

/* The structure to keep track of memory usage. */
typedef struct {
    float for_lu;
    float total_needed;
    int_t   expansions;
} superlu_memusage_t;

typedef struct {
     flops_t flops;
     int_t     nzs;
     double  fctime;
} stat_relax_t;

typedef struct {
     flops_t flops;
     int_t nzs;
     double fctime;
} stat_col_t;

typedef struct {
     int_t ncols;
     flops_t flops;
     int_t nzs;
     double fctime;
} stat_snode_t;

/* -------------------------------------------------------------------
   The definitions below are used to simulate parallel execution time.
   ------------------------------------------------------------------- */
typedef struct {
    float est;  /* earliest (possible) start time of the panel */
    float pdiv; /* time in flops spent in the (inner) panel factorization */
} cp_panel_t;

typedef struct {
    float eft;  /* earliest finishing time */
    float pmod; /* pmod update to the ancestor panel */
} desc_eft_t;
		   
/* All statistics. */
typedef struct {
    int_t     	*panel_histo;	/* Panel size distribution */
    double  	*utime;
    flops_t 	*ops;
    procstat_t 	*procstat;
    panstat_t	*panstat;
    int_t      	num_panels;
    float     	dom_flopcnt;
    float     	flops_last_P_panels;
    /**/
    stat_relax_t *stat_relax;
    stat_col_t *stat_col; 
    stat_snode_t *stat_snode; 
    int_t *panhows;
    cp_panel_t *cp_panel; /* panels on the critical path */
    desc_eft_t *desc_eft; /* all we need to know from descendants */
    int_t        *cp_firstkid, *cp_nextkid; /* linked list of children */
    int_t        *height;
    float      *flops_by_height;
} Gstat_t;

struct Branch {
    int_t root, first_desc, which_bin;
    struct Branch *next;
};


#if 0

/* Statistics for supernode and panel size */
int_t 	no_panels;
float   sum_w;          /* Sum (Wi) */
float 	sum_np_w;       /* Sum (Npi*Wi) */
int_t 	max_np;          
int_t     no_sups;
float   sum_sup;        /* Sum (Supi) */
int_t     max_sup;     
flops_t reuse_flops;    /* Triangular solve and matrix vector multiply */
float   reuse_data;     /* Doubles in updating supernode */

/* Statistics for blas operations */
int_t     num_blas;       /* no of BLAS2 operations, including trsv/gemv */
int_t     max_blas_n;     /* max dimension n in tri-solve and mat-vec */
int_t     min_blas_n;     /* min dimension n in tri-solve and mat-vec */
float   sum_blas_n;     /* sum of "        "        " */
int_t     max_gemv_m;     /* max dimension m in mat-vec */
int_t     min_gemv_m;     /* max dimension m in mat-vec */
float   sum_gemv_m;     /* sum of "        "        " */
int_t     lda_blas_m;
int_t     lda_blas_n;
flops_t *gemv_ops;      /* flops distribution on (m,n) */
flops_t *trsv_ops;      /* flops distribution on n */

#define i_trsv_ops(i)      trsv_ops[i]
#define ij_gemv_ops(i,j)   gemv_ops[j*lda_blas_m + i]

#endif


/* *********************
   Function prototypes
   *********************/

#ifdef __cplusplus
extern "C" {
#endif

extern int  cpp_defs();
extern int  xerbla_ (char *, int *);
extern void superlu_abort_and_exit(char*);
extern void *superlu_malloc (size_t);
extern void superlu_free (void*);
extern void PrintStat(Gstat_t *);
extern int_t  ParallelProfile(const int_t, const int_t, const int_t,
			    const int_t procs, Gstat_t *);

#ifdef __cplusplus
	   }
#endif

#endif /* __SUPERLU_UTIL */

