#ifndef __VBLOCK_HEADER_H__
#define __VBLOCK_HEADER_H__

#define MAX_BLOCK_SIZE   100

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define B_DIM(bs,i)      (bs[i+1]-bs[i])

typedef struct SpaFmt {
/*--------------------------------------------- 
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
  int n;
  int *nzcount;  /* length of each row */
  int **ja;      /* pointer-to-pointer to store column indices  */
  double **ma;   /* pointer-to-pointer to store nonzero entries */
} SparMat, *csptr;

typedef double *BData;

typedef struct VBSpaFmt {
    int n;        /* the block row dimension of the matrix      */
    int *bsz;     /* the row/col of the first element of each   */
                  /* diagonal block                             */
    int *nzcount;  /* length of each row                         */
    int **ja;     /* pointer-to-pointer to store column indices */
    BData **ba;   /* pointer-to-pointer to store nonzero blocks */
    BData *D;     /* to store inversion of diagonals            */
} VBSparMat, *vbsptr;

typedef struct VBILUfac {
    int n;
    int *bsz;     /* the row/col of the first element of each   */
                  /* diagonal block                             */

    BData *D;     /* diagonal blocks                            */
    vbsptr L;     /* L part blocks                              */
    vbsptr U;     /* U part blocks                              */
    int *work;    /* working buffer                             */
    BData bf;     /* buffer of a temp block                     */
    int DiagOpt;  /* Option for diagonal inversion/solutiob     */
                  /* opt =  1 -->> call luinv                   */
                  /* opt == 2 -->> block inverted call dgemv    */
} VBILUSpar, *vbiluptr; 

typedef struct ILUfac {
    int n;
    csptr L;      /* L part elements                            */
    double *D;    /* diagonal elements                          */
    csptr U;      /* U part elements                            */
    int *work;    /* working buffer */
} ILUSpar, LDUmat, *iluptr;

typedef struct PerMat4 *p4ptr;
typedef struct PerMat4 {
/*------------------------------------------------------------
| struct for storing the block LU factorization 
| contains all the block factors except the 
| data related to the last block. 
| n       = size of current block
| symperm = whether or not permutations are symmetric.
|           used only in cleanP4..
| nB      = size of B-block
| L, U    = ILU factors of B-block
| F, E    = sparse matrices in (1,2) and (2,1) 
|           parts of matrix. 
| perm    = (symmetric) permutation used in factorization
|           comes from the independent set ordering
| rperm   = unsymmetric permutation (rows) not used in this
|           version -- but left here for compatibility..
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|----------------------------------------------------------*/ 
  int n;                  
  int nB; 
  int symperm;
/*   LU factors  */
  struct SpaFmt *L;
  struct SpaFmt *U;
/* E, F blocks   */
  struct SpaFmt *E;
  struct SpaFmt *F;
  int *rperm;       /* row permutation         */ 
  int *perm;        /* col. permutation        */ 
  double *D1 ;      /* diagonal scaling row    */  
  double *D2 ;      /* diagonal scaling columns*/  
  double *wk;       /* work array              */
/* pointer to next and previous struct         */
  p4ptr prev; 
  p4ptr next;
} Per4Mat; 
/* -------------------------------------------------------------------*/
typedef struct ILUTfac *ilutptr;
typedef struct ILUTfac {
/*------------------------------------------------------------
| struct for storing data related to the last schur complement 
| we need to store the C matrix associated with the last block
| and the ILUT factorization of the related Schur complement.
| 
| n       = size of C block = size of Schur complement
| C       = C block of last level matrix. 
| L, U    = ILU factors of last schur complement. 
|
| meth[4] = parameters for defining variants in factorization 
|           - see function readin for details
| rperm    = row permutation used for very nonsymmetric matrices 
|            [such as bottleneck transversal] -- NOT IN THIS VERSION
| perm2     = unsymmetric permutation (columns) - used primarily
|           for the ILUTP version of ILUT/.. 
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|-----------------------------------------------------------*/
   int n;                  
 /*-------------------- C matrix of last block */
   struct SpaFmt *C;
 /* LU factorization       */
   struct SpaFmt *L;
   struct SpaFmt *U;
 /*--------------------  related to variants and methods */
 /*    int meth[4];   */
   int *rperm;   /* row single-sinded permutation */
   int *perm;    /* column perm .                */
   int *perm2;   /* column permutation coming from pivoting in ILU */ 
   double *D1;
   double *D2;
   double *wk;
} IluSpar;

typedef struct arms_st *arms;
typedef struct arms_st {
  /* this is the arms preconditioner struct 
  | it consists of a linked list of p4mat structs
  | and the ILUt factorization (in the  form of an 
  | IluSpar struct  
  |---------------------------------------------- */
  int n;                   /* dimension of matrix */
  int nlev;                /* number of levels    */
  ilutptr ilus;            /* ILU for last level  */
  p4ptr levmat;            /* level structure     */
} armsMat;

typedef struct __CompressType
{
  int grp;   /* -1: begin new group, >=0: belong to grp-th row */
  int count; /* block size, valid only if grp = -1 */
} CompressType;

typedef struct _SMat {
  /*-------------------- 3 types of matrices so far */
  int n; 
  int Mtype;           /*--  type 1 = CSR, 2 = VBCSR, 3 = LDU    */
  csptr CS;            /* place holder for a CSR/CSC type matrix */
  iluptr LDU;          /* struct for an LDU type matrix          */
  vbsptr VBCSR;        /* place holder for a block matrix        */
  void (*matvec)(struct _SMat*, double *, double *);
} SMat, *SMatptr;

typedef struct _SPre {
  /*-------------------- 3 types of matrices so far */
  int Ptype;           /*-- Ptype 1 = ILU, 2 = VBILU, 3 = Crout */
  iluptr   ILU;        /* struct for an ILU type preconditioner */
  vbiluptr VBILU;      /* struct for a block preconditioner */
  arms ARMS;           /* struct for a block preconditioner */
  int (*precon) (double *, double *, struct _SPre*); 
} SPre, *SPreptr;
  

#endif  /* __VBLOCK_HEADER_H__ */

