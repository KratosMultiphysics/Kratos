#define MAX_MAT	   100
#define MAX_LINE        256
#define MaxNamLen 64
#define HB   1
#define MM0  2
#define MM1  3
#define UNK  4

typedef struct _io_t {
    FILE *fout;                 /* output file handle              */
    char outfile[MAX_LINE];     /* output filename                 */
    char Fname[MAX_LINE];       /* full matrix path name           */
    char MatNam[MaxNamLen];     /* short name                      */
    char PrecMeth[MAX_LINE];    /* preconditioner being tested     */
    char type[4];               /* type for HB matrices            */
    int Fmt;                    /* matrix format type              */
    int ndim;                   /* matrix size                     */
    int nnz;                    /* number of nonzero               */
/* parameters from inputs -----------------------------------------*/
    int im;                     /* Dim of Krylov subspace [fgmr]   */
    int maxits;                 /* maximum number of fgmres iters  */
    double tol;                 /* tolerance for stopping fgmres   */
    double eps;   /* for checking how close two rows of matrix are */
    int nparam;         /* number of tests for each preconditioner */
    int lfil0;                  /* initial lfil                    */
    int lfilInc;                /* increment for lfil              */
    double tol0;                /* initial drop tolerance          */
    double tolMul;              /* multiplier for tol              */
    int fill_lev;               /* initial level of fill for ILUK  */
    int fill_lev_inc;           /* increment for level of fill for ILUK */
                                /* value always set to 1           */
   int perm_type;               /* indset perms (0) or PQ perms (1)*/
                                /*                  or coarsen (2) */
   int Bsize;                   /* block size - dual role. see input file */
                                /* for explanations */
/* result for output ----------------------------------------------*/
    double rt_v;                /* compression rate of vertices    */
    double rt_e;                /* compression rate of edges       */
    double ceff;                /* compression efficiency          */
    double tm_h;                /* time for hash method  [vbilu]   */
    double tm_a;                /* time for angle method [vbilu]   */
    double tm_b;                /* time for initial blocks (s)     */
    double tm_p;                /* time for preconditioner (s)     */
    double tm_i;                /* time for iteration (s)          */
    double fillfact;            /* memory used during precondition */
    int its;                    /* number of iterations            */
    double enorm;               /* error norm:          || x- x0|| */
    double rnorm;               /* final residual norm: ||Ax-Ax0|| */
} io_t;

