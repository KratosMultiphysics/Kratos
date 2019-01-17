/**
 *
 * @file out.h
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX memory tracking function.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 */
#ifndef _out_h_
#define _out_h_

#define OUT_HEADER                                              \
    "+-------------------------------------------------+\n"     \
    "+     PaStiX : Parallel Sparse matriX package     +\n"     \
    "+-------------------------------------------------+\n"     \
    "  Version:                                   %d.%d.%d\n"   \
    "  Schedulers:\n"                                           \
    "    sequential:                           %8s\n"           \
    "    thread static:                        %8s\n"           \
    "    thread dynamic:                       %8s\n"           \
    "    PaRSEC:                               %8s\n"           \
    "    StarPU:                               %8s\n"           \
    "  Number of MPI processes:                %8d\n"           \
    "  Number of threads per process:          %8d\n"           \
    "  Number of GPUs:                         %8d\n"           \
    "  MPI communication support:              %8s\n"           \
    "  Distribution level:               %8s(%4ld)\n"           \
    "  Blocking size (min/max):             %4ld / %4ld\n"      \
    "  Computational models\n"                                  \
    "    CPU: %41s\n"                                           \
    "    GPU: %41s\n"                                           \
    "  Low rank parameters:\n"                                  \
    "    Strategy                      %16s\n"

#define OUT_HEADER_LR                                           \
    "    Tolerance                             %8.0e\n"         \
    "    Compress method                       %8s\n"           \
    "    Compress minimal width                %8ld\n"          \
    "    Compress minimal height               %8ld\n"          \
    "    Compress min ratio                    %8f\n"           \
    "    Orthogonalization method              %8s\n"

#define OUT_STEP_ORDER                                          \
    "+-------------------------------------------------+\n"     \
    "  Ordering step :\n"
#define OUT_SUBSTEP_GRAPH                       \
    "    Prepare graph structure:\n"
#define OUT_ORDER_SYMGRAPH                      \
    "      Symmetrizing graph\n"
#define OUT_ORDER_NODIAG                        \
    "      Removing diagonal elements\n"
#define OUT_ORDER_SORT                          \
    "      Sort row indexes in each column\n"
#define OUT_ORDER_INIT                          \
    "    Compute ordering\n"
#define OUT_ORDER_METHOD                        \
    "    Ordering method is: %s\n"
#define OUT_ORDER_TIME                                  \
    "    Time to compute ordering              %e s\n"

#define OUT_STEP_FAX                                            \
    "+-------------------------------------------------+\n"     \
    "  Symbolic factorization step:\n"
#define OUT_FAX_METHOD                          \
    "    Symbol factorization using: %s\n"
#define OUT_FAX_SUMMARY                                                 \
    "    Number of nonzeroes in L structure    %8ld\n"                  \
    "    Fill-in of L                          %8lf\n"                  \
    "    Time to compute symbol matrix         %e s\n"


#define OUT_STEP_REORDER                                        \
    "+-------------------------------------------------+\n"     \
    "  Reordering step:\n"                                      \
    "    Split level                           %8ld\n"          \
    "    Stoping criteria                      %8ld\n"
#define OUT_REORDERING_TIME                             \
    "    Time for reordering                   %e s\n"
#define OUT_REORDERING_OPS                                              \
    "    Iops for the last supernode           %8ld ( %5.2lf%% )\n"     \
    "    Iops for the reordering               %8ld\n"

#define OUT_STEP_BLEND                                          \
    "+-------------------------------------------------+\n"     \
    "  Analyse step:\n"
#define OUT_BLEND_CONF                                  \
    "    Number of cluster                     %8ld\n"  \
    "    Number of processor per cluster       %8ld\n"  \
    "    Number of thread per MPI process      %8ld\n"

#define OUT_BLEND_CHKSMBMTX                     \
    "    Check the symbol matrix\n"
#define OUT_BLEND_CHKSOLVER                     \
    "    Check the solver matrix\n"
#define OUT_BLEND_ELIMTREE                      \
    "    Building elimination tree\n"
#define OUT_BLEND_ELIMTREE_TIME                         \
    "    Elimination tree built in             %e s\n"
#define OUT_BLEND_COSTMATRIX                    \
    "    Building cost matrix\n"
#define OUT_BLEND_COSTMATRIX_TIME                       \
    "    Cost matrix built in                  %e s\n"
#define OUT_BLEND_ELIMTREE_TOTAL_COST                   \
    "    Total estimated cost of the etree     %e s\n"
#define OUT_BLEND_PROPMAP                       \
    "    Perform proportional mapping\n"
#define OUT_BLEND_PROPMAP_TIME                          \
    "    Proportional mapping done in          %e s\n"
#define OUT_BLEND_SPLITSYMB                     \
    "    Split large symbolic blocks\n"
#define OUT_BLEND_SPLITSYMB_TIME                        \
    "    Symbol split done in                  %e s\n"
#define OUT_BLEND_BUILDSIMU                     \
    "    Build simulation structures\n"
#define OUT_BLEND_BUILDSIMU_TIME                        \
    "    Simulation structures built in        %e s\n"  \
    "    Number of tasks found                 %8ld\n"
#define OUT_BLEND_SIMU                                  \
    "    Start simulation (Data distribution)\n"
#define OUT_BLEND_SIMU_TIME                             \
    "    Simulation done in                    %e s\n"
#define OUT_BLEND_ELIMGRAPH                     \
    "    Building elimination graph\n"
#define OUT_BLEND_ELIMGRAPH_TIME                        \
    "    Elimination graph built in            %e s\n"
#define OUT_BLEND_SOLVER                        \
    "    Building solver structure\n"
#define OUT_BLEND_SOLVER_TIME                           \
    "    Solver built in                       %e s\n"
#define OUT_BLEND_TIME                                  \
    "    Time for analyze                      %e s\n"

#define OUT_BLEND_SUMMARY                                               \
    "    Number of non-zeroes in blocked L     %8ld\n"                  \
    "    Fill-in                               %8lf\n"                  \
    "    Number of operations in full-rank: %-5s    %5.2lf %cFlops\n"   \
    "    Prediction:\n"                                                 \
    "      Model                       %20s\n"                          \
    "      Time to factorize                   %e s\n"                  \
    "    Time for analyze                      %e s\n"

#define OUT_STEP_SOPALIN                                          \
    "+-------------------------------------------------+\n"     \
    "  Factorization step:\n"                                   \
    "    Factorization used: %s\n"

#define OUT_BCSC_TIME                                   \
    "    Time to initialize internal csc       %e s\n"

#define OUT_COEFTAB_TIME                                \
    "    Time to initialize coeftab            %e s\n"

#define OUT_SOPALIN_TIME                                                \
    "    Time to factorize                     %e s (%5.2lf %cFlop/s)\n" \
    "    Number of operations                       %5.2lf %cFlops\n"   \
    "    Number of static pivots               %8ld\n"

#define OUT_LOWRANK_SUMMARY                                     \
    "    Compression:\n"                                        \
    "      Elements removed             %8ld / %8ld\n"          \
    "      Memory saved              %.3g %co / %.3g %co\n"

#define OUT_STARPU_TP         " StarPU : Thread policy : %s\n"
#define OUT_STARPU_STP        " StarPU : No thread policy, setting thread policy to : %s\n"
#define OUT_MATRIX_SIZE       "  Matrix size                                   %ld x %ld\n"
#define OUT_NNZ               "  Number of nonzeros in A                       %ld\n"

#define OUT_GLOBAL_NNZL       "   Number of nonzeroes in L structure      %ld\n"
#define OUT_GLOBAL_FILLIN     "   Fill-in                                 %lf\n"
#define OUT_GLOBAL_THFLOPCNT  "   Number of theoretical flop            %.5g %cflops\n"
#define OUT_GLOBAL_RLFLOPCNT  "   Number of performed flop              %.5g %cflops\n"

#define TIME_TO_ANALYSE       "   Time to analyze                              %.3g s\n"
#define NNZERO_WITH_FILLIN_TH "   Number of nonzeros in factorized matrix      %ld\n"
#define NNZERO_WITH_FILLIN    "%d : Number of nonzeros (local block structure) %ld\n"
#define SOLVMTX_WITHOUT_CO    "%d : SolverMatrix size (without coefficients)   %.3g %s\n"
#define OUT_FILLIN_TH         "   Fill-in                                      %lg\n"
#define NUMBER_OP_LU          "   Number of operations (LU)                    %g\n"
#define NUMBER_OP_LLT         "   Number of operations (LLt)                   %g\n"
#define TIME_FACT_PRED        "   Prediction Time to factorize (%s) %.3g s\n"
#define OUT_COEFSIZE          "   Maximum coeftab size (cefficients)           %.3g %co\n"
#define OUT_REDIS_CSC         "   Redistributing user CSC into PaStiX distribution\n"
#define OUT_REDIS_RHS         "   Redistributing user RHS into PaStiX distribution\n"
#define OUT_REDIS_SOL         "   Redistributing solution into Users' distribution\n"
#define OUT2_SOP_BINITG       "   --- Sopalin : Allocation de la structure globale ---\n"
#define OUT2_SOP_EINITG       "   --- Fin Sopalin Init                             ---\n"
#define OUT2_SOP_TABG         "   --- Initialisation des tableaux globaux          ---\n"
#define OUT2_SOP_BINITL       "   --- Sopalin : Local structure allocation         ---\n"
#define OUT2_SOP_NOTBIND      "   --- Sopalin : Threads are NOT binded             ---\n"
#define OUT2_SOP_BIND         "   --- Sopalin : Threads are binded                 ---\n"
#define OUT2_FUN_STATS        "     - %3ld : Envois %5ld - Receptions %5ld          -\n"
#define OUT2_SOP_BSOP         "   --- Sopalin Begin                                ---\n"
#define OUT2_SOP_ESOP         "   --- Sopalin End                                  ---\n"
#define OUT4_UPDO_TIME_INIT   " [%d][%d] Solve initialization time : %lg s\n"
#define OUT4_UPDO_COMM_TIME   " [%d][%d] Solve communication time : %lg s\n"
#define OUT4_FACT_COMM_TIME   " [%d][%d] Factorization communication time : %lg s\n"
#define OUT2_SOP_DOWN         "   --- Down Step                                    ---\n"
#define OUT2_SOP_DIAG         "   --- Diag Step                                    ---\n"
#define OUT2_SOP_UP           "   --- Up Step                                      ---\n"
#define GEN_RHS_1             "   Generate RHS for X=1\n"
#define GEN_RHS_I             "   Generate RHS for X=i\n"
#define GEN_SOL_0             "   Generate X0=0\n"
#define OUT_ITERREFINE_GMRES    "   GMRES :\n"
#define OUT_ITERREFINE_PIVOT    "   Simple refinement :\n"
#define OUT_ITERREFINE_BICGSTAB "   BICGSTAB :\n"
#define OUT_ITERREFINE_GRAD     "   Conjuguate gradient :\n"
#define OUT_ITERREFINE_ITER     "    - iteration %d :\n"
#define OUT_ITERREFINE_TTS      "         time to solve                          %.3g s\n"
#define OUT_ITERREFINE_TTT      "         total iteration time                   %.3g s\n"
#define OUT_ITERREFINE_ERR      "         error                                  %.5g\n"
#define OUT_ITERREFINE_NORMA    "         ||A||                                  %.5g\n"
#define OUT_ITERREFINE_NORMR    "         ||r||                                  %.5g\n"
#define OUT_ITERREFINE_NORMB    "         ||b||                                  %.5g\n"
#define OUT_ITERREFINE_BDIVR    "         ||r||/||b||                            %.5g\n"
#define OUT_REDISCSCDTIME     "   Time to redistribute cscd                    %.3g s\n"
#define OUT_FILLCSCTIME       "   Time to fill internal csc                    %.3g s\n"
#define OUT_MAX_MEM_AF_SOP    "   Max memory used after factorization          %.3g %s\n"
#define OUT_MEM_USED_AF_SOP   "   Memory used after factorization              %.3g %s\n"
#define MAX_MEM_AF_CL         "   Max memory used after clean                  %.3g %s\n"
#define MEM_USED_AF_CL        "   Memory used after clean                      %.3g %s\n"
#define OUT_STATIC_PIVOTING   "   Static pivoting                              %ld\n"
#define OUT_ESP_NBTASKS       "   Number of tasks added by esp                 %ld\n"
#define OUT_TIME_FACT         "   Time to factorize                            %.3g s  (%.3g %s)\n"
#define OUT_FLOPS_FACT        "   FLOPS during factorization                   %.5g %s\n"
#define OUT_TIME_SOLV         "    Time to solve                         %e s\n"
#define OUT_REFINE_ITER_NORM  "    Refinement                            %ld iterations, norm=%e\n"
#define OUT_PREC1             "    ||b-Ax||/||b||                        %e\n"
#define OUT_PREC2             "    max_i(|b-Ax|_i/(|b| + |A||x|)_i       %e\n"
#define OUT_TIME_REFINE       "    Time for refinement                   %e s\n"
#define OUT_END               " +--------------------------------------------------------------------+\n"

/*
 * Printing function to redirect to the correct file
 */
#if defined(__GNUC__)
static inline void pastix_print( int mpirank, int thrdrank, const char *fmt, ...) __attribute__((format(printf,3,4)));
static inline void pastix_print_error  ( const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void pastix_print_warning( const char *fmt, ...) __attribute__((format(printf,1,2)));
#endif

static inline void
pastix_print( int mpirank, int thrdrank, const char *fmt, ...)
{
    va_list ap;

    if( (mpirank == 0) && (thrdrank == 0) )
    {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap );
        va_end(ap);
    }
}

static inline void
pastix_print_error( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

static inline void
pastix_print_warning( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

#define errorPrint  pastix_print_error
#define errorPrintW pastix_print_warning

static inline double
pastix_print_value( double flops )
{
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return flops;
}

static inline char
pastix_print_unit( double flops )
{
    static char units[9] = { ' ', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' };
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return units[unit];
}

static inline char *
pastixFactotypeStr( pastix_factotype_t ft ) {
    switch( ft ) {
    case PastixFactLLT:
        return "LL^t";
    case PastixFactLDLT:
        return "LDL^t";
    case PastixFactLU:
        return "LU";
    case PastixFactLLH:
       return "LL^h";
    case PastixFactLDLH:
       return "LDL^h";
     default:
        return "None";
    }
}

void   pastix_gendirtemp( char **dirtemp );
FILE * pastix_fopenw( char       **directory,
                      const char  *filename,
                      const char  *mode );
FILE * pastix_fopen ( const char  *filename );

#endif /* _out_h_ */
