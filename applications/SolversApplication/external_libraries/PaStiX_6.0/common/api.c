/**
 *
 * @file api.c
 *
 * PaStiX API routines
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#define _GNU_SOURCE 1
#include "common.h"
#if defined(PASTIX_ORDERING_METIS)
#include <metis.h>
#endif
#include "graph.h"
#include "pastix/order.h"
#include "solver.h"
#include "bcsc.h"
#include "isched.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "models.h"
#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif
#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Generate a unique temporary directory to store output files
 *
 *******************************************************************************
 *
 * @param[inout] dirtemp
 *          On entry, if dirtemp is not initizalized, the unique string is
 *          generated, otherwise nothing is done.
 *
 *******************************************************************************/
void
pastix_gendirtemp( char **dirtemp )
{
    if ( *dirtemp == NULL ) {
        mode_t old_mask = umask(S_IWGRP | S_IWOTH);

        *dirtemp = strdup( "pastix-XXXXXX" );
        *dirtemp = mkdtemp( *dirtemp );
        (void)umask(old_mask);

        if ( *dirtemp == NULL ) {
            errorPrint("pastix_gendirtemp: Couldn't not generate the tempory directory to store the output files");
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Open a file in the unique directory of the pastix instance
 *
 *******************************************************************************
 *
 * @param[inout] dirtemp
 *          The pointer to the directory string associated to the instance.
 *          If not initizalized, initialized on exit.
 *
 * @param[in] filename
 *          The filename to create in the unique directory.
 *
 * @param[in] mode
 *          Opening mode of the file as described by the fopen() function.
 *
 *******************************************************************************
 *
 * @return The FILE structure of the opened file. NULL if failed.
 *
 *******************************************************************************/
FILE *
pastix_fopenw( char       **dirtemp,
               const char  *filename,
               const char  *mode )
{
    char *fullname;
    FILE *f = NULL;
    int rc;

    pastix_gendirtemp( dirtemp );
    if ( *dirtemp == NULL ) {
        return NULL;
    }

    rc = asprintf( &fullname, "%s/%s", *dirtemp, filename );
    if (rc <= 0 ) {
        errorPrint("pastix_fopenw: Couldn't not generate the tempory filename for the output file");
        return NULL;
    }

    f = fopen(fullname, mode);
    free( fullname );

    if (NULL == f)
    {
        perror("pastix_fopenw");
        errorPrint( "pastix_fopenw: Couldn't open file: %s with mode %s\n",
                    filename, mode );
        return NULL;
    }

    return f;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Open a file in the current directory in read only mode.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename of the file to open.
 *
 *******************************************************************************
 *
 * @return The FILE structure of the opened file. NULL if failed.
 *
 *******************************************************************************/
FILE *
pastix_fopen( const char *filename )
{
    FILE *f = fopen(filename, "r");

    if (NULL == f)
    {
        perror("pastix_fopen");
        errorPrint( "pastix_fopen: Couldn't open file: %s with mode r\n",
                    filename );
        return NULL;
    }

    return f;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Print information about the solver configuration
 *
 *******************************************************************************
 *
 * @param[in] pastix
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixWelcome( const pastix_data_t *pastix )
{
    pastix_print( pastix->procnum, 0, OUT_HEADER,
                  /* Version    */ PASTIX_VERSION_MAJOR, PASTIX_VERSION_MINOR, PASTIX_VERSION_MICRO,
                  /* Sched. seq */ "Enabled",
                  /* Sched. sta */ (pastix->isched ? "Started" : "Disabled"),
                  /* Sched. dyn */ "Disabled",
                  /* Sched. PaR */
#if defined(PASTIX_WITH_PARSEC)
                  (pastix->parsec == NULL ? "Enabled" : "Started" ),
#else
                  "Disabled",
#endif
                  /* Sched. SPU */
#if defined(PASTIX_WITH_STARPU)
                  (pastix->starpu == NULL ? "Enabled" : "Started" ),
#else
                  "Disabled",
#endif
                  /* MPI nbr    */ pastix->procnbr,
                  /* Thrd nbr   */ (int)(pastix->iparm[IPARM_THREAD_NBR]),
                  /* GPU nbr    */ (int)(pastix->iparm[IPARM_GPU_NBR]),
#if defined(PASTIX_WITH_MPI)
                  /* MPI mode   */ ((pastix->iparm[IPARM_THREAD_COMM_MODE] == PastixThreadMultiple) ? "Multiple" : "Funneled"),
#else
                  "Disabled",
#endif
                  /* Distrib    */ ((pastix->iparm[IPARM_TASKS2D_LEVEL] == 0) ? "1D" : "2D"),
                                   ((pastix->iparm[IPARM_TASKS2D_LEVEL] <  0) ? ((long)pastix->iparm[IPARM_TASKS2D_WIDTH]) :
                                                                               -((long)pastix->iparm[IPARM_TASKS2D_LEVEL])),
                  /* Block size */ (long)pastix->iparm[IPARM_MIN_BLOCKSIZE],
                                   (long)pastix->iparm[IPARM_MAX_BLOCKSIZE],
                  /* Model CPU  */ pastix->cpu_models->name,
                  /* Model GPU  */ pastix->gpu_models->name,
                  /* Strategy   */ ((pastix->iparm[IPARM_COMPRESS_WHEN] == PastixCompressNever) ? "No compression" :
                                    (pastix->iparm[IPARM_COMPRESS_WHEN] == PastixCompressWhenBegin) ? "Memory Optimal" : "Just-In-Time") );


    if ( pastix->iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever ) {
        pastix_print( pastix->procnum, 0, OUT_HEADER_LR,
                      /* Tolerance       */ (double)pastix->dparm[DPARM_COMPRESS_TOLERANCE],
                      /* Compress method */ ((pastix->iparm[IPARM_COMPRESS_METHOD] == PastixCompressMethodSVD) ? "SVD" : "RRQR"),
                      /* Compress width  */ (long)pastix->iparm[IPARM_COMPRESS_MIN_WIDTH],
                      /* Compress height */ (long)pastix->iparm[IPARM_COMPRESS_MIN_HEIGHT],
                      /* Min ratio       */ (double)pastix->dparm[DPARM_COMPRESS_MIN_RATIO],
                      /* Ortho    method */ ((pastix->iparm[IPARM_COMPRESS_ORTHO] == PastixCompressOrthoCGS) ? "CGS" : (pastix->iparm[IPARM_COMPRESS_ORTHO] == PastixCompressOrthoQR) ? "QR" : "partialQR") );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Print summary information
 *
 *******************************************************************************
 *
 * @param[in] pastix
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixSummary( const pastix_data_t *pastix )
{
    (void)pastix;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the iparm and dparm arrays to their default
 * values.
 *
 * This is performed only if iparm[IPARM_MODIFY_PARAMETER] is set to 0.
 *
 *******************************************************************************
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 *******************************************************************************/
void
pastixInitParam( pastix_int_t *iparm,
                 double       *dparm )
{
    memset( iparm, 0, IPARM_SIZE * sizeof(pastix_int_t) );
    memset( dparm, 0, DPARM_SIZE * sizeof(double) );

    iparm[IPARM_VERBOSE]               = PastixVerboseNo;
    iparm[IPARM_IO_STRATEGY]           = PastixIONo;

    /* Stats */
    iparm[IPARM_NNZEROS]               = 0;
    iparm[IPARM_NNZEROS_BLOCK_LOCAL]   = 0;
    iparm[IPARM_ALLOCATED_TERMS]       = 0;
    iparm[IPARM_PRODUCE_STATS]         = 0;

    /* Scaling */
    iparm[IPARM_MC64]                  = 0;

    /*
     * Ordering parameters
     */
    iparm[IPARM_ORDERING]              = PastixOrderScotch;
    iparm[IPARM_ORDERING_DEFAULT]      = 1;

    /* Scotch */
    {
        iparm[IPARM_SCOTCH_SWITCH_LEVEL] = 120;
        iparm[IPARM_SCOTCH_CMIN]         = 0;
        iparm[IPARM_SCOTCH_CMAX]         = 100000;
        iparm[IPARM_SCOTCH_FRAT]         = 8;
    }

    /* Metis */
    {
#if defined(PASTIX_ORDERING_METIS)
        iparm[IPARM_METIS_CTYPE   ] = METIS_CTYPE_SHEM;
        iparm[IPARM_METIS_RTYPE   ] = METIS_RTYPE_SEP1SIDED;
#else
        iparm[IPARM_METIS_CTYPE   ] = -1;
        iparm[IPARM_METIS_RTYPE   ] = -1;
#endif
        iparm[IPARM_METIS_NO2HOP  ] = 0;
        iparm[IPARM_METIS_NSEPS   ] = 1;
        iparm[IPARM_METIS_NITER   ] = 10;
        iparm[IPARM_METIS_UFACTOR ] = 200;
        iparm[IPARM_METIS_COMPRESS] = 1;
        iparm[IPARM_METIS_CCORDER ] = 0;
        iparm[IPARM_METIS_PFACTOR ] = 0;
        iparm[IPARM_METIS_SEED    ] = 3452;
        iparm[IPARM_METIS_DBGLVL  ] = 0;
    }

    /* Symbolic factorization */
    iparm[IPARM_SF_KASS]               = 0;
    iparm[IPARM_AMALGAMATION_LVLCBLK]  = 5;
    iparm[IPARM_AMALGAMATION_LVLBLAS]  = 5;

    /* Reordering */
    iparm[IPARM_REORDERING_SPLIT]      = 0;
    iparm[IPARM_REORDERING_STOP]       = PASTIX_INT_MAX;

    /* Analyze */
    iparm[IPARM_MIN_BLOCKSIZE]         = 160;
    iparm[IPARM_MAX_BLOCKSIZE]         = 320;
    iparm[IPARM_TASKS2D_LEVEL]         = -1;
    iparm[IPARM_TASKS2D_WIDTH]         = iparm[IPARM_MIN_BLOCKSIZE];
    iparm[IPARM_ABS]                   = 0;
    iparm[IPARM_ALLCAND]               = 0;

    /* Incomplete */
    iparm[IPARM_INCOMPLETE]            = 0;
    iparm[IPARM_LEVEL_OF_FILL]         = 0;

    /* Factorization */
    iparm[IPARM_FACTORIZATION]         = PastixFactLU;
    iparm[IPARM_STATIC_PIVOTING]       = 0;
    iparm[IPARM_FREE_CSCUSER]          = 0;
    iparm[IPARM_SCHUR_FACT_MODE]       = PastixFactModeLocal;

    /* Solve */
    iparm[IPARM_SCHUR_SOLV_MODE]       = PastixSolvModeLocal;
    iparm[IPARM_APPLYPERM_WS]          = 1;

    /* Refinement */
    iparm[IPARM_REFINEMENT]            = PastixRefineGMRES;
    iparm[IPARM_NBITER]                = 0;
    iparm[IPARM_ITERMAX]               = 250;
    iparm[IPARM_GMRES_IM]              = 25;

    /* Context */
    iparm[IPARM_SCHEDULER]             = PastixSchedStatic;
    iparm[IPARM_THREAD_NBR]            = -1;
    iparm[IPARM_AUTOSPLIT_COMM]        = 0;

    /* GPU */
    iparm[IPARM_GPU_NBR]               = 0;
    iparm[IPARM_GPU_MEMORY_PERCENTAGE] = 95;
    iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] = 32 * 1024;

    /* Compression */
    iparm[IPARM_COMPRESS_MIN_WIDTH]    = 120;
    iparm[IPARM_COMPRESS_MIN_HEIGHT]   = 20;
    iparm[IPARM_COMPRESS_WHEN]         = PastixCompressNever;
    iparm[IPARM_COMPRESS_METHOD]       = PastixCompressMethodRRQR;
    iparm[IPARM_COMPRESS_ORTHO]        = PastixCompressOrthoCGS;

    /* MPI modes */
#if defined(PASTIX_WITH_MPI)
    {
        int flag = 0;
        int provided = MPI_THREAD_SINGLE;
        MPI_Initialized(&flag);

        iparm[IPARM_THREAD_COMM_MODE] = 0;
        if (flag) {
            MPI_Query_thread(&provided);
            switch( provided ) {
            case MPI_THREAD_MULTIPLE:
                iparm[IPARM_THREAD_COMM_MODE] = PastixThreadMultiple;
                break;
            case MPI_THREAD_SERIALIZED:
            case MPI_THREAD_FUNNELED:
                iparm[IPARM_THREAD_COMM_MODE] = PastixThreadFunneled;
                break;
                /*
                 * In the folowing cases, we consider that any MPI implementation
                 * should provide enough level of parallelism to turn in Funneled mode
                 */
            case MPI_THREAD_SINGLE:
            default:
                iparm[IPARM_THREAD_COMM_MODE] = PastixThreadFunneled;
            }
        }
    }
#endif /* defined(PASTIX_WITH_MPI) */

    /* Subset for old pastix interface  */
    iparm[IPARM_MODIFY_PARAMETER] = 1;
    iparm[IPARM_START_TASK] = PastixTaskOrdering;
    iparm[IPARM_END_TASK]   = PastixTaskClean;
    iparm[IPARM_FLOAT]      = -1;
    iparm[IPARM_MTX_TYPE]   = -1;
    iparm[IPARM_DOF_NBR]    = 1;

    dparm[DPARM_FILL_IN]            =  0.;
    dparm[DPARM_EPSILON_REFINEMENT] = -1.;
    dparm[DPARM_RELATIVE_ERROR]     = -1.;
    dparm[DPARM_EPSILON_MAGN_CTRL]  =  0.;
    dparm[DPARM_ANALYZE_TIME]       =  0.;
    dparm[DPARM_PRED_FACT_TIME]     =  0.;
    dparm[DPARM_FACT_TIME]          =  0.;
    dparm[DPARM_SOLV_TIME]          =  0.;
    dparm[DPARM_FACT_FLOPS]         =  0.;
    dparm[DPARM_FACT_THFLOPS]       =  0.;
    dparm[DPARM_FACT_RLFLOPS]       =  0.;
    dparm[DPARM_SOLV_FLOPS]         =  0.;
    dparm[DPARM_SOLV_THFLOPS]       =  0.;
    dparm[DPARM_SOLV_RLFLOPS]       =  0.;
    dparm[DPARM_REFINE_TIME]        =  0.;
    dparm[DPARM_A_NORM]             = -1.;
    dparm[DPARM_COMPRESS_TOLERANCE] = 0.01;
    dparm[DPARM_COMPRESS_MIN_RATIO] =  1.;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Internal function that setups the multiple communicators in order
 * to perform the ordering step in MPI only mode, and the factorization in
 * MPI+Thread mode with the same amount of ressources.
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The pastix datat structure to initialize.
 *
 * @param[in] comm
 *          The MPI communicator associated to the pastix data structure
 *
 * @param[in] autosplit
 *          Enable automatic split when multiple processes run on the same node
 *
 *******************************************************************************/
static inline void
apiInitMPI( pastix_data_t *pastix,
            MPI_Comm       comm,
            int            autosplit )
{
    /*
     * Setup all communicators for autosplitmode and initialize number/rank of
     * processes.
     */
    pastix->pastix_comm = MPI_COMM_SELF; /* MPI is not available yet */
    MPI_Comm_size(pastix->pastix_comm, &(pastix->procnbr));
    MPI_Comm_rank(pastix->pastix_comm, &(pastix->procnum));
    assert(pastix->procnbr==1); /* MPI is not available yet */

#if defined(PASTIX_WITH_MPI)
    if ( autosplit )
    {
        int     i, len;
        char    procname[MPI_MAX_PROCESSOR_NAME];
        int     rc, key = pastix->procnum;
        int64_t color;
        (void)rc;

        /*
         * Get hostname to generate a hash that will be the color of each node
         * MPI_Get_processor_name is not used as it can returned different
         * strings for processes of a same physical node.
         */
        rc = gethostname(procname, MPI_MAX_PROCESSOR_NAME-1);
        assert(rc == 0);
        procname[MPI_MAX_PROCESSOR_NAME-1] = '\0';
        len = strlen( procname );

        /* Compute hash */
        color = 0;
        for (i = 0; i < len; i++) {
            color = color*256*sizeof(char) + procname[i];
        }

        /* Create intra-node communicator */
        MPI_Comm_split(pastix->pastix_comm, color, key, &(pastix->intra_node_comm));
        MPI_Comm_size(pastix->intra_node_comm, &(pastix->intra_node_procnbr));
        MPI_Comm_rank(pastix->intra_node_comm, &(pastix->intra_node_procnum));

        /* Create inter-node communicator */
        MPI_Comm_split(pastix->pastix_comm, pastix->intra_node_procnum, key, &(pastix->inter_node_comm));
        MPI_Comm_size(pastix->inter_node_comm, &(pastix->inter_node_procnbr));
        MPI_Comm_rank(pastix->inter_node_comm, &(pastix->inter_node_procnum));
    }
    else
#endif
    {
        pastix->intra_node_comm    = MPI_COMM_SELF;
        pastix->intra_node_procnbr = 1;
        pastix->intra_node_procnum = 0;
        pastix->inter_node_comm    = pastix->pastix_comm;
        pastix->inter_node_procnbr = pastix->procnbr;
        pastix->inter_node_procnum = pastix->procnum;
    }
    (void)autosplit;
    (void)comm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the solver instance with a bintab array to specify the
 * thread binding.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 * @param[in] pastix_comm
 *          The MPI communicator.
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 * @param[in]   bindtab
 *          Integer array of size iparm[IPARM_THREAD_NBR] that will specify the
 *          thread binding. NULL if let to the system.
 *
 *******************************************************************************/
void
pastixInitWithAffinity( pastix_data_t **pastix_data,
                        MPI_Comm        pastix_comm,
                        pastix_int_t   *iparm,
                        double         *dparm,
                        const int      *bindtab )
{
    pastix_data_t *pastix;

    /*
     * Allocate pastix_data structure when we enter PaStiX for the first time.
     */
    MALLOC_INTERN(pastix, 1, pastix_data_t);
    memset( pastix, 0, sizeof(pastix_data_t) );

    /*
     * Check if MPI is initialized
     */
    pastix->initmpi = 0;
#if defined(PASTIX_WITH_MPI)
    {
        int provided = MPI_THREAD_SINGLE;
        int flag = 0;
        MPI_Initialized(&flag);
        if ( !flag ) {
            MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided );
            pastix->initmpi = 1;
        }
        else {
            MPI_Query_thread( &provided );
        }
    }
#endif

    /*
     * Initialize iparm/dparm vectors and set them to default values if not set
     * by the user.
     */
    if ( iparm[IPARM_MODIFY_PARAMETER] == 0 ) {
        pastixInitParam( iparm, dparm );
    }

    pastix->iparm = iparm;
    pastix->dparm = dparm;

    pastix->steps = 0;

    pastix->isched = NULL;
#if defined(PASTIX_WITH_PARSEC)
    pastix->parsec = NULL;
#endif
#if defined(PASTIX_WITH_STARPU)
    pastix->starpu = NULL;
#endif

    apiInitMPI( pastix, pastix_comm, iparm[IPARM_AUTOSPLIT_COMM] );

    if ( (pastix->intra_node_procnbr > 1) &&
         (pastix->iparm[IPARM_THREAD_NBR] != -1 ) ) {
        pastix_print( pastix->procnum, 0,
                      "WARNING: Thread number forced by MPI autosplit feature\n" );
        iparm[IPARM_THREAD_NBR] = pastix->intra_node_procnbr;
    }

#if defined(PASTIX_GENERATE_MODEL)
    pastix_print( pastix->procnum, 0,
                  "WARNING: PaStiX compiled with -DPASTIX_GENERATE_MODEL forces single thread computations\n" );
    iparm[IPARM_THREAD_NBR] = 1;
#endif

    /*
     * Start the internal threads
     */
    pastix->isched = ischedInit( pastix->iparm[IPARM_THREAD_NBR], bindtab );
    pastix->iparm[IPARM_THREAD_NBR] = pastix->isched->world_size;

    /*
     * Start PaRSEC if compiled with it and scheduler set to PaRSEC
     */
#if defined(PASTIX_WITH_PARSEC)
    if ( (pastix->parsec == NULL) &&
         (iparm[IPARM_SCHEDULER] == PastixSchedParsec) ) {
        int argc = 0;
        pastix_parsec_init( pastix, &argc, NULL, bindtab );
    }
#endif /* defined(PASTIX_WITH_PARSEC) */

    /*
     * Start StarPU if compiled with it and scheduler set to StarPU
     */
#if defined(PASTIX_WITH_STARPU)
    if ( (pastix->starpu == NULL) &&
         (iparm[IPARM_SCHEDULER] == PastixSchedStarPU) ) {
        int argc = 0;
        pastix_starpu_init( pastix, &argc, NULL, bindtab );
    }
#endif /* defined(PASTIX_WITH_STARPU) */

    pastix->graph      = NULL;
    pastix->schur_n    = 0;
    pastix->schur_list = NULL;
    pastix->zeros_n    = 0;
    pastix->zeros_list = NULL;
    pastix->ordemesh   = NULL;

    pastix->symbmtx    = NULL;

    pastix->bcsc       = NULL;
    pastix->solvmatr   = NULL;

    pastix->cpu_models = NULL;
    pastix->gpu_models = NULL;

    pastix->dirtemp    = NULL;

    pastixModelsLoad( pastix );

    /* DIRTY Initialization for Scotch */
    srand(1);

    /* Environement variables */
    /* On Mac set VECLIB_MAXIMUM_THREADS if not setted */
    setenv("VECLIB_MAXIMUM_THREADS", "1", 0);

    if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
        pastixWelcome( pastix );
    }

    /* Initialization step done, overwrite anything done before */
    pastix->steps = STEP_INIT;

    *pastix_data = pastix;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the solver instance
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 * @param[in] pastix_comm
 *          The MPI communicator.
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 *******************************************************************************/
void
pastixInit( pastix_data_t **pastix_data,
            MPI_Comm        pastix_comm,
            pastix_int_t   *iparm,
            double         *dparm )
{
    pastixInitWithAffinity( pastix_data, pastix_comm,
                            iparm, dparm, NULL );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Finalize the solver instance
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixFinalize( pastix_data_t **pastix_data )
{
    pastix_data_t *pastix = *pastix_data;

    pastixSummary( *pastix_data );

    ischedFinalize( pastix->isched );

    if ( pastix->graph != NULL )
    {
        graphExit( pastix->graph );
        memFree_null( pastix->graph );
    }

    if ( pastix->ordemesh != NULL )
    {
        pastixOrderExit( pastix->ordemesh );
        memFree_null( pastix->ordemesh );
    }

    if ( pastix->symbmtx != NULL )
    {
        pastixSymbolExit( pastix->symbmtx );
        memFree_null( pastix->symbmtx );
    }

    if ( pastix->solvmatr != NULL )
    {
        solverExit( pastix->solvmatr );
        memFree_null( pastix->solvmatr );
    }

    if ( pastix->bcsc != NULL )
    {
        bcscExit( pastix->bcsc );
        memFree_null( pastix->bcsc );
    }

    if (pastix->schur_list != NULL )
    {
        memFree_null( pastix->schur_list );
    }
#if defined(PASTIX_WITH_PARSEC)
    if (pastix->parsec != NULL) {
        pastix_parsec_finalize( pastix );
    }
#endif /* defined(PASTIX_WITH_PARSEC) */
#if defined(PASTIX_WITH_STARPU)
    if (pastix->starpu != NULL) {
        pastix_starpu_finalize( pastix );
    }
#endif /* defined(PASTIX_WITH_STARPU) */

#if defined(PASTIX_WITH_MPI)
    if ( pastix->initmpi ) {
        MPI_Finalize();
    }
#endif

    if ( pastix->cpu_models != NULL ) {
        pastixModelsFree( pastix->cpu_models );
        pastix->cpu_models = NULL;
    }
    if ( pastix->gpu_models != NULL ) {
        pastixModelsFree( pastix->gpu_models );
        pastix->gpu_models = NULL;
    }

    if ( pastix->dirtemp != NULL ) {
        free( pastix->dirtemp );
    }
    memFree_null(*pastix_data);
}
