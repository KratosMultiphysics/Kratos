/**
 *
 * @file kernels_trace.h
 *
 * Wrappers to trace kernels.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup eztrace_dev
 * @{
 *
 **/

#ifndef _kernels_trace_h_
#define _kernels_trace_h_

#include "common.h"
#include "flops.h"

/**
 * @brief Main stop enum event for all the events in traces
 */
#define PastixKernelStop 0

/**
 * @brief List of the Level 0 events that may be traced in PaStiX
 *
 * This is only the high level steps.
 */
typedef enum pastix_ktype0_e {
    PastixKernelLvl0Facto,
    PastixKernelLvl0Solve,
    PastixKernelLvl0Diag,
    PastixKernelLvl0Nbr
} pastix_ktype0_t;

/**
 * @brief List of the Level 1 events that may be traced in PaStiX
 *
 * This is the main information that traces all the major kernels during the
 * factorization step.
 */
typedef enum pastix_ktype_e {
    PastixKernelGETRF,        /**< LU diagonal block kernel             */
    PastixKernelHETRF,        /**< LDLh diagonal block kernel           */
    PastixKernelPOTRF,        /**< Cholesky diagonal block kernel       */
    PastixKernelPXTRF,        /**< Complex LL^t diagonal block kernel   */
    PastixKernelSYTRF,        /**< LDLt diagonal block kernel           */
    PastixKernelSCALOCblk,    /**< Scaling out-of-place of a panel      */
    PastixKernelSCALOBlok,    /**< Scaling out-of-place of a block      */
    PastixKernelTRSMCblk1d,   /**< TRSM applied to a panel in 1d layout */
    PastixKernelTRSMCblk2d,   /**< TRSM applied to a panel in 2d layout */
    PastixKernelTRSMCblkLR,   /**< TRSM applied to a panel in low-rank  */
    PastixKernelTRSMBlok2d,   /**< TRSM applied to a block in 2d layout */
    PastixKernelTRSMBlokLR,   /**< TRSM applied to a block in low-rank  */
    PastixKernelGEMMCblk1d1d, /**< GEMM applied from a panel in 1d layout to a panel in 1d layout */
    PastixKernelGEMMCblk1d2d, /**< GEMM applied from a panel in 1d layout to a panel in 2d layout */
    PastixKernelGEMMCblk2d2d, /**< GEMM applied from a panel in 2d layout to a panel in 2d layout */
    PastixKernelGEMMCblkFRLR, /**< GEMM applied from a panel in full-rank to a panel in low-rank  */
    PastixKernelGEMMCblkLRLR, /**< GEMM applied from a panel in low-rank to a panel in low-rank   */
    PastixKernelGEMMBlok2d2d, /**< GEMM applied from a block in 2d layout to a block in 2d layout */
    PastixKernelGEMMBlokLRLR, /**< GEMM applied from a block in low-rank to a block in low-rank   */
    PastixKernelLvl1Nbr
} pastix_ktype_t;

/**
 * @brief List of the Level 2 events that may be traced in PaStiX
 *
 * This is the low-level information that traces all the individual calls to
 * blas/lapack routines in the code. It is used to compute the number of flops
 * in low-rank compression, and to distinguish the amount of flops spent in each
 * part of the low-rank updates.
 *
 */
typedef enum pastix_ktype2_e {

    /* General kernels: similar in low-rank and dense */
    PastixKernelLvl2GETRF,             /**< LU diagonal block kernel           */
    PastixKernelLvl2HETRF,             /**< LDLh diagonal block kernel         */
    PastixKernelLvl2POTRF,             /**< Cholesky diagonal block kernel     */
    PastixKernelLvl2PXTRF,             /**< Complex LL^t diagonal block kernel */
    PastixKernelLvl2SYTRF,             /**< LDLt diagonal block kernel         */

    /* Solve operations */
    PastixKernelLvl2_FR_TRSM,
    PastixKernelLvl2_LR_TRSM,

    /* Update operations */
    PastixKernelLvl2_FR_GEMM,

    /* Formation (and application) of A * B */
    PastixKernelLvl2_LR_FRFR2FR,
    PastixKernelLvl2_LR_FRLR2FR,
    PastixKernelLvl2_LR_LRFR2FR,
    PastixKernelLvl2_LR_LRLR2FR,
    PastixKernelLvl2_LR_FRFR2LR,
    PastixKernelLvl2_LR_FRLR2LR,
    PastixKernelLvl2_LR_LRFR2LR,
    PastixKernelLvl2_LR_LRLR2LR,
    PastixKernelLvl2_LR_FRFR2null,
    PastixKernelLvl2_LR_FRLR2null,
    PastixKernelLvl2_LR_LRFR2null,
    PastixKernelLvl2_LR_LRLR2null,

    /* Compression kernels */
    PastixKernelLvl2_LR_init_compress,
    PastixKernelLvl2_LR_add2C_uncompress,
    PastixKernelLvl2_LR_add2C_recompress,
    PastixKernelLvl2_LR_add2C_updateCfr,
    PastixKernelLvl2_LR_add2C_orthou,
    PastixKernelLvl2_LR_add2C_rradd_orthogonalize, /**<< CGS, partialQR or fullQR */
    PastixKernelLvl2_LR_add2C_rradd_recompression,
    PastixKernelLvl2_LR_add2C_rradd_computeNewU,

    PastixKernelLvl2Nbr
} pastix_ktype2_t;

/**
 * @brief Total number of kernel events
 */
#define PastixKernelsNbr (PastixKernelLvl0Nbr + PastixKernelLvl1Nbr + PastixKernelLvl2Nbr)

/**
 * @brief Global array to store the number of flops executed per kernel
 */
extern volatile double kernels_flops[PastixKernelLvl1Nbr];

/**
 * @brief Lock to accumulate flops
 */
extern pastix_atomic_lock_t lock_flops;

/**
 * @brief Overall number of flops
 */
extern double overall_flops;

#if defined(PASTIX_WITH_EZTRACE)

#include "eztrace_module/kernels_ev_codes.h"

/**
 * @brief Define the level traced by the EZTrace module
 */
extern int pastix_eztrace_level;

#else

static inline void kernel_trace_start_lvl0     ( pastix_ktype0_t ktype )  { (void)ktype; }
static inline void kernel_trace_stop_lvl0      ( double flops )           { (void)flops; }
static inline void kernel_trace_start_lvl2     ( pastix_ktype2_t ktype )  { (void)ktype; (void)index; }
static inline void kernel_trace_stop_lvl2      ( double flops )           { (void)flops; }
static inline void kernel_trace_stop_lvl2_rank ( double flops, int rank ) { (void)flops; (void)rank; }

#endif

#if defined(PASTIX_GENERATE_MODEL)

/**
 * @brief Structure to store information linked to a kernel in order to generate
 * the cost model
 */
typedef struct pastix_model_entry_s {
    pastix_ktype_t ktype; /**< The type of the kernel             */
    int m;                /**< The first diemension of the kernel */
    int n;                /**< The second dimension of the kernel if present, 0 otherwise */
    int k;                /**< The third dimension of the kernel, 0 otherwise             */
    double time;          /**< The time spent in the kernel (s)                           */
} pastix_model_entry_t;

extern pastix_model_entry_t *model_entries;     /**< Array to all entries                 */
extern volatile int32_t      model_entries_nbr; /**< Index of the last entry in the array */
extern int32_t               model_size;        /**< Size of the model_entries array      */

#endif

void   kernelsTraceStart( const pastix_data_t *pastix_data );
double kernelsTraceStop(  const pastix_data_t *pastix_data );

/**
 *******************************************************************************
 *
 * @brief Start the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PASTIX_EZTRACE_LEVEL is different from 1.
 *
 *******************************************************************************
 *
 * @return the starting time if PASTIX_GENERATE_MODEL is enabled, 0. otherwise.
 *
 *******************************************************************************/
static inline double
kernel_trace_start( pastix_ktype_t ktype )
{
    double time = 0.;

#if defined(PASTIX_WITH_EZTRACE)

    if (pastix_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_0( KERNELS_LVL1_CODE(ktype) );
    }

#endif

#if defined(PASTIX_GENERATE_MODEL)

    time = clockGet();

#endif

    (void)ktype;
    return time;
}

/**
 *******************************************************************************
 *
 * @brief Stop the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PASTIX_EZTRACE_LEVEL is different from 1.
 *
 * @param[in] m
 *          The m parameter of the kernel (used by xxTRF, TRSM, and GEMM)
 *
 * @param[in] n
 *          The n parameter of the kernel (used by TRSM, and GEMM)
 *
 * @param[in] k
 *          The k parameter of the kernel (used by GEMM)
 *
 * @param[in] flops
 *          The number of flops of the kernel
 *
 * @param[in] starttime
 *          The stating time of the kernel. Used only if PASTIX_GENERATE_MODEL
 *          is enabled.
 *
 *******************************************************************************/
static inline void
kernel_trace_stop( pastix_ktype_t ktype, int m, int n, int k, double flops, double starttime )
{

#if defined(PASTIX_WITH_EZTRACE)

    if (pastix_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_1( KERNELS_CODE( PastixKernelStop ), flops );
    }

#endif

#if defined(PASTIX_GENERATE_MODEL)

    {
        double  time = clockGet() - starttime;
        int32_t index = pastix_atomic_inc_32b( &model_entries_nbr );

        if ( index < model_size ) {
            model_entries[index].ktype = ktype;
            model_entries[index].m     = m;
            model_entries[index].n     = n;
            model_entries[index].k     = k;
            model_entries[index].time  = time;
        }
        else {
            fprintf(stderr, "WARNING: too many kernels to log %d\n", index);
        }
    }

#endif

    pastix_atomic_lock( &lock_flops );
    overall_flops += flops;
    pastix_atomic_unlock( &lock_flops );

    (void)ktype;
    (void)m;
    (void)n;
    (void)k;
    (void)flops;
    (void)starttime;
    return;
}

#endif /* _kernels_trace_h_ */
