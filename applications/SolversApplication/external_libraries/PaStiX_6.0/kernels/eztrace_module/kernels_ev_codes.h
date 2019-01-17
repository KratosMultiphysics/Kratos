/**
 *
 * @file kernels_ev_codes.h
 *
 * Wrappers to trace kernels with eztrace
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 * @addtogroup eztrace_dev
 * @{
 *
 **/

#ifndef _kernels_ev_codes_h_
#define _kernels_ev_codes_h_

#include "eztrace.h"
#include "eztrace_sampling.h"
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x51)
#define KERNELS_PREFIX       GENERATE_USER_MODULE_PREFIX( KERNELS_EVENTS_ID )

#define KERNELS_PREFIX_MASK  (((1 << NB_BITS_PREFIX) -1) << NB_BITS_EVENTS)
#define KERNELS_EVENTS_MASK  (~KERNELS_PREFIX_MASK)

#define IS_A_KERNELS_EV(ev)  ((LITL_READ_GET_CODE(ev) & KERNELS_PREFIX_MASK) == KERNELS_PREFIX)
#define KERNELS_GET_CODE(ev)  (LITL_READ_GET_CODE(ev) & KERNELS_EVENTS_MASK)

#define KERNELS_CODE(event) ( KERNELS_PREFIX | (event) )

#define KERNELS_LVL0_CODE(event) ( KERNELS_CODE( event + 1 ) )
#define KERNELS_LVL1_CODE(event) ( KERNELS_CODE( event + PastixKernelLvl0Nbr + 1 ) )
#define KERNELS_LVL2_CODE(event) ( KERNELS_CODE( event + PastixKernelLvl1Nbr + PastixKernelLvl0Nbr + 1 ) )

extern int pastix_eztrace_level;

/**
 *******************************************************************************
 *
 * @brief Start to trace a kernel of level 0
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          The kernel's type to trace
 *
 *******************************************************************************/
static inline void
kernel_trace_start_lvl0( pastix_ktype0_t ktype )
{
    if (pastix_eztrace_level == 0){
        EZTRACE_EVENT_PACKED_0(KERNELS_LVL0_CODE(ktype));
    }
}

/**
 *******************************************************************************
 *
 * @brief Stop to trace a kernel of level 0
 *
 *******************************************************************************
 *
 * @param[in] flops
 *          The number of operations performed during the call
 *
 *******************************************************************************/
static inline void
kernel_trace_stop_lvl0( double flops )
{
    if (pastix_eztrace_level == 0){
        EZTRACE_EVENT_PACKED_1(KERNELS_CODE(PastixKernelStop), flops);
    }
}

/**
 *******************************************************************************
 *
 * @brief Start to trace a kernel of level 2
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          The kernel's type to trace
 *
 *******************************************************************************/
static inline void
kernel_trace_start_lvl2( pastix_ktype2_t ktype )
{
    if (pastix_eztrace_level == 2){
        EZTRACE_EVENT_PACKED_0( KERNELS_LVL2_CODE(ktype) );
    }
}

/**
 *******************************************************************************
 *
 * @brief Stop to trace a kernel of level 2
 *
 *******************************************************************************
 *
 * @param[in] flops
 *          The number of operations performed during the call
 *
 *******************************************************************************/
static inline void
kernel_trace_stop_lvl2( double flops )
{
    if (pastix_eztrace_level == 2){
        EZTRACE_EVENT_PACKED_1(KERNELS_CODE(PastixKernelStop), flops);
    }
}

static inline void
kernel_trace_stop_lvl2_rank( double flops, int rank )
{
    if (pastix_eztrace_level == 2){
        EZTRACE_EVENT_PACKED_2(KERNELS_CODE(PastixKernelStop), flops, rank);
    }
}

#endif /* _kernels_ev_codes_h_ */
