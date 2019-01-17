/**
 *
 * @file eztrace_convert_kernels.c
 *
 * Module to convert eztrace events
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 */
#include <eztrace_convert.h>
#include <stdio.h>
#include <strings.h>
#include <GTG.h>
#include "kernels_trace.h"

static inline double kernels_dmin( double a, double b) {
    return ( a < b ) ? a : b;
}

static inline double kernels_dmax( double a, double b) {
    return ( a > b ) ? a : b;
}

/**
 * @brief EZTrace kernels module structure
 */
struct eztrace_convert_module kernels_module;

/**
 * @brief Kernels basic information structure for traces
 */
typedef struct kernels_s {
    char        *name;   /**< Name of the kernel in the traces          */
    gtg_color_t  color;  /**< Default color of the kernel in the traces */
} kernels_t;

/**
 * @brief Table of kernels default names and colors
 */
static kernels_t kernels_properties[PastixKernelsNbr];

/**
 * @brief Initialize the kernels_properties table
 */
void
define_kernels_properties()
{
    kernels_t *kernels_lvl1, *kernels_lvl2;
    int i;

    for (i=0; i<PastixKernelsNbr; i++){
        kernels_properties[i] = (kernels_t) {"unknown", GTG_BLACK};
    }

    /* Level 1 kernels */
    kernels_lvl1 = kernels_properties + PastixKernelLvl0Nbr;
    kernels_lvl1[PastixKernelGETRF]     = (kernels_t) {"lvl1_getrf",      GTG_RED};
    kernels_lvl1[PastixKernelHETRF]     = (kernels_t) {"lvl1_hetrf",      GTG_RED};
    kernels_lvl1[PastixKernelPOTRF]     = (kernels_t) {"lvl1_potrf",      GTG_RED};
    kernels_lvl1[PastixKernelPXTRF]     = (kernels_t) {"lvl1_pxtrf",      GTG_RED};
    kernels_lvl1[PastixKernelSYTRF]     = (kernels_t) {"lvl1_sytrf",      GTG_RED};
    kernels_lvl1[PastixKernelSCALOCblk] = (kernels_t) {"lvl1_scalo_cblk", GTG_SEABLUE};
    kernels_lvl1[PastixKernelSCALOBlok] = (kernels_t) {"lvl1_scalo_blok", GTG_GREEN};

    kernels_lvl1[PastixKernelTRSMCblk1d  ] = (kernels_t) {"lvl1_trsm_cblk_1d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMCblk2d  ] = (kernels_t) {"lvl1_trsm_cblk_2d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMCblkLR  ] = (kernels_t) {"lvl1_trsm_cblk_lr", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMBlok2d  ] = (kernels_t) {"lvl1_trsm_blok_2d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMBlokLR  ] = (kernels_t) {"lvl1_trsm_blok_lr", GTG_BLUE};
    kernels_lvl1[PastixKernelGEMMCblk1d1d] = (kernels_t) {"lvl1_gemm_cblk_1d1d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblk1d2d] = (kernels_t) {"lvl1_gemm_cblk_1d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblk2d2d] = (kernels_t) {"lvl1_gemm_cblk_2d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblkFRLR] = (kernels_t) {"lvl1_gemm_cblk_frlr", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblkLRLR] = (kernels_t) {"lvl1_gemm_cblk_lrlr", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMBlok2d2d] = (kernels_t) {"lvl1_gemm_blok_2d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMBlokLRLR] = (kernels_t) {"lvl1_gemm_blok_lrlr", GTG_GREEN};

    /* Level 2 kernels - Factorization kernels */
    kernels_lvl2 = kernels_lvl1 + PastixKernelLvl1Nbr;
    kernels_lvl2[PastixKernelLvl2GETRF] = (kernels_t) {"lvl2_getrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2HETRF] = (kernels_t) {"lvl2_hetrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2POTRF] = (kernels_t) {"lvl2_potrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2PXTRF] = (kernels_t) {"lvl2_pxtrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2SYTRF] = (kernels_t) {"lvl2_sytrf", GTG_RED};

    /* Level 2 kernels - Solve operations */
    kernels_lvl2[PastixKernelLvl2_FR_TRSM] = (kernels_t) {"lvl2_fr_trsm", GTG_SEABLUE};
    kernels_lvl2[PastixKernelLvl2_LR_TRSM] = (kernels_t) {"lvl2_lr_trsm", GTG_SEABLUE   };

    /* Level 2 kernels - Update operations */
    kernels_lvl2[PastixKernelLvl2_FR_GEMM]      = (kernels_t) {"lvl2_fr_gemm",   GTG_GREEN};

    kernels_lvl2[PastixKernelLvl2_LR_FRFR2FR]   = (kernels_t) {"lvl2_frfr2fr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_FRLR2FR]   = (kernels_t) {"lvl2_frlr2fr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRFR2FR]   = (kernels_t) {"lvl2_lrfr2fr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRLR2FR]   = (kernels_t) {"lvl2_lrlr2fr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_FRFR2LR]   = (kernels_t) {"lvl2_frfr2lr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_FRLR2LR]   = (kernels_t) {"lvl2_frlr2lr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRFR2LR]   = (kernels_t) {"lvl2_lrfr2lr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRLR2LR]   = (kernels_t) {"lvl2_lrlr2lr",   GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_FRFR2null] = (kernels_t) {"lvl2_frfr2null", GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_FRLR2null] = (kernels_t) {"lvl2_frlr2null", GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRFR2null] = (kernels_t) {"lvl2_lrfr2null", GTG_GREEN};
    kernels_lvl2[PastixKernelLvl2_LR_LRLR2null] = (kernels_t) {"lvl2_lrlr2null", GTG_GREEN};

    /* Level 2 kernels - Compression kernels */
    kernels_lvl2[PastixKernelLvl2_LR_init_compress]       = (kernels_t) {"lvl2_init_compress", GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_uncompress]    = (kernels_t) {"lvl2_uncompress",    GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_recompress]    = (kernels_t) {"lvl2_recompress",    GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_updateCfr]     = (kernels_t) {"lvl2_updateCfr",     GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_orthou]        = (kernels_t) {"lvl2_orthou",        GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_rradd_orthogonalize] = (kernels_t) {"lvl2_orthogonalize", GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_rradd_recompression] = (kernels_t) {"lvl2_recompression", GTG_PURPLE};
    kernels_lvl2[PastixKernelLvl2_LR_add2C_rradd_computeNewU]   = (kernels_t) {"lvl2_computeNewQ",   GTG_PURPLE};
}

/**
 * @brief Enum of the statistic values kept per data
 */
typedef enum counter_e {
    CounterMin, /**< Minimal value of the counter */
    CounterMax, /**< Maximal value of the counter */
    CounterSum, /**< Total value of the counter   */
    CounterSum2 /**< Sum of the square of each value to compute standard deviation */
} counter_t;

/**
 * @brief Information structure associated to one kernel on one thread to
 * compute statistics.
 */
typedef struct kernels_info_s {
    int    nb;          /**< Number of calls to the kernel       */
    int    uselessnb;   /**< Number of useless compression calls */
    double uselesstime; /**< Number of useless compression calls */
    double flop[4];     /**< Flops counters for statistics       */
    double time[4];     /**< Timing counters for statistics      */
    double perf[4];     /**< Performance counters for statistics */
} kernels_info_t;

/**
 * @brief Information structure associated to each thread to store temporary
 * information between two events while computing statistics.
 */
typedef struct kernels_thread_info_s {
    struct thread_info_t *p_thread;   /**< Common EZTrace thread information                                 */
    double                time_start; /**< Starting time of the last non PastixKernelStop event              */
    int                   current_ev; /**< Id of the last non PastixKernelStop event                         */
    kernels_info_t       *kernels;    /**< Table of size PastixKernelsNbr to store statistics on each kernel */
} kernels_thread_info_t;

/**
 * @brief Macro to initialize and eventually register a thread_info data
 * @param p_thread The current thread
 * @param var The name of the output kernels_thread_info_t variable
 * @param stats Boolean to specify if we are in stats or convert mode
 */
#define INIT_KERNELS_THREAD_INFO(p_thread, var, stats)                  \
    kernels_thread_info_t *var = (kernels_thread_info_t *)              \
        ezt_hook_list_retrieve_data(&p_thread->hooks, KERNELS_EVENTS_ID); \
    if(!(var)) {                                                        \
        var = kernels_register_thread_hook(p_thread, stats);            \
    }

/**
 *******************************************************************************
 *
 * @brief Initialize data associated with each kernel for a thread
 *
 *******************************************************************************
 *
 * @param[in] p_thread
 *          The reference to the thread managed by eztrace
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************
 *
 * @return  The structure containing the thread and statistics for each kernel
 *
 *******************************************************************************/
static kernels_thread_info_t *
kernels_register_thread_hook( struct thread_info_t *p_thread,
                              int stats )
{
    kernels_thread_info_t *p_info = (kernels_thread_info_t*) malloc(sizeof(kernels_thread_info_t));

    p_info->p_thread = p_thread;

    if (stats == 1){
        kernels_info_t *kernels = calloc( PastixKernelsNbr, sizeof(kernels_info_t) );
        int i;

        for(i=0; i<PastixKernelsNbr; i++) {
            kernels[i].flop[CounterMin] = 1.e99;
            kernels[i].time[CounterMin] = 1.e99;
            kernels[i].perf[CounterMin] = 1.e99;
        }
        p_info->kernels    = kernels;
        p_info->current_ev = -1;
    }

    ezt_hook_list_add(&p_info->p_thread->hooks, p_info, KERNELS_EVENTS_ID);
    return p_info;
}

/**
 *******************************************************************************
 *
 * @brief Init the kernels module in case of conversion mode.
 *
 * Meta-information of the events such as name or color are initialized by this
 * function.
 * They are also added to the trace file.
 *
 *******************************************************************************
 *
 * @retval 0 The events where correcly initialized
 * @retval 1 The events where not correcly initialized
 *
 *******************************************************************************/
int
eztrace_convert_kernels_init()
{
    if (get_mode() == EZTRACE_CONVERT) {
        int k;

        define_kernels_properties();

        for (k=0; k<PastixKernelsNbr; k++) {
            addEntityValue( kernels_properties[k].name, "ST_Thread",
                            kernels_properties[k].name, kernels_properties[k].color );
        }
    }
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Handle associated to the beginning of an event
 *
 *******************************************************************************
 *
 * @param[in] event_id
 *          The kernel id as one of the pastix enum for which statistics will be
 *          updated.
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************/
static inline void
handle_start( int event_id, int stats )
{
    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info, stats);

    assert( event_id < PastixKernelsNbr );

    if (stats == 1)
    {
        assert( p_info->current_ev == -1 );

        p_info->kernels[event_id].nb ++;
        p_info->current_ev = event_id;
        p_info->time_start = CURRENT;
    }
    else {
        pushState(CURRENT, "ST_Thread", thread_id, kernels_properties[event_id].name);
    }
}

/**
 *******************************************************************************
 *
 * @brief Handle associated to the end of an event
 *
 *******************************************************************************
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************/
static inline void
handle_stop( int stats )
{
    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info, stats);

    if (stats == 1)
    {
        int eventid = p_info->current_ev;
        int shift = PastixKernelLvl0Nbr + PastixKernelLvl1Nbr;
        kernels_info_t *kernel = p_info->kernels + eventid;
        double flop, time, perf;

        assert( eventid >= 0 );
        p_info->current_ev = -1;

        eventid -= shift;
        if ( (eventid == PastixKernelLvl2_LR_init_compress)    ||
             (eventid == PastixKernelLvl2_LR_add2C_recompress) ||
             (eventid == PastixKernelLvl2_LR_add2C_rradd_recompression) )
        {
            int rank;
            GET_PARAM_PACKED_2(CUR_EV, flop, rank);
            if ( rank == -1 ) {
                kernel->uselessnb ++;
                kernel->uselesstime += CURRENT - p_info->time_start;
            }
        }
        else
        {
            GET_PARAM_PACKED_1(CUR_EV, flop);
        }
        flop = flop * 1.e-9;
        kernel->flop[CounterMin ]  = kernels_dmin( kernel->flop[CounterMin], flop );
        kernel->flop[CounterMax ]  = kernels_dmax( kernel->flop[CounterMax], flop );
        kernel->flop[CounterSum ] += flop;
        kernel->flop[CounterSum2] += flop * flop;

        time = CURRENT - p_info->time_start;
        assert( time >= 0.0 );
        kernel->time[CounterMin ]  = kernels_dmin( kernel->time[CounterMin], time );
        kernel->time[CounterMax ]  = kernels_dmax( kernel->time[CounterMax], time );
        kernel->time[CounterSum ] += time;
        kernel->time[CounterSum2] += time * time;

        perf = flop / (time * 1.e-3);
        kernel->perf[CounterMin ]  = kernels_dmin( kernel->perf[CounterMin], perf );
        kernel->perf[CounterMax ]  = kernels_dmax( kernel->perf[CounterMax], perf );
        kernel->perf[CounterSum ] += perf;
        kernel->perf[CounterSum2] += perf * perf;
    }
    else{
        popState(CURRENT, "ST_Thread", thread_id);
    }
}

/**
 *******************************************************************************
 *
 * @brief Fonction called by eztrace_convert() to handle a single event
 *
 *******************************************************************************
 *
 * @param[in] ev
 *          The event to be handled
 *
 *******************************************************************************
 *
 * @retval 0 The event was not handled by the kernels module
 * @retval 1 The event was correclty handled by the kernels module
 *
 *******************************************************************************/
int
handle_kernels_events(eztrace_event_t *ev)
{
    int event_id;

    if( !(CUR_TRACE->start) ||
        !IS_A_KERNELS_EV(ev) )
    {
        return 0;
    }

    event_id = KERNELS_GET_CODE( ev );
    switch (event_id) {
    case PastixKernelStop:
        handle_stop( 0 );
        break;
    default:
        handle_start( event_id-1, 0 );
        break;
    }
    return 1;
}

/**
 *******************************************************************************
 *
 * @brief Fonction called by eztrace_stats() to handle a single event
 *
 *******************************************************************************
 *
 * @param[in] ev
 *          The event to be handled
 *
 *******************************************************************************
 *
 * @retval 0 The event was not handled by the kernels module
 * @retval 1 The event was correclty handled by the kernels module
 *
 *******************************************************************************/
int
handle_kernels_stats(eztrace_event_t *ev)
{
    int event_id;

    if( !(CUR_TRACE->start) ||
        !IS_A_KERNELS_EV(ev) )
    {
        return 0;
    }

    event_id = KERNELS_GET_CODE( ev );
    switch (event_id) {
    case PastixKernelStop:
        handle_stop( 1 );
        break;
    default:
        handle_start( event_id-1, 1 );
        break;
    }
    return 1;
}

/**
 *******************************************************************************
 *
 * @brief Print one line of statistics associated to a kernels_info_t structure.
 *
 *******************************************************************************
 *
 * @param[in] title
 *          The title of the line
 *
 * @param[in] kernel
 *          The structure for which the statistics need to be printed
 *
 * @param[inout] acc
 *          The accumulator structure.
 *          If acc != NULL, on exit acc = merge( acc, kernel )
 *
 *******************************************************************************/
void
print_kernel_stats( const char           *title,
                    const kernels_info_t *kernel,
                    kernels_info_t       *acc )
{
    const double *data;

    printf( "      %-23s | %8d |", title, kernel->nb );

    data = kernel->flop;
    printf( " %e %8.3g %8.3g %8.3g +-%8.3g |",
            data[CounterSum], data[CounterMin], data[CounterMax],
            data[CounterSum] / (double)kernel->nb,
            sqrt( (data[CounterSum2] - (data[CounterSum] * data[CounterSum]) / (double)kernel->nb) / (double)kernel->nb ) );

    data = kernel->time;
    printf( " %e %8.3g %8.3g %8.3g +-%8.3g |",
            data[CounterSum], data[CounterMin], data[CounterMax],
            data[CounterSum] / (double)kernel->nb,
            sqrt( (data[CounterSum2] - (data[CounterSum] * data[CounterSum]) / (double)kernel->nb) / (double)kernel->nb ) );

    data = kernel->perf;
    printf( " %8.3g %8.3g %8.3g +-%8.3g |\n",
            data[CounterMin], data[CounterMax],
            data[CounterSum] / (double)kernel->nb,
            sqrt( (data[CounterSum2] - (data[CounterSum] * data[CounterSum]) / (double)kernel->nb) / (double)kernel->nb ) );

    if ( acc != NULL ) {
        acc->nb += kernel->nb;
        acc->uselessnb   += kernel->uselessnb;
        acc->uselesstime += kernel->uselesstime;

        acc->flop[CounterMin ]  = kernels_dmin( acc->flop[CounterMin], kernel->flop[CounterMin] );
        acc->flop[CounterMax ]  = kernels_dmax( acc->flop[CounterMax], kernel->flop[CounterMax] );
        acc->flop[CounterSum ] += kernel->flop[CounterSum ];
        acc->flop[CounterSum2] += kernel->flop[CounterSum2];

        acc->time[CounterMin ]  = kernels_dmin( acc->time[CounterMin], kernel->time[CounterMin] );
        acc->time[CounterMax ]  = kernels_dmax( acc->time[CounterMax], kernel->time[CounterMax] );
        acc->time[CounterSum ] += kernel->time[CounterSum ];
        acc->time[CounterSum2] += kernel->time[CounterSum2];

        acc->perf[CounterMin ]  = kernels_dmin( acc->perf[CounterMin], kernel->perf[CounterMin] );
        acc->perf[CounterMax ]  = kernels_dmax( acc->perf[CounterMax], kernel->perf[CounterMax] );
        acc->perf[CounterSum ] += kernel->perf[CounterSum ];
        acc->perf[CounterSum2] += kernel->perf[CounterSum2];
    }
}

/**
 * @brief Print the full summary statistics table
 */
void
print_kernels_stats()
{
    kernels_info_t total[PastixKernelsNbr];
    kernels_info_t final;
    int i, j, k;
    int main_header   = 1;
    int thread_header = 1;

    define_kernels_properties();
    memset( total, 0, PastixKernelsNbr * sizeof(kernels_info_t) );
    memset( &final, 0, sizeof(kernels_info_t) );

    for(i=0; i<PastixKernelsNbr; i++) {
        total[i].flop[CounterMin] = 1.e99;
        total[i].time[CounterMin] = 1.e99;
        total[i].perf[CounterMin] = 1.e99;
    }
    final.flop[CounterMin] = 1.e99;
    final.time[CounterMin] = 1.e99;
    final.perf[CounterMin] = 1.e99;

    /* Browse the list of processes */
    for (i=0; i<NB_TRACES; i++) {
        struct eztrace_container_t *p_process = GET_PROCESS_CONTAINER(i);

        /* For each process, browse the list of threads */
        for(j=0; j<(int)(p_process->nb_children); j++) {

            struct eztrace_container_t *thread_container;
            struct thread_info_t       *p_thread;
            kernels_info_t *kernel;

            thread_container = p_process->children[j];
            p_thread = (struct thread_info_t*)(thread_container->container_info);

            if( !p_thread ) {
                continue;
            }

            INIT_KERNELS_THREAD_INFO(p_thread, p_info, 1);

            thread_header = 1;
            kernel = p_info->kernels;
            for (k=0; k<PastixKernelsNbr; k++, kernel++)
            {
                /*
                 * Print kernel statistics on the current thread
                 */
                if (kernel->nb > 0)
                {
                    if ( main_header ) {
                        printf( "         %20s | Nb calls | %-50s | %-50s | %-37s |\n"
                                "         %20s |          | %12s %8s %8s %8s %10s |"
                                " %12s %8s %8s %8s %10s | %8s %8s %8s %10s |\n",
                                "", "Number of flop (GFlop)", "Timing (ms)", "Performance (GFlop/s)",
                                "", "total", "min", "max", "avg", "sd",
                                    "total", "min", "max", "avg", "sd",
                                             "min", "max", "avg", "sd" );
                        main_header = 0;
                    }

                    if ( thread_header ) {
                        printf( "  Thread %20s |\n", thread_container->name );
                        thread_header = 0;
                    }

                    print_kernel_stats( kernels_properties[k].name,
                                        kernel, &(total[k]) );
                }
            }
            if ( !thread_header ) {
                printf("\n");
            }
        }
    }

    /*
     * Print summary per kernel
     */
    thread_header = 1;
    for (k=0; k<PastixKernelsNbr; k++)
    {
        if ( total[k].nb > 0 ) {
            if ( thread_header ) {
                printf( "  Total  %20s |\n", "" );
                thread_header = 0;
            }

            print_kernel_stats( kernels_properties[k].name,
                                &(total[k]), &final );
        }
    }
    if ( !thread_header ) {
        printf("\n");
    }

    thread_header = 1;
    for (k=0; k<PastixKernelsNbr; k++)
    {
        if ( total[k].uselessnb > 0 ) {
            if ( thread_header ) {
                printf( "  %-27s |\n", "Useless calls" );
                thread_header = 0;
            }

            printf( "      %-23s | %8d | %e (%8.3g %%) |\n",
                    kernels_properties[k].name, total[k].uselessnb,
                    total[k].uselesstime, (total[k].uselesstime / final.time[CounterSum]) * 100. );
        }
    }
    if ( !thread_header ) {
        printf("\n");
    }

    /*
     * Print final summary
     */
    if ( final.nb > 0 ) {
        print_kernel_stats("Summary", &final, NULL );
        printf("\n");
        printf( "      %-23s | %8s | %8s | %8s |\n",
                "", "Calls", "Flops", "Time" );
        for (k=0; k<PastixKernelsNbr; k++)
        {
            if ( total[k].nb > 0 ) {
                printf( "      %-23s | %6.2g %% | %6.2g %% | %6.2g %% |\n",
                        kernels_properties[k].name,
                        ((double)total[k].nb / (double)final.nb) * 100.,
                        (total[k].flop[CounterSum] / final.flop[CounterSum]) * 100.,
                        (total[k].time[CounterSum] / final.time[CounterSum]) * 100. );
            }
        }
    }
}

/**
 * @brief Kernels module initialization function called by EZTrace to load
 * the module
 */
void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
    int rc;

    kernels_module.api_version   = EZTRACE_API_VERSION;
    kernels_module.init          = eztrace_convert_kernels_init;
    kernels_module.handle        = handle_kernels_events;
    kernels_module.handle_stats  = handle_kernels_stats;
    kernels_module.print_stats   = print_kernels_stats;
    kernels_module.module_prefix = KERNELS_EVENTS_ID;

    rc = asprintf(&kernels_module.name,        "kernels"       );
    rc = asprintf(&kernels_module.description, "PaStiX kernels");

    kernels_module.token.data = &kernels_module;
    eztrace_convert_register_module(&kernels_module);

    (void)rc;
}

/**
 * @brief Kernels module finalization function called by EZTrace to unload the
 * module
 */
void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
    printf("unloading module kernels\n");
}

