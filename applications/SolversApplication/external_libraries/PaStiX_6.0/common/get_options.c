/**
 *
 * @file get_options.c
 *
 * @copyright 2006-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2018-07-16
 *
 */
#include "common.h"
#include <unistd.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse iparm keywords to return its associated index in the iparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          The iparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the iparm in the array.
 * @retval -1 if the string is not an iparm parameter.
 *
 *******************************************************************************/
static inline pastix_iparm_t
parse_iparm( const char *iparm )
{
    if(0 == strcasecmp("iparm_verbose",               iparm)) { return IPARM_VERBOSE; }
    if(0 == strcasecmp("iparm_io_strategy",           iparm)) { return IPARM_IO_STRATEGY; }

    if(0 == strcasecmp("iparm_produce_stats",         iparm)) { return IPARM_PRODUCE_STATS; }

    if(0 == strcasecmp("iparm_mc64",                  iparm)) { return IPARM_MC64; }

    if(0 == strcasecmp("iparm_ordering",              iparm)) { return IPARM_ORDERING; }
    if(0 == strcasecmp("iparm_ordering_default",      iparm)) { return IPARM_ORDERING_DEFAULT; }

    if(0 == strcasecmp("iparm_scotch_switch_level",   iparm)) { return IPARM_SCOTCH_SWITCH_LEVEL; }
    if(0 == strcasecmp("iparm_scotch_cmin",           iparm)) { return IPARM_SCOTCH_CMIN; }
    if(0 == strcasecmp("iparm_scotch_cmax",           iparm)) { return IPARM_SCOTCH_CMAX; }
    if(0 == strcasecmp("iparm_scotch_frat",           iparm)) { return IPARM_SCOTCH_FRAT; }

    if(0 == strcasecmp("iparm_metis_ctype",           iparm)) { return IPARM_METIS_CTYPE; }
    if(0 == strcasecmp("iparm_metis_rtype",           iparm)) { return IPARM_METIS_RTYPE; }
    if(0 == strcasecmp("iparm_metis_no2hop",          iparm)) { return IPARM_METIS_NO2HOP; }
    if(0 == strcasecmp("iparm_metis_nseps",           iparm)) { return IPARM_METIS_NSEPS; }
    if(0 == strcasecmp("iparm_metis_niter",           iparm)) { return IPARM_METIS_NITER; }
    if(0 == strcasecmp("iparm_metis_ufactor",         iparm)) { return IPARM_METIS_UFACTOR; }
    if(0 == strcasecmp("iparm_metis_compress",        iparm)) { return IPARM_METIS_COMPRESS; }
    if(0 == strcasecmp("iparm_metis_ccorder",         iparm)) { return IPARM_METIS_CCORDER; }
    if(0 == strcasecmp("iparm_metis_pfactor",         iparm)) { return IPARM_METIS_PFACTOR; }
    if(0 == strcasecmp("iparm_metis_seed",            iparm)) { return IPARM_METIS_SEED; }
    if(0 == strcasecmp("iparm_metis_dbglvl",          iparm)) { return IPARM_METIS_DBGLVL; }

    if(0 == strcasecmp("iparm_sf_kass",               iparm)) { return IPARM_SF_KASS; }
    if(0 == strcasecmp("iparm_amalgamation_lvlcblk",  iparm)) { return IPARM_AMALGAMATION_LVLCBLK; }
    if(0 == strcasecmp("iparm_amalgamation_lvlblas",  iparm)) { return IPARM_AMALGAMATION_LVLBLAS; }

    if(0 == strcasecmp("iparm_reordering_split",      iparm)) { return IPARM_REORDERING_SPLIT; }
    if(0 == strcasecmp("iparm_reordering_stop",       iparm)) { return IPARM_REORDERING_STOP; }

    if(0 == strcasecmp("iparm_min_blocksize",         iparm)) { return IPARM_MIN_BLOCKSIZE; }
    if(0 == strcasecmp("iparm_max_blocksize",         iparm)) { return IPARM_MAX_BLOCKSIZE; }
    if(0 == strcasecmp("iparm_tasks2d_level",         iparm)) { return IPARM_TASKS2D_LEVEL; }
    if(0 == strcasecmp("iparm_tasks2d_width",         iparm)) { return IPARM_TASKS2D_WIDTH; }
    if(0 == strcasecmp("iparm_abs",                   iparm)) { return IPARM_ABS; }
    if(0 == strcasecmp("iparm_allcand",               iparm)) { return IPARM_ALLCAND; }

    if(0 == strcasecmp("iparm_incomplete",            iparm)) { return IPARM_INCOMPLETE; }
    if(0 == strcasecmp("iparm_level_of_fill",         iparm)) { return IPARM_LEVEL_OF_FILL; }

    if(0 == strcasecmp("iparm_factorization",         iparm)) { return IPARM_FACTORIZATION; }
    if(0 == strcasecmp("iparm_free_cscuser",          iparm)) { return IPARM_FREE_CSCUSER; }
    if(0 == strcasecmp("iparm_schur_fact_mode",       iparm)) { return IPARM_SCHUR_FACT_MODE; }

    if(0 == strcasecmp("iparm_schur_solv_mode",       iparm)) { return IPARM_SCHUR_SOLV_MODE; }
    if(0 == strcasecmp("iparm_applyperm_ws",          iparm)) { return IPARM_APPLYPERM_WS; }

    if(0 == strcasecmp("iparm_refinement",            iparm)) { return IPARM_REFINEMENT; }
    if(0 == strcasecmp("iparm_itermax",               iparm)) { return IPARM_ITERMAX; }
    if(0 == strcasecmp("iparm_gmres_im",              iparm)) { return IPARM_GMRES_IM; }

    if(0 == strcasecmp("iparm_scheduler",             iparm)) { return IPARM_SCHEDULER; }
    if(0 == strcasecmp("iparm_thread_nbr",            iparm)) { return IPARM_THREAD_NBR; }
    if(0 == strcasecmp("iparm_autosplit_comm",        iparm)) { return IPARM_AUTOSPLIT_COMM; }

    if(0 == strcasecmp("iparm_gpu_nbr",               iparm)) { return IPARM_GPU_NBR; }
    if(0 == strcasecmp("iparm_gpu_memory_percentage", iparm)) { return IPARM_GPU_MEMORY_PERCENTAGE; }
    if(0 == strcasecmp("iparm_gpu_memory_block_size", iparm)) { return IPARM_GPU_MEMORY_BLOCK_SIZE; }

    if(0 == strcasecmp("iparm_compress_min_width",    iparm)) { return IPARM_COMPRESS_MIN_WIDTH; }
    if(0 == strcasecmp("iparm_compress_min_height",   iparm)) { return IPARM_COMPRESS_MIN_HEIGHT; }
    if(0 == strcasecmp("iparm_compress_when",         iparm)) { return IPARM_COMPRESS_WHEN; }
    if(0 == strcasecmp("iparm_compress_method",       iparm)) { return IPARM_COMPRESS_METHOD; }
    if(0 == strcasecmp("iparm_compress_ortho",        iparm)) { return IPARM_COMPRESS_ORTHO; }

    if(0 == strcasecmp("iparm_modify_parameter",      iparm)) { return IPARM_MODIFY_PARAMETER; }
    if(0 == strcasecmp("iparm_start_task",            iparm)) { return IPARM_START_TASK; }
    if(0 == strcasecmp("iparm_end_task",              iparm)) { return IPARM_END_TASK; }
    if(0 == strcasecmp("iparm_mtx_type",              iparm)) { return IPARM_MTX_TYPE; }
    if(0 == strcasecmp("iparm_dof_nbr",               iparm)) { return IPARM_DOF_NBR; }

    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse dparm keywords to return its associated index in the dparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] dparm
 *          The dparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the dparm in the array.
 * @retval -1 if the string is not a dparm parameter.
 *
 *******************************************************************************/
static inline pastix_dparm_t
parse_dparm( const char *dparm )
{
    if(0 == strcasecmp("dparm_epsilon_refinement", dparm)) { return DPARM_EPSILON_REFINEMENT; }
    if(0 == strcasecmp("dparm_epsilon_magn_ctrl",  dparm)) { return DPARM_EPSILON_MAGN_CTRL;  }
    if(0 == strcasecmp("dparm_compress_tolerance", dparm)) { return DPARM_COMPRESS_TOLERANCE; }
    if(0 == strcasecmp("dparm_compress_min_ratio", dparm)) { return DPARM_COMPRESS_MIN_RATIO; }
    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse enum values for iparm array, and return the enum value
 * associated to it.
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The enum name to convert to its value
 *
 *******************************************************************************
 *
 * @retval The value if the enum associated to the string
 * @retval -1 if the string is not an enum in the pastix API.
 *
 *******************************************************************************/
static inline int
parse_enums( const char *string )
{
    if(0 == strcasecmp("pastix_task_init", string))
    {
        return PastixTaskInit;
    }
    if( (0 == strcasecmp("pastixtaskordering", string)) ||
        (0 == strcasecmp("pastixtaskscotch",   string)) )
    {
        return PastixTaskOrdering;
    }
    if( (0 == strcasecmp("pastixtasksymbfact", string)) ||
        (0 == strcasecmp("pastixtaskfax",      string)) )
    {
        return PastixTaskSymbfact;
    }
    if( (0 == strcasecmp("pastixtaskanalyse", string)) ||
        (0 == strcasecmp("pastixtaskblend",   string)) )
    {
        return PastixTaskAnalyze;
    }
    if( (0 == strcasecmp("pastixtasknumfact", string)) ||
        (0 == strcasecmp("pastixtasksopalin", string)) )
    {
        return PastixTaskNumfact;
    }
    if( (0 == strcasecmp("pastixtasksolve",  string)) ||
        (0 == strcasecmp("pastixtaskupdown", string)) )
    {
        return PastixTaskSolve;
    }
    if( (0 == strcasecmp("pastixtaskrefine",     string)) ||
        (0 == strcasecmp("pastixtaskrefinement", string)) )
    {
        return PastixTaskRefine;
    }
    if(0 == strcasecmp("pastixtaskclean", string))
    {
        return PastixTaskClean;
    }

    if(0 == strcasecmp("pastixverbosenot", string)) { return PastixVerboseNot; }
    if(0 == strcasecmp("pastixverboseno",  string)) { return PastixVerboseNo;  }
    if(0 == strcasecmp("pastixverboseyes", string)) { return PastixVerboseYes; }

    if(0 == strcasecmp("pastixiono",        string)) { return PastixIONo;        }
    if(0 == strcasecmp("pastixioload",      string)) { return PastixIOLoad;      }
    if(0 == strcasecmp("pastixiosave",      string)) { return PastixIOSave;      }
    if(0 == strcasecmp("pastixioloadgraph", string)) { return PastixIOLoadGraph; }
    if(0 == strcasecmp("pastixiosavegraph", string)) { return PastixIOSaveGraph; }
    if(0 == strcasecmp("pastixioloadcsc",   string)) { return PastixIOLoadCSC;   }
    if(0 == strcasecmp("pastixiosavecsc",   string)) { return PastixIOSaveCSC;   }

    if(0 == strcasecmp("pastixfactmodelocal", string)) { return PastixFactModeLocal; }
    if(0 == strcasecmp("pastixfactmodeschur", string)) { return PastixFactModeSchur; }
    if(0 == strcasecmp("pastixfactmodeboth",  string)) { return PastixFactModeBoth;  }

    if(0 == strcasecmp("pastixsolvmodelocal",     string)) { return PastixSolvModeLocal;     }
    if(0 == strcasecmp("pastixsolvmodeinterface", string)) { return PastixSolvModeInterface; }
    if(0 == strcasecmp("pastixsolvmodeschur",     string)) { return PastixSolvModeSchur;     }

    if(0 == strcasecmp("pastixrefinegmres",    string)) { return PastixRefineGMRES;    }
    if(0 == strcasecmp("pastixrefinecg",       string)) { return PastixRefineCG;       }
    if(0 == strcasecmp("pastixrefinesr",       string)) { return PastixRefineSR;       }
    if(0 == strcasecmp("pastixrefinebicgstab", string)) { return PastixRefineBiCGSTAB; }

    if(0 == strcasecmp("pastixorderscotch",   string)) { return PastixOrderScotch;   }
    if(0 == strcasecmp("pastixordermetis",    string)) { return PastixOrderMetis;    }
    if(0 == strcasecmp("pastixorderpersonal", string)) { return PastixOrderPersonal; }
    if(0 == strcasecmp("pastixorderptscotch", string)) { return PastixOrderPtScotch; }
    if(0 == strcasecmp("pastixorderparmetis", string)) { return PastixOrderParMetis; }

    if(0 == strcasecmp("pastixfactllt",  string)) { return PastixFactLLT;  }
    if(0 == strcasecmp("pastixfactldlt", string)) { return PastixFactLDLT; }
    if(0 == strcasecmp("pastixfactlu",   string)) { return PastixFactLU;   }
    if(0 == strcasecmp("pastixfactllh",  string)) { return PastixFactLLH;  }
    if(0 == strcasecmp("pastixfactldlh", string)) { return PastixFactLDLH; }

    if(0 == strcasecmp("pastixfactpotrf", string)) { return PastixFactPOTRF; }
    if(0 == strcasecmp("pastixfactpxtrf", string)) { return PastixFactPXTRF; }
    if(0 == strcasecmp("pastixfactgetrf", string)) { return PastixFactGETRF; }
    if(0 == strcasecmp("pastixfactsytrf", string)) { return PastixFactSYTRF; }
    if(0 == strcasecmp("pastixfacthetrf", string)) { return PastixFactHETRF; }

    if(0 == strcasecmp("pastixgeneral",   string)) { return PastixGeneral;   }
    if(0 == strcasecmp("pastixhermitian", string)) { return PastixHermitian; }
    if(0 == strcasecmp("pastixsymmetric", string)) { return PastixSymmetric; }

    if(0 == strcasecmp("pastixrefinegmres",    string)) { return PastixRefineGMRES;    }
    if(0 == strcasecmp("pastixrefinecg",       string)) { return PastixRefineCG;       }
    if(0 == strcasecmp("pastixrefinesr",       string)) { return PastixRefineSR;       }
    if(0 == strcasecmp("pastixrefinebicgstab", string)) { return PastixRefineBiCGSTAB; }

    if(0 == strcasecmp("pastixschedsequential", string)) { return PastixSchedSequential; }
    if(0 == strcasecmp("pastixschedstatic",     string)) { return PastixSchedStatic;     }
    if(0 == strcasecmp("pastixscheddynamic",    string)) { return PastixSchedDynamic;    }
    if(0 == strcasecmp("pastixschedparsec",     string)) { return PastixSchedParsec;     }
    if(0 == strcasecmp("pastixschedstarpu",     string)) { return PastixSchedStarPU;     }

    if(0 == strcasecmp("pastixcompressnever",      string)) { return PastixCompressNever;      }
    if(0 == strcasecmp("pastixcompresswhenbegin",  string)) { return PastixCompressWhenBegin;  }
    if(0 == strcasecmp("pastixcompresswhenend",    string)) { return PastixCompressWhenEnd;    }
    if(0 == strcasecmp("pastixcompresswhenduring", string)) { return PastixCompressWhenDuring; }

    if(0 == strcasecmp("pastixcompressmethodsvd",  string)) { return PastixCompressMethodSVD;  }
    if(0 == strcasecmp("pastixcompressmethodrrqr", string)) { return PastixCompressMethodRRQR; }

    /* If the value is directly given without string */
    {
        int value = atoi(string);
        if (value == 0)
        {
            if(0 == strcmp("0", string)) {
                return 0;
            }
            else {
                return -1;
            }
        }
        else {
            return value;
        }
    }
}

/**
 * @brief Print default usage for PaStiX binaries
 */
static inline void
pastix_usage(void)
{
    fprintf(stderr,
            "Matrix input (mandatory):\n"
            " -0 --rsa          : RSA/Matrix Market Fortran driver (only real)\n"
            " -1 --hb           : Harwell Boeing C driver\n"
            " -2 --ijv          : IJV coordinate C driver\n"
            " -3 --mm           : Matrix Market C driver\n"
            " -4 --spm          : SPM Matrix driver\n"
            " -9 --lap          : Generate a Laplacian (5-points stencil)\n"
            " -x --xlap         : Generate an extended Laplacian (9-points stencil)\n"
            " -G --graph        : SCOTCH Graph file\n"
            "\n"
            "Architecture arguments:\n"
            " -t --threads      : Number of threads per node (default: -1 to use the number of cores available)\n"
            " -g --gpus         : Number of gpus per node (default: 0)\n"
            " -s --sched        : Set the default scheduler (default: 1)\n"
            "                     0: Sequential, 1: Static, 2: PaRSEC, 3: StarPU\n"
            "\n"
            "Optional arguments:\n"
            " -f --fact                     : Choose factorization method (default: LU)\n"
            "                                 0: Cholesky, 1: LDL^[th], 2: LU, 3:LL^t, 4:LDL^t\n"
            "                                 3 and 4 are for complex matrices only\n"
            " -c --check                    : Choose the level of check to perform (default: 1)\n"
            "                                 0: None, 1: Backward error, 2: Backward and forward errors\n"
            " -o --ord                      : Choose between ordering libraries (default: scotch)\n"
            "                                 scotch, ptscotch, metis, parmetis\n"
            " -i --iparm <IPARM_ID> <value> : set any given integer parameter\n"
            " -d --dparm <DPARM_ID> <value> : set any given floating parameter\n"
            "\n"
            " -v --verbose[=lvl] : extra verbose output\n"
            " -h --help          : this message\n"
            "\n"
            );
}

/**
 * @brief Define the options and their requirement used by PaStiX
 */
#define GETOPT_STRING "0:1:2:3:4:9:x:G:t:g:s:o:f:c:i:d:v::h"

#if defined(HAVE_GETOPT_LONG)
/**
 * @brief Define the long options when getopt_long is available
 */
static struct option long_options[] =
{
    {"rsa",         required_argument,  0, '0'},
    {"hb",          required_argument,  0, '1'},
    {"ijv",         required_argument,  0, '2'},
    {"mm",          required_argument,  0, '3'},
    {"spm",         required_argument,  0, '4'},
    {"lap",         required_argument,  0, '9'},
    {"xlap",        required_argument,  0, 'x'},
    {"graph",       required_argument,  0, 'G'},

    {"threads",     required_argument,  0, 't'},
    {"gpus",        required_argument,  0, 'g'},
    {"sched",       required_argument,  0, 's'},

    {"ord",         required_argument,  0, 'o'},
    {"fact",        required_argument,  0, 'f'},
    {"check",       required_argument,  0, 'c'},
    {"iparm",       required_argument,  0, 'i'},
    {"dparm",       required_argument,  0, 'd'},

    {"verbose",     optional_argument,  0, 'v'},
    {"help",        no_argument,        0, 'h'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

/**
 *******************************************************************************
 *
 * @ingroup pastix_examples
 *
 * @brief PaStiX helper function to read command line options in examples.
 *
 * This function takes the command line arguments, and read the given parameters
 * (integers and doubles), as well as the matrix filename and the driver to read
 * it.
 *
 *******************************************************************************
 *
 * @param[in] argc
 *          The number of input parameters
 *
 * @param[in] argv
 *          The NULL terminated list of parameters
 *
 * @param[inout] iparam
 *          The integer array of parameters.
 *          On entry, must be initialized to the default value with pastixInitParam(),
 *          On exit, is updated with any option that matches the pastix parameters.
 *
 * @param[inout] dparam
 *          The double array of parameters.
 *          On entry, must be initialized to the default value with pastixInitParam(),
 *          On exit, is updated with any option that matches the pastix parameters.
 *
 * @param[inout] check
 *          On exit, the value is updated by the value of the -c option.
 *
 * @param[inout] driver
 *          On exit, contains the driver type give as option. -1, if no driver
 *          is specified.
 *
 * @param[out] filename
 *          The allocated string of the filename given with the driver.
 *
 *******************************************************************************/
void
pastixGetOptions( int argc, char **argv,
                  pastix_int_t *iparam, double *dparam,
                  int *check, spm_driver_t *driver, char **filename )
{
    int c;
    (void)dparam;

    if (argc == 1) {
        pastix_usage(); exit(0);
    }

    *driver = -1;
    do
    {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long( argc, argv, GETOPT_STRING,
                         long_options, NULL );
#else
        c = getopt( argc, argv, GETOPT_STRING );
#endif  /* defined(HAVE_GETOPT_LONG) */

        switch(c)
        {
        case '0':
            fprintf(stderr, "RSA driver is no longer supported and is replaced by the HB driver\n");
            pastix_attr_fallthrough;

        case '1':
            *driver = SpmDriverHB;
            *filename = strdup( optarg );
            break;

        case '2':
            *driver = SpmDriverIJV;
            *filename = strdup( optarg );
            break;

        case '3':
            *driver = SpmDriverMM;
            *filename = strdup( optarg );
            break;

        case '4':
            *driver = SpmDriverSPM;
            *filename = strdup( optarg );
            break;

        case '9':
            *driver = SpmDriverLaplacian;
            *filename = strdup( optarg );
            break;

        case 'x':
            *driver = SpmDriverXLaplacian;
            *filename = strdup( optarg );
            break;

        case 'G':
            *driver = SpmDriverGraph;
            *filename = strdup( optarg );
            break;

        case 't': iparam[IPARM_THREAD_NBR] = atoi(optarg); break;
        case 'g': iparam[IPARM_GPU_NBR] = atoi(optarg); break;

        case 'o':
            if (strncasecmp(optarg, "scotch", 6) == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderScotch;
            }
            else if (strncasecmp(optarg, "metis", 5) == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderMetis;
            }
            else if (strncasecmp(optarg, "ptscotch", 8) == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderPtScotch;
            }
            else if (strncasecmp(optarg, "parmetis", 8) == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderParMetis;
            }
            else if (strncasecmp(optarg, "personal", 8) == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderPersonal;
            }
            else {
                fprintf(stderr, "\nInvalid value for ordering option: %s\n\n", optarg);
                goto unknown_option;
            }
            break;

        case 'f': {
            int factotype = atoi( optarg );
            if ( (factotype >= 0) && (factotype <= 4)){
                iparam[IPARM_FACTORIZATION] = factotype;
            }
            else {
                fprintf(stderr, "\nInvalid value for factorization option: %s\n\n", optarg);
                goto unknown_option;
            }
        }
            break;

        case 'c': {
            int checkvalue = atoi( optarg );
            if ( (checkvalue >= 0) && (checkvalue < 6) ) {
                if ( check != NULL ) {
                    *check = checkvalue;
                }
            }
            else {
                fprintf(stderr, "\nInvalid value for check option: %s\n\n", optarg);
                goto unknown_option;
            }
        }
            break;

        case 's': {
            int schedtype = atoi( optarg );
            if ( (schedtype >= 0) && (schedtype <= 3)){
                iparam[IPARM_SCHEDULER] = schedtype;
            }
            else {
                fprintf(stderr, "\nInvalid value for scheduler option: %s\n\n", optarg);
                goto unknown_option;
            }
        }
            break;

        case 'i':
        {
            pastix_iparm_t iparm_idx;
            int iparm_val;

            /* Get iparm index */
            iparm_idx = parse_iparm( optarg );
            if ( iparm_idx == (pastix_iparm_t)-1 ) {
                fprintf(stderr, "\n%s is not a correct iparm parameter\n\n", optarg );
                goto unknown_option;
            }

            /* Get iparm value */
            iparm_val = parse_enums( argv[optind] );
            if ( iparm_val == -1 ){
                fprintf(stderr, "\n%s is not a correct value for the iparm parameters\n\n", argv[optind] );
                goto unknown_option;
            }
            iparam[iparm_idx] = iparm_val;
        }
        break;

        case 'd':
        {
            pastix_dparm_t dparm_idx;
            double dparm_val;

            /* Get iparm index */
            dparm_idx = parse_dparm( optarg );
            if ( dparm_idx == (pastix_dparm_t)-1 ) {
                fprintf(stderr, "\n%s is not a correct dparm parameter\n\n", optarg );
                goto unknown_option;
            }

            /* Get iparm value */
            dparm_val = atof( argv[optind] );
            dparam[dparm_idx] = dparm_val;
        }
        break;

        case 'v':
            if(optarg)  iparam[IPARM_VERBOSE] = atoi(optarg);
            else        iparam[IPARM_VERBOSE] = 2;
            break;

        case 'h':
            pastix_usage(); exit(EXIT_FAILURE);

        case ':':
            fprintf(stderr, "\nOption %c is missing an argument\n\n", c );
            goto unknown_option;

        case '?': /* getopt_long already printed an error message. */
            pastix_usage(); exit(EXIT_FAILURE);
        default:
            break;
        }
    } while(-1 != c);

    return;

  unknown_option:
    pastix_usage(); exit(EXIT_FAILURE);
}
