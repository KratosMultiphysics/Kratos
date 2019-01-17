/**
 *
 * @file get_options.c
 *
 * @copyright 2006-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2017-06-15
 *
 */
#include <spm_tests.h>
#include <unistd.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

/**
 * @brief Print default usage for SpM binaries
 */
static inline void
spm_usage(void)
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
            " -h --help          : this message\n"
            "\n"
            );
}

/**
 * @brief Define the options and their requirement used by SpM
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

    {"help",        no_argument,        0, 'h'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

/**
 *******************************************************************************
 *
 * @ingroup spm_examples
 *
 * @brief SpM helper function to read command line options in examples.
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
 * @param[out] driver
 *          On exit, contains the driver type give as option. -1, if no driver
 *          is specified.
 *
 * @param[out] filename
 *          The allocated string of the filename given with the driver.
 *
 *******************************************************************************/
void
spmGetOptions( int argc, char **argv,
               spm_driver_t *driver, char **filename )
{
    int c;

    if (argc == 1) {
        spm_usage();
        exit(0);
    }

    *driver = (spm_driver_t)-1;
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
#if defined(SPM_WITH_FORTRAN)
            *driver = SpmDriverRSA;
            *filename = strdup( optarg );
#else
            fprintf(stderr, "spmGetOptions: Please compile with SPM_WITH_FORTRAN option to enable RSA driver or use HB driver instead\n");
            goto unknown_option;
#endif
            break;

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

        case 'h':
            spm_usage();
            exit(EXIT_FAILURE);

        case ':':
            fprintf(stderr, "\nOption %c is missing an argument\n\n", c );
            goto unknown_option;

        case '?': /* getopt_long already printed an error message. */
            spm_usage();
            exit(EXIT_FAILURE);
        default:
            break;
        }
    } while(-1 != c);

    return;

  unknown_option:
    spm_usage();
    exit(EXIT_FAILURE);
}
