/**
 *
 * @file models.c
 *
 * PaStiX performance models routines
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "models.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Convert a kernel string name found in a model file to its kernel Id
 *
 *******************************************************************************
 *
 * @param[in] kernelstr
 *            The kernel string name
 *
 * @param[out] nbcoef
 *            The number of coefficient that this kernel will use. Set to 0 on
 *            failure.
 *
 *******************************************************************************
 *
 * @retval The kernel Id on success
 * @retval -1 on failure
 *
 *******************************************************************************/
int
modelsGetKernelId( const char *kernelstr,
                   int        *nbcoef )
{
    if(0 == strcasecmp("getrf",  kernelstr)) { *nbcoef = 4; return PastixKernelGETRF; }
    if(0 == strcasecmp("hetrf",  kernelstr)) { *nbcoef = 4; return PastixKernelHETRF; }
    if(0 == strcasecmp("potrf",  kernelstr)) { *nbcoef = 4; return PastixKernelPOTRF; }
    if(0 == strcasecmp("pxtrf",  kernelstr)) { *nbcoef = 4; return PastixKernelPXTRF; }
    if(0 == strcasecmp("sytrf",  kernelstr)) { *nbcoef = 4; return PastixKernelSYTRF; }

    if(0 == strcasecmp("trsmcblk1d", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMCblk1d; }
    if(0 == strcasecmp("trsmcblk2d", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMCblk2d; }
    if(0 == strcasecmp("trsmcblklr", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMCblkLR; }

    if(0 == strcasecmp("trsmblok2d", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMBlok2d; }
    if(0 == strcasecmp("trsmbloklr", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMBlokLR; }

    if(0 == strcasecmp("gemmcblk1d1d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblk1d1d; }
    if(0 == strcasecmp("gemmcblk1d2d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblk1d2d; }
    if(0 == strcasecmp("gemmcblk2d2d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblk2d2d; }
    if(0 == strcasecmp("gemmcblkfrlr", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblkFRLR; }
    if(0 == strcasecmp("gemmcblklrlr", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblkLRLR; }

    if(0 == strcasecmp("gemmblok2d2d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMBlok2d2d; }
    if(0 == strcasecmp("gemmbloklrlr", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMBlokLRLR; }

    *nbcoef = 0;
    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Propagate a given model to all other similare cases to be sure
 * everything is initialized.
 *
 * The given model coefficients defined by the couple (arithm, kernelid) is
 * first extended to all the kernels of the same family in the same arithmetic,
 * and it is then propagated to the other arithmetic by applying a computation
 * ratio on the coefficients.
 *     - Single real costs 1
 *     - Double real costs 2
 *     - Single complex costs 3
 *     - Double complex costs 4
 *
 *******************************************************************************
 *
 * @param[inout] model
 *          The pointer to the allocated model to complete.
 *
 * @param[in] arithm
 *          The arithmetic of the initial coefficients to replicate.
 *
 * @param[in] kernelid
 *          The kernel Id of the initial coefficients to replicate.
 *
 *******************************************************************************/
void
modelsPropagate( pastix_model_t *model,
                 int arithm, pastix_ktype_t kernelid )
{
    double *coefs0 = model->coefficients[arithm][kernelid];
    double ratio;
    int a, k;
    int kstart = 0;
    int kend = -1;

    /* Look for loaded information about factorization kernels */
    if ( kernelid < PastixKernelSCALOCblk ) {
        for( k=PastixKernelGETRF; k<=PastixKernelSYTRF; k++) {
            if ( (k == (int)kernelid) || (model->coefficients[arithm][k][0] != 0xdeadbeef) ) {
                continue;
            }

            ratio = (( k == (int)PastixKernelGETRF ) ? 2. : 1. ) / (( kernelid == PastixKernelGETRF ) ? 2. : 1. );

            model->coefficients[arithm][k][0] =         coefs0[0];
            model->coefficients[arithm][k][1] =         coefs0[1];
            model->coefficients[arithm][k][2] = ratio * coefs0[2];
            model->coefficients[arithm][k][3] = ratio * coefs0[3];
        }

        for( a=0; a<4; a++) {
            if (a == arithm) {
                continue;
            }
            ratio = (0.5 * a + 0.5) / (0.5 * arithm + 0.5);

            for( k=PastixKernelGETRF; k<=PastixKernelSYTRF; k++) {
                if ( model->coefficients[a][k][0] != 0xdeadbeef ) {
                    continue;
                }

                model->coefficients[a][k][0] = ratio * model->coefficients[arithm][k][0];
                model->coefficients[a][k][1] = ratio * model->coefficients[arithm][k][1];
                model->coefficients[a][k][2] = ratio * model->coefficients[arithm][k][2];
                model->coefficients[a][k][3] = ratio * model->coefficients[arithm][k][3];
            }
        }
    }
    else if ( kernelid < PastixKernelTRSMCblk1d ) {
    }
    else if ( kernelid < PastixKernelGEMMCblk1d1d ) {
        kstart = PastixKernelTRSMCblk1d;
        kend   = PastixKernelTRSMBlok2d;
    }
    else {
        kstart = PastixKernelGEMMCblk1d1d;
        kend   = PastixKernelGEMMBlok2d2d;
    }

    /*
     * Propagate to other kernels of the same arithmetic
     */
    for( k=kstart; k<=kend; k++) {
        if ( (k == (int)kernelid) || (model->coefficients[arithm][k][0] != 0xdeadbeef) ) {
            continue;
        }

        model->coefficients[arithm][k][0] = coefs0[0];
        model->coefficients[arithm][k][1] = coefs0[1];
        model->coefficients[arithm][k][2] = coefs0[2];
        model->coefficients[arithm][k][3] = coefs0[3];
        model->coefficients[arithm][k][4] = coefs0[4];
        model->coefficients[arithm][k][5] = coefs0[5];
        model->coefficients[arithm][k][6] = coefs0[6];
        model->coefficients[arithm][k][7] = coefs0[7];
    }

    /*
     * Propagate to other arithmetics
     */
    for( a=0; a<4; a++) {
        if (a == arithm) {
            continue;
        }
        ratio = (0.5 * a + 0.5) / (0.5 * arithm + 0.5);

        for( k=kstart; k<=kend; k++) {
            if ( model->coefficients[a][k][0] != 0xdeadbeef ) {
                continue;
            }

            model->coefficients[a][k][0] = ratio * model->coefficients[arithm][k][0];
            model->coefficients[a][k][1] = ratio * model->coefficients[arithm][k][1];
            model->coefficients[a][k][2] = ratio * model->coefficients[arithm][k][2];
            model->coefficients[a][k][3] = ratio * model->coefficients[arithm][k][3];
            model->coefficients[a][k][4] = ratio * model->coefficients[arithm][k][4];
            model->coefficients[a][k][5] = ratio * model->coefficients[arithm][k][5];
            model->coefficients[a][k][6] = ratio * model->coefficients[arithm][k][6];
            model->coefficients[a][k][7] = ratio * model->coefficients[arithm][k][7];
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Initialize the given model with the file given in parameters.
 *
 *******************************************************************************
 *
 * @param[inout] model
 *          The pointer to the allocated model to initialize.
 *
 * @param[in] modelfilename
 *          The name of the file in which the coefficient values are stored.
 *
 *******************************************************************************
 *
 * @return 0 on success.
 * @return -1 on failure.
 *
 *******************************************************************************/
int
modelsRead( pastix_model_t *model,
            const char *modelfilename )
{
    FILE *f = pastix_fopen( modelfilename );
    char *str, *strcoef;
    char kernelstr[13];
    int  rc, arithm, nbcoef;
    size_t strsize = 256;
    pastix_ktype_t kernelid;
    double *coefs;

    if ( f == NULL ) {
        fprintf(stderr, "Can't open model file\n");
        return -1;
    }

    str = malloc( strsize * sizeof(char) );
    do {
        rc = getline( &str, &strsize, f );
        if ( rc == -1 ) {
            perror( "modelsRead(getline header)" );
            return -1;
        }
    }
    while( str[0] == '#' );

    /* Read the model name */
    if ( rc == -1 ) {
        perror( "modelsRead(getline model name)" );
        return -1;
    }
    model->name = strdup( str );

    /* Read the model values */
    while( getline( &str, &strsize, f ) != -1 ) {

        /* Skip commented lines */
        if ( str[0] == '#' ) {
            continue;
        }

        /* Read the arithmetic, and the kernel name */
        if ( sscanf( str, "%d;%12[a-z0-9];", &arithm, kernelstr ) != 2 ) {
            fprintf(stderr, "modelRead: %s - Error reading line (%s)\n", model->name, str );
            continue;
        }

        if ( (arithm < 0) || (arithm > 3) ) {
            fprintf(stderr, "modelRead: %s - Incorrect arithmetic %d in line:\n\t%s\n",
                    model->name, arithm, str );
            continue;
        }

        kernelid = modelsGetKernelId( kernelstr, &nbcoef );
        if ( (int)kernelid == -1 ) {
            fprintf(stderr, "modelRead: %s - Incorrect kernel type %s in line:\n\t%s\n",
                    model->name, kernelstr, str );
            continue;
        }

        /* Read the corrrect number of coefficients and store them */
        coefs = model->coefficients[arithm][kernelid];
        strcoef = str + 3 + strlen( kernelstr );

        switch ( nbcoef ) {
        case 4:
            if ( sscanf( strcoef, "%le;%le;%le;%le",
                         coefs, coefs+1, coefs+2, coefs+3 ) != 4 )
            {
                fprintf(stderr, "modelRead: %s - Pb reading the 4 coefficients in line:\n\t%s\n", model->name, str );
                continue;
            }
            break;
        case 6:
            if ( sscanf( strcoef, "%le;%le;%le;%le;%le;%le",
                         coefs,   coefs+1, coefs+2,
                         coefs+3, coefs+4, coefs+5 ) != 6 )
            {
                fprintf(stderr, "modelRead: %s - Pb reading the 6 coefficients in line:\n\t%s\n", model->name, str );
                continue;
            }
            break;
        case 8:
            if ( sscanf( strcoef, "%le;%le;%le;%le;%le;%le;%le;%le",
                         coefs,   coefs+1, coefs+2, coefs+3,
                         coefs+4, coefs+5, coefs+6, coefs+7 ) != 8 )
            {
                fprintf(stderr, "modelRead: %s - Pb reading the 8 coefficients in line:\n\t%s\n", model->name, str );
                continue;
            }
            break;
        default:
            ;
        }

        modelsPropagate( model, arithm, kernelid );
    }

    fclose(f);
    free(str);

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Initialize the CPU model with default values.
 *
 *******************************************************************************
 *
 * @param[inout] model
 *          The pointer to the allocated model to initialize.
 *
 *******************************************************************************
 *
 * @return 0 on success.
 *
 *******************************************************************************/
int
modelsInitDefaultCPU( pastix_model_t *model )
{
    int a = 1; /* Real double */
    int ktype;
    double *coefs;

    assert( model != NULL );

    /*
     * All coefficiensts given are for double real arithmetic
     */
    model->name = strdup("AMD Opteron 6180 - Intel MKL");

    /* POTRF */
    ktype = PastixKernelPOTRF;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] =  4.071507e-07;
    coefs[1] = -1.469893e-07;
    coefs[2] =  1.707006e-08;
    coefs[3] =  2.439599e-11;
    modelsPropagate( model, a, ktype );

    /* TRSM Cblk */
    ktype = PastixKernelTRSMCblk2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] = 3.255168e-06;
    coefs[1] = 3.976198e-08;
    coefs[2] = 0.;
    coefs[3] = 0.;
    coefs[4] = 0.;
    coefs[5] = 2.626177e-10;
    modelsPropagate( model, a, ktype );

    /* TRSM Blok */
    /*
     * We don't have a TRSM blok model for this old architecture, so we use the
     * TRSM Cblk
     */

    /* GEMM Cblk */
    ktype = PastixKernelGEMMCblk2d2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] =  1.216278e-06;
    coefs[1] =  0.;
    coefs[2] = -2.704179e-10;
    coefs[3] =  1.148989e-07;
    coefs[4] =  2.724804e-10;
    coefs[5] =  1.328900e-09;
    coefs[6] =  0.;
    coefs[7] =  2.429169e-10;
    modelsPropagate( model, a, ktype );

    /* GEMM Blok */
    ktype = PastixKernelGEMMBlok2d2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] = 0.0;
    coefs[1] = 0.0;
    coefs[2] = 0.0;
    coefs[3] = 0.0;
    coefs[4] = 0.0;
    coefs[5] = 0.0;
    coefs[6] = 0.0;
    coefs[7] = 2. / 24.e9;
    modelsPropagate( model, a, ktype );

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Initialize the GPU model with default values.
 *
 *******************************************************************************
 *
 * @param[inout] model
 *          The pointer to the allocated model to initialize.
 *
 *******************************************************************************
 *
 * @return 0 on success.
 *
 *******************************************************************************/
int
modelsInitDefaultGPU( pastix_model_t *model )
{
    int a = 1; /* Real double */
    int ktype;
    double *coefs;

    assert( model != NULL );

    /*
     * All coefficiensts given are for double real arithmetic
     */
    model->name = strdup("Nvidia K40 GK1108L - CUDA 8.0");

    /* TRSM Blok */
    ktype = PastixKernelTRSMBlok2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] = -3.16663635648446e-05;
    coefs[1] =  2.63809317549331e-06;
    coefs[2] =  5.86447245256688e-07;
    coefs[3] = -1.57859559108480e-09;
    coefs[4] = -4.74303242824929e-09;
    coefs[5] =  5.36284121953867e-12;
    modelsPropagate( model, a, ktype );

    /* GEMM Cblk */
    ktype = PastixKernelGEMMCblk2d2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] =  1.216278e-06;
    coefs[1] =  0.;
    coefs[2] = -2.704179e-10;
    coefs[3] =  1.148989e-07;
    coefs[4] =  2.724804e-10;
    coefs[5] =  1.328900e-09;
    coefs[6] =  0.;
    coefs[7] =  2.429169e-10;
    modelsPropagate( model, a, ktype );

    /* GEMM Blok */
    ktype = PastixKernelGEMMBlok2d2d;
    coefs = &(model->coefficients[a][ktype][0]);
    coefs[0] = 0.0;
    coefs[1] = 0.0;
    coefs[2] = 0.0;
    coefs[3] = 0.0;
    coefs[4] = 0.0;
    coefs[5] = 0.0;
    coefs[6] = 0.0;
    coefs[7] = 2. /  1.2e12;
    modelsPropagate( model, a, ktype );

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Create a new model data structure and initialize the values to their
 * default.
 *
 *******************************************************************************
 *
 * @return The pointer to the allocated and initialized data structure.
 *
 *******************************************************************************/
pastix_model_t *
pastixModelsNew()
{
    pastix_model_t *model = malloc(sizeof(pastix_model_t));

    int a, k;

    memset( model, 0, sizeof( pastix_model_t ) );

    for(a=0; a<4; a++) {
        for(k=0; k<PastixKernelLvl1Nbr; k++) {
            model->coefficients[a][k][0] = 0xdeadbeef;
        }
    }
    return model;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Free a model data structure.
 *
 *******************************************************************************
 *
 * @param[inout] model
 *          The model structure to free.
 *
 *******************************************************************************/
void
pastixModelsFree( pastix_model_t *model )
{
    if ( model != NULL ) {
        if ( model->name != NULL ) {
            free(model->name);
        }
        free(model);
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Load the performance models that will be used by the solver
 *
 * This function initializes the model coefficients with the values stored in
 * the files defined by the environment variables PASTIX_MODELS_CPU and
 * PASTIX_MODELS_GPU. If they are not defined, models are initialized with the
 * embedded default models.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure in which to store the CPU and GPU models.
 *
 *******************************************************************************/
void
pastixModelsLoad( pastix_data_t *pastix_data )
{
    char *filename = NULL;
    int rc = 0;

    /*
     * Get the model filename for the CPUs
     */
    pastix_data->cpu_models = pastixModelsNew();
    filename = pastix_getenv( "PASTIX_MODELS_CPU" );

    if ( filename == NULL ) {
        rc = modelsInitDefaultCPU( pastix_data->cpu_models );
    }
    else {
        rc = modelsRead( pastix_data->cpu_models,
                         filename );
        pastix_cleanenv( filename );
    }
    if ( rc == -1 ) {
        pastixModelsFree( pastix_data->cpu_models );
        pastix_data->cpu_models = NULL;
    }

    /*
     * Get the model filename for the GPUs
     */
    pastix_data->gpu_models = pastixModelsNew();
    filename = pastix_getenv( "PASTIX_MODELS_GPU" );

    if ( filename == NULL ) {
        rc = modelsInitDefaultGPU( pastix_data->gpu_models );
    }
    else {
        rc = modelsRead( pastix_data->gpu_models,
                         filename );
        pastix_cleanenv( filename );
    }
    if ( rc == -1 ) {
        pastixModelsFree( pastix_data->gpu_models );
        pastix_data->gpu_models = NULL;
    }
}
