/**
 *
 * @file models.h
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
 * @addtogroup pastix_models
 * @{
 *    This module contains all the function to load performance models that will
 *    be used afterwards by simulation, and runtime schedulers.
 *
 **/
#ifndef _models_h_
#define _models_h_

#include "kernels/kernels_trace.h"

/**
 * @brief Model structure to store the coefficients table and its name
 */
typedef struct pastix_model_s {
    char  *name;                                     /**< Name of the computational unit considered by the model */
    double coefficients[4][PastixKernelLvl1Nbr][8];  /**< Coefficients table of the model                        */
} pastix_model_t;

/**
 * @brief Return the time in s of factorization kernels using a single size parameter.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] N The size parameter
 * @return The estimated time in s.
 */
static inline double
modelsGetCost1Param( const double *coefs, pastix_int_t N )
{
    /* a3 * N^3 + a2 * N^2 + a1 * N + a0 */
    double time = ((coefs[3] * N + coefs[2]) * N + coefs[1]) * N + coefs[0];
    return (time < 0.) ? 0. : time;
}

/**
 * @brief Return the time in s of TRSM kernels using two size parameters.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] M The first size parameter (number of rows)
 * @param[in] N The second size parameter (number of columns, size of the triangular matrix)
 * @return The estimated time in s.
 */
static inline double
modelsGetCost2Param( const double *coefs, pastix_int_t M, pastix_int_t N )
{
    /* a5 * M*N^2 + a4 * M*N + a3 * N^2 + a2 * M + a1 * N + a0 */
    double time = ((coefs[5] * (double)M + coefs[3]) * (double)N + coefs[4] * (double)M + coefs[1]) * (double)N + coefs[2] * (double)M + coefs[0];
    return (time < 0.) ? 0. : time;
}

/**
 * @brief Return the time in s of GEMM kernels using three size parameters.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] M The number of rows in the C matrix
 * @param[in] N The number of columns in the C matrix
 * @param[in] K The third dimension of the matrix product
 * @return The estimated time in s.
 */
static inline double
modelsGetCost3Param( const double *coefs, pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* a7 * M * N * K + a6 * M * K + a5 * K * N + a4 * M * N + a3 * M + a2 * N + a1 * K + a0 */
    double time = (coefs[7] * (double)M * (double)N * (double)K +
                   coefs[6] * (double)M * (double)K +
                   coefs[5] * (double)K * (double)N +
                   coefs[4] * (double)M * (double)N +
                   coefs[3] * (double)M +
                   coefs[2] * (double)N +
                   coefs[1] * (double)K +
                   coefs[0]);
    return (time < 0.) ? 0. : time;
}

pastix_model_t *pastixModelsNew();
void pastixModelsFree( pastix_model_t *model );
void pastixModelsLoad( pastix_data_t *pastix_data );

#endif /* _models_h_ */

/**
 * @}
 */
