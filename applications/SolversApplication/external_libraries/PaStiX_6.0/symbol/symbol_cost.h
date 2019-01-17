/**
 *
 * @file symbol_cost.h
 *
 * PaStiX cost functions for the symbol structure
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @addtogroup symbol_dev_cost
 * @{
 **/
#ifndef _symbol_cost_h_
#define _symbol_cost_h_

/**
 * @brief Cost functions to compute statistics on the symbolic structure
 */
typedef struct symbol_function_s {
    double (*diag     )(pastix_int_t);                             /**< Return a statistic based on the diagonal block */
    double (*trsm     )(pastix_int_t, pastix_int_t);               /**< Return a statistic based on the sum of all
                                                                        off-diagonal of each column-block              */
    double (*update   )(pastix_int_t, pastix_int_t);               /**< Return a statistic for a large accumulated
                                                                        update per column-block                        */
    double (*blkupdate)(pastix_int_t, pastix_int_t, pastix_int_t); /**< Return a statistic for each individual
                                                                        off-diagonal block                             */
} symbol_function_t;

/**
 * @brief array of pointer to the flops functions per factorization and arithmetic
 */
extern symbol_function_t flopstable[2][5];
/**
 * @brief array of pointer to the performance functions per factorization and arithmetic
 */
extern symbol_function_t perfstable[2][5];

/**
 * @}
 */
#endif /* _symbol_cost_h_ */
