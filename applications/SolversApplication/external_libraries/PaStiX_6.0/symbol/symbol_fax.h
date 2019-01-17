/**
 *
 * @file symbol_fax.h
 *
 * PaStiX fax symbol structure routines issued from Scotch esmumps library.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup symbol_dev_fax
 * @{
 *
 */
#ifndef _symbol_fax_h_
#define _symbol_fax_h_
/**
 * @brief Prime number for hashing vertex numbers.
 */
#define SYMBOL_FAX_HASHPRIME        17

/**
 * @brief The chained column block structure.
 *
 * These blocks are chained in a single linked list
 * for block merge with blocks of left columns.
 *
 */
typedef struct symbol_faxtlok_s {
  pastix_int_t frownum; /**< First row index            */
  pastix_int_t lrownum; /**< Last row index (inclusive) */
  pastix_int_t fcblknm; /**< Facing column block        */
  pastix_int_t nextnum; /**< Index of next block        */
} SymbolFaxTlok;

/**
 * @}
 */
#endif /* _symbol_fax_h_ */
