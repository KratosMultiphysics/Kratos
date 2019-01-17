/**
 *
 * @file pastix_zparsec.h
 *
 * Pastix PaRSEC codelets header
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @precisions normal z -> z c d s
 *
 **/
#ifndef _pastix_zparsec_h_
#define _pastix_zparsec_h_

void
parsec_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

void
parsec_zpxtrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

void
parsec_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

void
parsec_zhetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

void
parsec_zsytrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

#endif /* _pastix_zparsec_h_ */
