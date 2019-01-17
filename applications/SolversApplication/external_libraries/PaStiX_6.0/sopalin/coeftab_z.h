/**
 * @file coeftab_z.h
 *
 * Precision dependent coeficient array header.
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 * @addtogroup coeftab
 * @{
 *
 **/
#ifndef _coeftab_z_h_
#define _coeftab_z_h_

/**
 *    @name PastixComplex64 compression/uncompression routines
 *    @{
 */
pastix_int_t coeftab_zcompress  ( SolverMatrix *solvmtx );
void         coeftab_zuncompress( SolverMatrix *solvmtx );
void         coeftab_zmemory    ( const SolverMatrix *solvmtx );

/**
 *    @}
 *    @name PastixComplex64 initialization routines
 *    @{
 */
/**
 *    @}
 *    @name PastixComplex64 Extraction routines
 *    @{
 */
void coeftab_zgetschur( const SolverMatrix *solvmtx,
                        pastix_complex64_t *S, pastix_int_t lds );

void coeftab_zgetdiag( const SolverMatrix *solvmtx,
                       pastix_complex64_t *D, pastix_int_t incD );

/**
 *    @}
 *    @name PastixComplex64 debug routines
 *    @{
 */
void coeftab_zdump    ( pastix_data_t      *pastix_data,
                        const SolverMatrix *solvmtx,
                        const char         *filename );
int  coeftab_zdiff    ( pastix_coefside_t   side,
                        const SolverMatrix *solvA,
                        SolverMatrix       *solvB );

/**
 *    @}
 */
#endif /* _coeftab_z_h_ */

/**
 * @}
 */
