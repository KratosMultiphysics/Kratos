/**
 *
 * @file symbol_kass.h
 *
 * PaStiX kass symbol factorization routines
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup symbol_dev_kass
 * @{
 *
 **/
#ifndef _symbol_kass_h_
#define _symbol_kass_h_

/**
 * @brief Kass blocked csr structure
 */
typedef struct kass_csr_s {
    pastix_int_t   n;
    pastix_int_t  *nnz;
    pastix_int_t **rows;
} kass_csr_t;

void         kass_csrInit(   pastix_int_t n,
                             kass_csr_t *csr );
void         kass_csrClean(  kass_csr_t *csr );
pastix_int_t kass_csrGetNNZ( const kass_csr_t *csr );
int          kass_csrGenPA(  const pastix_graph_t *graphA,
                             const pastix_int_t   *perm,
                             kass_csr_t *graphPA );
void         kass_csrCompact(kass_csr_t *csr );

void         kassBuildSymbol(kass_csr_t   *P,
                             pastix_int_t  cblknbr,
                             const pastix_int_t *rangtab,
                             symbol_matrix_t *symbmtx);
void         kassPatchSymbol( symbol_matrix_t *symbmtx );

pastix_int_t kassFactDirect(const kass_csr_t *graphA,
                                  pastix_int_t  cblknbr,
                            const pastix_int_t *rangtab,
                                  pastix_int_t *treetab,
                            kass_csr_t   *graphL);
pastix_int_t kassFactLevel( const kass_csr_t   *graphA,
                                  pastix_int_t  level,
                                  kass_csr_t   *graphL);

void amalgamate(double rat_cblk, double rat_blas,
                kass_csr_t    *graphL,
                pastix_int_t   nnzL,
                pastix_int_t   snodenbr,
                const pastix_int_t  *snodetab,
                pastix_int_t  *treetab,
                pastix_int_t  *cblknbr,
                pastix_int_t **newrangtab,
                pastix_int_t **newtreetab,
                pastix_int_t  *nodetab,
                MPI_Comm pastix_comm );

/**
 * @}
 */
#endif /* _symbol_kass_h_ */
