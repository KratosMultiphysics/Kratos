/**
 *
 * @file starpu_sparse_matrix.c
 *
 * PaStiX sparse matrix descriptor for StarPU.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "solver.h"
#include "pastix_starpu.h"
#include <starpu_data.h>

void
starpu_pastix_filter_list( void *father_interface,
                           void *child_interface,
                           struct starpu_data_filter *f,
                           unsigned id,
                           unsigned nchunks )
{
    struct starpu_vector_interface *vector_father = (struct starpu_vector_interface *) father_interface;
    struct starpu_vector_interface *vector_child = (struct starpu_vector_interface *) child_interface;

    size_t *length_tab = (size_t *) f->filter_arg_ptr;
    size_t  elemsize = vector_father->elemsize;
    size_t  chunk_size = length_tab[id+1] - length_tab[id];

    STARPU_ASSERT_MSG(vector_father->id == STARPU_VECTOR_INTERFACE_ID, "%s can only be applied on a vector data", __func__);
    vector_child->id = vector_father->id;
    vector_child->nx = chunk_size;
    vector_child->elemsize = elemsize;

    if (vector_father->dev_handle)
    {
        size_t current_pos = length_tab[id] * elemsize;

        if (vector_father->ptr) {
            vector_child->ptr = vector_father->ptr + current_pos;
        }
        vector_child->offset = vector_father->offset + current_pos;
        vector_child->dev_handle = vector_father->dev_handle;
    }

    (void)nchunks;
}

/**
 *******************************************************************************
 *
 * @brief Generate the StarPU descriptor of the sparse matrix.
 *
 * This function creates the StarPU descriptor that will provide tha data
 * mapping and memory location to StarPU for the computation.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure that describes the sparse matrix for
 *          PaStiX.
 *
 * @param[in] typesize
 *          The memory size of the arithmetic used to store the matrix
 *          coefficients.
 *
 * @param[in] mtxtype
 *          The type of sparse matrix to describe.
 *          @arg PastixGeneral:   The sparse matrix is general.
 *          @arg PastixSymmetric: The sparse matrix is lower triangular symmetric.
 *          @arg PastixHermitian: The sparse matrix is lower triangular hermitian.
 *
 * @param[in] nodes
 *          The number of processes used to solve the problem.
 *
 * @param[in] myrank
 *          The rank of the calling process.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_init( SolverMatrix *solvmtx,
                           int typesize, int mtxtype,
                           int nodes, int myrank )
{
    pastix_int_t   cblknbr, cblkmin2d;
    size_t key1, key2;
    SolverCblk *cblk;
    SolverBlok *blok, *lblok;
    pastix_int_t n=0, cblknum;
    pastix_int_t nbcol, nbrow;
    size_t size;

    starpu_sparse_matrix_desc_t *spmtx = solvmtx->starpu_desc;
    if ( spmtx != NULL ) {
        starpu_sparse_matrix_destroy( spmtx );
    }
    else {
        spmtx = (starpu_sparse_matrix_desc_t*)malloc(sizeof(starpu_sparse_matrix_desc_t));
    }

    spmtx->typesze = typesize;
    spmtx->mtxtype = mtxtype;
    spmtx->solvmtx = solvmtx;
    spmtx->cblktab_handle = NULL;
    spmtx->d_blocktab     = NULL;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    key1      = 2 * cblknbr;

    /* Initialize 1D cblk handlers */
    cblk = spmtx->solvmtx->cblktab;
    for(cblknum = 0;
        cblknum < cblkmin2d;
        cblknum++, n++, cblk++ )
    {
        starpu_data_handle_t *handler = (starpu_data_handle_t*)(cblk->handler);
        nbrow = cblk->stride;
        nbcol = cblk_colnbr( cblk );

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            pastix_lrblock_t *LRblocks = cblk->fblokptr->LRblock;
            pastix_int_t nbbloks = cblk[1].fblokptr - cblk[0].fblokptr;

            starpu_vector_data_register( handler, STARPU_MAIN_RAM,
                                         (uintptr_t)(LRblocks), nbbloks * 2, sizeof(pastix_lrblock_t) );

            if ( mtxtype == PastixGeneral ) {
                starpu_vector_data_register( handler + 1, STARPU_MAIN_RAM,
                                             (uintptr_t)(LRblocks + 1), nbbloks * 2, sizeof(pastix_lrblock_t) );
            }
        }
        else {
            starpu_vector_data_register( handler, STARPU_MAIN_RAM,
                                         (uintptr_t)(cblk->lcoeftab), nbrow * nbcol, spmtx->typesze );

            if ( mtxtype == PastixGeneral ) {
                starpu_vector_data_register( handler + 1, STARPU_MAIN_RAM,
                                             (uintptr_t)(cblk->ucoeftab), nbrow * nbcol, spmtx->typesze );
            }
        }
    }

    /* Initialize 2D cblk handlers */
    if ( cblkmin2d < cblknbr ) {
        struct starpu_data_filter filter = {
            .filter_func = starpu_pastix_filter_list
        };
        starpu_cblk_t *cblkhandle;
        size_t        *sizetab = NULL;
        pastix_int_t   nchildren, sizenbr = 0;

        spmtx->cblktab_handle = (starpu_cblk_t*)malloc( (cblknbr-cblkmin2d) * sizeof(starpu_cblk_t) );

        cblk = spmtx->solvmtx->cblktab + cblkmin2d;
        cblkhandle = spmtx->cblktab_handle;

        for(cblknum = cblkmin2d, n = 0;
            cblknum < cblknbr;
            cblknum++, n++, cblk++, cblkhandle++ )
        {
            starpu_data_handle_t *handler = (starpu_data_handle_t*)(cblk->handler);

            nbrow = cblk->stride;
            nbcol = cblk_colnbr( cblk );

            if ( cblk->cblktype & CBLK_COMPRESSED ) {
                pastix_lrblock_t *LRblocks = cblk->fblokptr->LRblock;
                pastix_int_t nbbloks = cblk[1].fblokptr - cblk[0].fblokptr;

                starpu_vector_data_register( handler, STARPU_MAIN_RAM,
                                             (uintptr_t)LRblocks, nbbloks * 2, sizeof(pastix_lrblock_t) );

                if ( mtxtype == PastixGeneral ) {
                    starpu_vector_data_register( handler + 1, STARPU_MAIN_RAM,
                                                 (uintptr_t)(LRblocks + 1), nbbloks * 2, sizeof(pastix_lrblock_t) );
                }
            }
            else {
                starpu_vector_data_register( handler, STARPU_MAIN_RAM,
                                             (uintptr_t)(cblk->lcoeftab), nbrow * nbcol, spmtx->typesze );

                if ( mtxtype == PastixGeneral ) {
                    starpu_vector_data_register( handler + 1, STARPU_MAIN_RAM,
                                                 (uintptr_t)(cblk->ucoeftab), nbrow * nbcol, spmtx->typesze );
                }
            }

            if ( !(cblk->cblktype & CBLK_TASKS_2D) )
                continue;

            /* Let's build the sizetab array */
            blok  = cblk[0].fblokptr;
            lblok = cblk[1].fblokptr;

            if ( (lblok - blok) >= sizenbr ) {
                sizenbr = (lblok - blok) + 1;
                free( sizetab );
                sizetab = malloc( sizenbr * sizeof(size_t) );
            }
            nchildren  = 0;
            sizetab[0] = 0;

            if ( cblk->cblktype & CBLK_COMPRESSED ) {
                /*
                 * Diagonal block
                 */
                sizetab[nchildren+1] = 2;
                nchildren++;

                /*
                 * Off-diagonal blocks
                 */
                blok++;
                for( ; blok < lblok; blok++ )
                {
                    nbrow = 1;

                    while( (blok+1 < lblok) &&
                           (blok[0].fcblknm == blok[1].fcblknm) &&
                           (blok[0].lcblknm == blok[1].lcblknm) )
                    {
                        blok++;
                        nbrow ++;
                    }
                    size = nbrow * 2;

                    sizetab[nchildren+1] = sizetab[nchildren] + size;
                    nchildren++;
                }
            }
            else {
                /*
                 * Diagonal block
                 */
                size = blok_rownbr( blok ) * cblk_colnbr( cblk );
                sizetab[nchildren+1] = sizetab[nchildren] + size;
                nchildren++;

                /*
                 * Off-diagonal blocks
                 */
                blok++;
                for( ; blok < lblok; blok++ )
                {
                    nbrow = blok_rownbr( blok );

                    while( (blok+1 < lblok) &&
                           (blok[0].fcblknm == blok[1].fcblknm) &&
                           (blok[0].lcblknm == blok[1].lcblknm) )
                    {
                        blok++;
                        nbrow += blok_rownbr( blok );
                    }
                    size = nbrow * cblk_colnbr( cblk );

                    sizetab[nchildren+1] = sizetab[nchildren] + size;
                    nchildren++;
                }
            }
            filter.nchildren = nchildren;
            filter.filter_arg_ptr = sizetab;

            cblkhandle->handlenbr = nchildren;
            if ( mtxtype == PastixGeneral ) {
                cblkhandle->handletab = (starpu_data_handle_t*)malloc( 2 * nchildren * sizeof(starpu_data_handle_t) );

                starpu_data_partition_plan( cblk->handler[0],
                                            &filter, cblkhandle->handletab );

                starpu_data_partition_plan( cblk->handler[1],
                                            &filter, cblkhandle->handletab + nchildren );
            }
            else {
                cblkhandle->handletab = (starpu_data_handle_t*)malloc( nchildren * sizeof(starpu_data_handle_t) );

                starpu_data_partition_plan( cblk->handler[0],
                                            &filter, cblkhandle->handletab );
            }

            nchildren = 0;
            blok  = cblk[0].fblokptr;
            lblok = cblk[1].fblokptr;

            /*
             * Diagonal block
             */
            blok->handler[0] = cblkhandle->handletab[ nchildren ];
            if ( mtxtype == PastixGeneral ) {
                blok->handler[1] = cblkhandle->handletab[ cblkhandle->handlenbr + nchildren ];
            }
            else {
                blok->handler[1] = NULL;
            }
            nchildren++;

            /*
             * Off-diagonal blocks
             */
            blok++;
            for( ; blok < lblok; blok++ )
            {
                blok->handler[0] = cblkhandle->handletab[ nchildren ];
                if ( mtxtype == PastixGeneral ) {
                    blok->handler[1] = cblkhandle->handletab[ cblkhandle->handlenbr + nchildren ];
                }
                else {
                    blok->handler[1] = NULL;
                }
                nchildren++;

                while( (blok < lblok) &&
                       (blok[0].fcblknm == blok[1].fcblknm) &&
                       (blok[0].lcblknm == blok[1].lcblknm) )
                {
                    blok++;
                    blok->handler[0] = NULL;
                    blok->handler[1] = NULL;
                }
            }
        }

        if (sizetab != NULL) {
            free( sizetab );
        }
    }
    solvmtx->starpu_desc = spmtx;

    (void)key1; (void)key2;
    (void)nodes; (void)myrank;
}

/**
 *******************************************************************************
 *
 * @brief Submit asynchronous calls to retrieve the data on main memory.
 *
 *******************************************************************************
 *
 * @param[inout] spmtx
 *          The sparse matrix descriptor to retrieve on main memory.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_getoncpu( starpu_sparse_matrix_desc_t *spmtx )
{
    SolverCblk *cblk;
    pastix_int_t i;

    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        assert( cblk->handler[0] );
        starpu_data_acquire_cb( cblk->handler[0], STARPU_R,
                                (void (*)(void*))&starpu_data_release,
                                cblk->handler[0] );

        if ( spmtx->mtxtype == PastixGeneral ) {
            assert( cblk->handler[1] );
            starpu_data_acquire_cb( cblk->handler[1], STARPU_R,
                                    (void (*)(void*))&starpu_data_release,
                                    cblk->handler[1] );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Free the StarPU descriptor of the sparse matrix.
 *
 * This function destroys the StarPU descriptor, but do not free the matrix data
 * that are managed by PaStiX.
 *
 *******************************************************************************
 *
 * @param[inout] spmtx
 *          The descriptor to free.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *spmtx )
{
    starpu_cblk_t *cblkhandle;
    SolverCblk *cblk;
    pastix_int_t i, cblkmin2d;

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<cblkmin2d; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            starpu_data_unregister( cblk->handler[0] );

            if ( spmtx->mtxtype == PastixGeneral ) {
                starpu_data_unregister( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;
    }

    cblkhandle = spmtx->cblktab_handle;
    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++, cblkhandle++)
    {
        if ( cblk->cblktype & CBLK_TASKS_2D ) {
             if ( cblk->handler[0] ) {
                 starpu_data_partition_clean( cblk->handler[0],
                                              cblkhandle->handlenbr,
                                              cblkhandle->handletab );
                 if ( spmtx->mtxtype == PastixGeneral ) {
                     starpu_data_partition_clean( cblk->handler[1],
                                                  cblkhandle->handlenbr,
                                                  cblkhandle->handletab + cblkhandle->handlenbr);
                 }
                 free( cblkhandle->handletab );
             }
        }

        if ( cblk->handler[0] ) {
            starpu_data_unregister( cblk->handler[0] );
            if ( spmtx->mtxtype == PastixGeneral ) {
                starpu_data_unregister( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;
    }

    if ( spmtx->cblktab_handle != NULL ) {
        free( spmtx->cblktab_handle );
    }
}

/**
 *@}
 */
