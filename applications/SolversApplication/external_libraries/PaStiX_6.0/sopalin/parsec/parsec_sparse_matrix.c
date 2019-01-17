/**
 *
 * @file parsec_sparse_matrix.c
 *
 * PaStiX sparse matrix descriptor for parsec.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @addtogroup pastix_parsec
 * @{
 *
 **/
#include "common.h"
#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_distribution.h>
#if defined(PASTIX_CUDA_FERMI)
#include <parsec/devices/cuda/dev_cuda.h>
#endif

#include "solver.h"
#include "pastix_parsec.h"

/**
 *******************************************************************************
 *
 * @brief Compute the triplet (uplo, cblknum, bloknum) from the key.
 *
 * This function convert the unique key identifier of each piece of data, to its
 * original information.
 *
 *******************************************************************************
 *
 * @param[in] key
 *          The key identifier to decode.
 *
 * @param[in] solvmtx
 *          The solver matrix structure to access information about the blocks
 *          and cblks.
 *
 * @param[out] uplo
 *          On exit, uplo is 0 if the key corresponds to the lower part, or 1
 *          for the upper part.
 *
 * @param[out] cblknum
 *          On exit, the index of the cblk encoded in the key.
 *
 * @param[out] bloknum
 *          On exit, the index of the blok encoded in the key.
 *          bloknum = 0, if the piece of data is the full cblk for 1D kernels.
 *          bloknum > 0, if the piece of data is the (bloknum-1)^th block in the cblk.
 *
 ******************************************************************************/
static inline void
spm_data_key_to_value( parsec_data_key_t   key,
                       const SolverMatrix *solvmtx,
                       int                *uplo,
                       pastix_int_t       *cblknum,
                       pastix_int_t       *bloknum )
{
    parsec_data_key_t key2;
    pastix_int_t cblkmin2d, cblknbr;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    key2 = 2 * cblknbr;

    /* This is a block */
    if ( key >= key2 ) {
        pastix_int_t m, n, ld;

        key2 = key - key2;
        ld   = solvmtx->cblkmaxblk * 2;

        m = key2 % ld;
        n = key2 / ld;

        *uplo    = m % 2;
        *bloknum = m / 2;
        *cblknum = cblkmin2d + n;
    }
    /* This is a cblk */
    else {
        *uplo    = key % 2;
        *cblknum = key / 2;
        *bloknum = -1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the unique key from the triplet (uplo, cblknum, bloknum).
 *
 * This function convert the unique triplet in a unique key identifier for each
 * piece of data.
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] ...
 *          This function must receive three int:
 *           - uplo: 0 for the lower part, 1 for the upper part
 *           - cblknum: the index of the cblk
 *           - bloknum: =0 if this is the full cblk,
                        >0 for the (bloknum-1)^th block in the cblk.
 *
 *******************************************************************************
 *
 * @return The unique key identifier for the piece of data.
 *
 ******************************************************************************/
static uint32_t
parsec_sparse_matrix_data_key( parsec_data_collection_t *mat, ... )
{
    va_list ap;
    parsec_sparse_matrix_desc_t *spmtx = (parsec_sparse_matrix_desc_t*)mat;
    int uplo;
    pastix_int_t cblknum, bloknum;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    uplo = uplo ? 1 : 0;

    if ( bloknum == -1 ) {
        return cblknum * 2 + uplo;
    }
    else {
        pastix_int_t offset, ld, cblknbr;
        pastix_int_t cblkmin2d, n;

        cblknbr   = spmtx->solvmtx->cblknbr;
        cblkmin2d = spmtx->solvmtx->cblkmin2d;
        ld        = spmtx->solvmtx->cblkmaxblk * 2;
        offset    = cblknbr * 2;
        n         = cblknum - cblkmin2d;

        return offset + n * ld + bloknum * 2 + uplo;
    }
}

/**
 *******************************************************************************
 *
 * @brief Return the rank of the owner of the piece of data (uplo, cblknum,
 * bloknum).
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] ...
 *          This function must receive three int:
 *           - uplo: 0 for the lower part, 1 for the upper part
 *           - cblknum: the index of the cblk
 *           - bloknum: =0 if this is the full cblk,
 *                      >0 for the (bloknum-1)^th block in the cblk.
 *
 *******************************************************************************
 *
 * @return The rank index of the owner of the data.
 *
 ******************************************************************************/
static uint32_t
parsec_sparse_matrix_rank_of( parsec_data_collection_t *mat, ... )
{
    (void)mat;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Return the rank of the owner of the piece of data (key)
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] key
 *          The unique key idenifier of a piece of data.
 *
 *******************************************************************************
 *
 * @return The rank index of the owner of the data.
 *
 ******************************************************************************/
static uint32_t
parsec_sparse_matrix_rank_of_key( parsec_data_collection_t    *mat,
                                  parsec_data_key_t  key )
{
    (void)mat; (void)key;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Return the rank of the virtual process owner of the piece of data
 * (uplo, cblknum, bloknum).
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] ...
 *          This function must receive three int:
 *           - uplo: 0 for the lower part, 1 for the upper part
 *           - cblknum: the index of the cblk
 *           - bloknum: =0 if this is the full cblk,
 *                      >0 for the (bloknum-1)^th block in the cblk.
 *
 *******************************************************************************
 *
 * @return The rank index of the virtual process owner of the data.
 *
 ******************************************************************************/
static int32_t
parsec_sparse_matrix_vpid_of( parsec_data_collection_t *mat, ... )
{
    (void)mat;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Return the rank of the virtual process owner of the piece of data (key)
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] key
 *          The unique key idenifier of a piece of data.
 *
 *******************************************************************************
 *
 * @return The rank index of the virtual process owner of the data.
 *
 ******************************************************************************/
static int32_t
parsec_sparse_matrix_vpid_of_key( parsec_data_collection_t    *mat,
                                  parsec_data_key_t  key )
{
    (void)mat; (void)key;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Return the data handler associated to the piece of data (uplo,
 * cblknum, bloknum).
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] ...
 *          This function must receive three int:
 *           - uplo: 0 for the lower part, 1 for the upper part
 *           - cblknum: the index of the cblk
 *           - bloknum: =0 if this is the full cblk,
 *                      >0 for the (bloknum-1)^th block in the cblk.
 *
 *******************************************************************************
 *
 * @return The pointer to the data handler of the data.
 *
 ******************************************************************************/
static parsec_data_t *
parsec_sparse_matrix_data_of( parsec_data_collection_t *mat, ... )
{
    parsec_sparse_matrix_desc_t *spmtx = (parsec_sparse_matrix_desc_t*)mat;
    SolverCblk *cblk;
    va_list ap;
    int uplo;
    pastix_int_t cblknum, bloknum;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    uplo = uplo ? 1 : 0;

    cblk = spmtx->solvmtx->cblktab + cblknum;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        assert( cblk->handler[uplo] );
        return (parsec_data_t*)(cblk->handler[uplo]);
    }
    /* This is a blok */
    else {
        SolverBlok *blok = cblk->fblokptr + bloknum;

        assert( blok->handler[uplo] );
        return (parsec_data_t*)(blok->handler[uplo]);
    }
}

/**
 *******************************************************************************
 *
 * @brief Return the data handler associated to the piece of data (key).
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] key
 *          The unique key idenifier of a piece of data.
 *
 *******************************************************************************
 *
 * @return The pointer to the data handler of the data.
 *
 ******************************************************************************/
static parsec_data_t *
parsec_sparse_matrix_data_of_key( parsec_data_collection_t    *mat,
                           parsec_data_key_t  key )
{
    parsec_sparse_matrix_desc_t *spmtx = (parsec_sparse_matrix_desc_t*)mat;
    SolverMatrix *solvmtx = spmtx->solvmtx;
    SolverCblk *cblk;
    int uplo;
    pastix_int_t cblknum, bloknum;

    spm_data_key_to_value( key, solvmtx,
                           &uplo, &cblknum, &bloknum );

    cblk = solvmtx->cblktab + cblknum;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        assert( cblk->handler[uplo] );
        return (parsec_data_t*)(cblk->handler[uplo]);
    }
    /* This is a blok */
    else {
        SolverBlok *blok = cblk->fblokptr + bloknum;

        assert( blok->handler[uplo] );
        return (parsec_data_t*)(blok->handler[uplo]);
    }
}

#if defined(PARSEC_PROF_TRACE)
/**
 *******************************************************************************
 *
 * @brief Convert the uinque key identifier to a human readable string.
 *
 *******************************************************************************
 *
 * @param[in] mat
 *          The sparse matrix descriptor.
 *
 * @param[in] key
 *          The unique key idenifier of a piece of data.
 *
 * @param[inout] buffer
 *          The char buffer of size buffer_size that will receive the string
 *          describing the piece of data.
 *
 * @param[in] buffer_size
 *          The size of the buffer buffer. buffer_size > 0.
 *
 *******************************************************************************
 *
 * @return The number of characters printed (excluding the null byte used to end
 * output to strings) as returned by snprintf.
 *
 ******************************************************************************/
static int
parsec_sparse_matrix_key_to_string( parsec_data_collection_t *mat,
                             uint32_t key,
                             char *buffer, uint32_t buffer_size )
{
    parsec_sparse_matrix_desc_t *spmtx = (parsec_sparse_matrix_desc_t*)mat;
    int uplo;
    pastix_int_t cblknum, bloknum;
    int res;

    spm_data_key_to_value( key, spmtx->solvmtx,
                           &uplo, &cblknum, &bloknum );

    res = snprintf(buffer, buffer_size, "(%d, %ld, %ld)",
                   uplo, (long int)cblknum, (long int)bloknum);
    if (res < 0)
    {
        printf("error in key_to_string for tile (%d, %ld, %ld) key: %u\n",
               uplo, (long int)cblknum, (long int)bloknum, key);
    }
    return res;
}
#endif

#if defined(PASTIX_CUDA_FERMI)
void
parsec_sparse_matrix_init_fermi( parsec_sparse_matrix_desc_t *spmtx,
                          const SolverMatrix   *solvmtx )
{
    gpu_device_t* gpu_device;
    SolverBlok *blok;
    pastix_int_t i, b, bloknbr, ndevices;
    size_t size;
    int *tmp, *bloktab;

    ndevices = parsec_devices_enabled();
    if ( ndevices <= 2 )
        return;

    bloknbr = solvmtx->bloknbr;
    size = 2 * bloknbr * sizeof(int);

	/**
     * Initialize array on CPU
     */
    bloktab = (int*)malloc( size );
    tmp = bloktab;
    for (b=0, blok = solvmtx->bloktab;
         b < bloknbr;
         b++, blok++, tmp+=2)
    {
        tmp[0] = blok->frownum;
        tmp[1] = blok->lrownum;
    }

    ndevices -= 2;
    spmtx->d_blocktab = calloc(ndevices, sizeof(void*));

    fprintf(stderr, "ndevices = %ld\n", ndevices );
    for(i = 0; i < ndevices; i++) {
        if( NULL == (gpu_device = (gpu_device_t*)parsec_devices_get(i+2)) ) continue;

        fprintf(stderr, "cuda index = %d\n", gpu_device->cuda_index );
        cudaSetDevice( gpu_device->cuda_index );

        cudaMalloc( &(spmtx->d_blocktab[i]),
                    size );

        cudaMemcpy( spmtx->d_blocktab[i],
                    bloktab, size,
                    cudaMemcpyHostToDevice );
    }
    free(bloktab);
}

void
parsec_sparse_matrix_destroy_fermi( parsec_sparse_matrix_desc_t *spmtx )
{
    gpu_device_t* gpu_device;
    pastix_int_t i, ndevices;

    ndevices = parsec_devices_enabled();
    if ( ndevices <= 2 )
        return;

    ndevices -= 2;
    for(i = 0; i < ndevices; i++) {
        if( NULL == (gpu_device = (gpu_device_t*)parsec_devices_get(i+2)) ) continue;

        cudaSetDevice( gpu_device->cuda_index );
        cudaFree( spmtx->d_blocktab[i] );
    }

    free( spmtx->d_blocktab );
}
#endif /*defined(PASTIX_CUDA_FERMI)*/

/**
 *******************************************************************************
 *
 * @brief Generate the PaRSEC descriptor of the sparse matrix.
 *
 * This function creates the PaRSEC descriptor that will provide tha data
 * mapping and memory location to PaRSEC for the computation.
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
parsec_sparse_matrix_init( SolverMatrix *solvmtx,
                           int typesize, int mtxtype,
                           int nodes, int myrank )
{
    parsec_sparse_matrix_desc_t *spmtx = solvmtx->parsec_desc;
    parsec_data_collection_t *o;
    pastix_int_t   cblknbr, cblkmin2d, ld;
    parsec_data_key_t key1, key2;
    SolverCblk *cblk;
    SolverBlok *blok, *fblok, *lblok;
    pastix_int_t m=0, n=0, cblknum;
    size_t size, offset;
    char *ptrL, *ptrU;

    if ( spmtx != NULL ) {
        parsec_sparse_matrix_destroy( spmtx );
    }
    else {
        spmtx = (parsec_sparse_matrix_desc_t*)malloc(sizeof(parsec_sparse_matrix_desc_t));
    }

    o = (parsec_data_collection_t*)spmtx;
    parsec_data_collection_init( o, nodes, myrank );

    o->data_key      = parsec_sparse_matrix_data_key;
#if defined(PARSEC_PROF_TRACE)
    o->key_to_string = parsec_sparse_matrix_key_to_string;
#endif

    o->rank_of     = parsec_sparse_matrix_rank_of;
    o->rank_of_key = parsec_sparse_matrix_rank_of_key;
    o->vpid_of     = parsec_sparse_matrix_vpid_of;
    o->vpid_of_key = parsec_sparse_matrix_vpid_of_key;
    o->data_of     = parsec_sparse_matrix_data_of;
    o->data_of_key = parsec_sparse_matrix_data_of_key;

    spmtx->typesze = typesize;
    spmtx->mtxtype = mtxtype;
    spmtx->solvmtx = solvmtx;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    ld        = solvmtx->cblkmaxblk * 2;
    key1      = 2 * cblknbr;

    /* Initialize 1D cblk handlers */
    cblk = spmtx->solvmtx->cblktab;
    for(cblknum = 0;
        cblknum < cblkmin2d;
        cblknum++, n++, cblk++ )
    {
        parsec_data_t **handler = (parsec_data_t**)(cblk->handler);
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;

        parsec_data_create( handler,
                            o, cblknum * 2, cblk->lcoeftab, size );

        if ( mtxtype == PastixGeneral ) {
            parsec_data_create( handler+1,
                                o, cblknum * 2 + 1, cblk->ucoeftab, size );
        }
    }

    /* Initialize 2D cblk handlers */
    cblk = spmtx->solvmtx->cblktab + cblkmin2d;
    for(cblknum = cblkmin2d, n = 0;
        cblknum < cblknbr;
        cblknum++, n++, cblk++ )
    {
        parsec_data_t **handler = (parsec_data_t**)(cblk->handler);
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;

        parsec_data_create( handler,
                            o, cblknum * 2, cblk->lcoeftab, size );

        if ( mtxtype == PastixGeneral ) {
            parsec_data_create( handler+1,
                                o, cblknum * 2 + 1, cblk->ucoeftab, size );
        }

        if ( !(cblk->cblktype & CBLK_TASKS_2D) )
            continue;

        /*
         * Diagonal block
         */
        ptrL   = cblk->lcoeftab;
        ptrU   = cblk->ucoeftab;
        blok   = cblk->fblokptr;
        size   = blok_rownbr( blok ) * cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        offset = blok->coefind * (size_t)spmtx->typesze;
        key2   = n * ld;

        assert(offset == 0);
        parsec_data_create( (parsec_data_t**)&(blok->handler[0]),
                            o, key1 + key2,
                            ptrL + offset, size );

        if ( mtxtype == PastixGeneral ) {
            parsec_data_create( (parsec_data_t**)&(blok->handler[1]),
                                o, key1 + key2 + 1,
                                ptrU + offset, size );
        }
        else {
            blok->handler[1] = NULL;
        }

        /*
         * Off-diagonal blocks
         */
        blok++; key2 += 2;
        lblok = cblk[1].fblokptr;
        for( ; blok < lblok; blok++, key2+=2 )
        {
            fblok = blok;
            m = 0;
            size   = blok_rownbr( blok );
            offset = blok->coefind * (size_t)spmtx->typesze;

            while( (blok < lblok) &&
                   (blok[0].fcblknm == blok[1].fcblknm) &&
                   (blok[0].lcblknm == blok[1].lcblknm) )
            {
                blok++; m++;
                size += blok_rownbr( blok );
            }
            size *= cblk_colnbr( cblk )
                *  (size_t)spmtx->typesze;

            parsec_data_create( (parsec_data_t**)&(fblok->handler[0]),
                                &spmtx->super, key1 + key2,
                                ptrL + offset, size );

            if ( mtxtype == PastixGeneral ) {
                parsec_data_create( (parsec_data_t**)&(fblok->handler[1]),
                                    &spmtx->super, key1 + key2 + 1,
                                    ptrU + offset, size );
            }
            else {
                fblok->handler[1] = NULL;
            }

            key2 += m * 2;
        }
    }

#if defined(PASTIX_CUDA_FERMI)
    parsec_sparse_matrix_init_fermi( spmtx, solvmtx );
#endif

    solvmtx->parsec_desc = spmtx;
}

/**
 *******************************************************************************
 *
 * @brief Free the PaRSEC descriptor of the sparse matrix.
 *
 * This function destroys the PaRSEC descriptor, but do not free the matrix data
 * that are managed by PaStiX.
 *
 *******************************************************************************
 *
 * @param[inout] spmtx
 *          The descriptor to free.
 *
 ******************************************************************************/
void
parsec_sparse_matrix_destroy( parsec_sparse_matrix_desc_t *spmtx )
{
    SolverCblk *cblk;
    SolverBlok *blok;
    pastix_int_t i, cblkmin2d;

#if defined(PASTIX_CUDA_FERMI)
    parsec_sparse_matrix_destroy_fermi( spmtx );
#endif

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<cblkmin2d; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            parsec_data_destroy( cblk->handler[0] );

            if ( spmtx->mtxtype == PastixGeneral ) {
                parsec_data_destroy( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;
    }

    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            parsec_data_destroy( cblk->handler[0] );
            if ( spmtx->mtxtype == PastixGeneral ) {
                parsec_data_destroy( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;

        blok = cblk->fblokptr;
        while( blok < cblk[1].fblokptr )
        {
            if ( blok->handler[0] ) {
                parsec_data_destroy( blok->handler[0] );
                if ( spmtx->mtxtype == PastixGeneral ) {
                    parsec_data_destroy( blok->handler[1] );
                }
            }

            blok->handler[0] = NULL;
            blok->handler[1] = NULL;

            blok++;
        }
    }

    parsec_data_collection_destroy( (parsec_data_collection_t*)spmtx );
}

/**
 *@}
 */
