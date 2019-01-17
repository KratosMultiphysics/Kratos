/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "parsec_internal.h"
#include "dsparse.h"
#include "data_dist/sparse-matrix/pastix_internal/pastix_internal.h"

#include "zcsc2cblk.h"

parsec_object_t*
dsparse_zcsc2cblk_New(sparse_matrix_desc_t *A)
{
    parsec_zcsc2cblk_object_t *parsec_zcsc2cblk = NULL;

    parsec_zcsc2cblk = parsec_zcsc2cblk_new(A, (parsec_data_collection_t *)A );

    /* dsparse_add2arena_tile(((parsec_zcsc2cblk_Url_object_t*)parsec_zcsc2cblk)->arenas[PARSEC_zcsc2cblk_Url_DEFAULT_ARENA],  */
    /*                        A->mb*A->nb*sizeof(parsec_complex64_t), */
    /*                        PARSEC_ARENA_ALIGNMENT_SSE, */
    /*                        MPI_DOUBLE_COMPLEX, A->mb); */

    return (parsec_object_t*)parsec_zcsc2cblk;
}

void
dsparse_zcsc2cblk_Destruct( parsec_object_t *o )
{
    /* parsec_zcsc2cblk_object_t *parsec_zcsc2cblk = NULL; */
    /* parsec_zcsc2cblk = (parsec_zcsc2cblk_object_t *)o; */

    /*dsparse_datatype_undefine_type( &(opotrf->arenas[PARSEC_zcsc2cblk_Url_DEFAULT_ARENA]->opaque_dtt) );*/

    PARSEC_INTERNAL_OBJECT_DESTRUCT(o);
}

int dsparse_zcsc2cblk( parsec_context_t *parsec, sparse_matrix_desc_t *A)
{
    parsec_object_t *parsec_zcsc2cblk = NULL;
    int info = 0;

    parsec_zcsc2cblk = dsparse_zcsc2cblk_New( A );

    if ( parsec_zcsc2cblk != NULL )
    {
        parsec_enqueue( parsec, (parsec_object_t*)parsec_zcsc2cblk);
        parsec_progress( parsec );
        dsparse_zcsc2cblk_Destruct( parsec_zcsc2cblk );
    }
    return info;
}
