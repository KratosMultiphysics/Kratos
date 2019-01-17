/**
 *
 * @file solver_copy.c
 *
 * PaStiX solver matrix copy and reallocation functions.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "queue.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver_null
 *
 * @brief Copy the solver matrix data structure from solvin to solvout.
 *
 * Every data is copied, event the coefficient if they are allocated and
 * initialized.
 * It is also used to reallocate the data in a contiguous way after the
 * initialization that allocates all internal arrays in multiple step which
 * might results in fragmentation.
 * @warning This function is not able to copy a solver matrix with low rank
 * blocks yet.
 *
 *******************************************************************************
 *
 * @param[in] solvin
 *          The solver matrix structure to duplicate.
 *
 * @param[out] solvout
 *          The allocated pointer to the solver matrix structure that will
 *          contain the copy.
 *
 * @param[in] flttype
 *          The floating point arithmetic sued in the input solver matrix to
 *          know the size of the memory space to duplicate for the coefficients.
 *
 *******************************************************************************/
static inline void
solver_copy( const SolverMatrix *solvin,
             SolverMatrix       *solvout,
             int                 flttype )
{
    SolverCblk *solvcblk;
    SolverBlok *solvblok;
    pastix_int_t i;

    /* Copy tasktab */
    MALLOC_INTERN(solvout->tasktab, solvout->tasknbr, Task);
    memcpy(solvout->tasktab, solvin->tasktab, solvout->tasknbr*sizeof(Task));

    /* Copy cblktab and bloktab */
    MALLOC_INTERN(solvout->cblktab, solvout->cblknbr+1, SolverCblk);
    memcpy(solvout->cblktab, solvin->cblktab,
           (solvout->cblknbr+1)*sizeof(SolverCblk));

    MALLOC_INTERN(solvout->bloktab, solvout->bloknbr, SolverBlok);
    memcpy(solvout->bloktab, solvin->bloktab,
           solvout->bloknbr*sizeof(SolverBlok));

    MALLOC_INTERN(solvout->browtab, solvout->brownbr, pastix_int_t);
    memcpy(solvout->browtab, solvin->browtab,
           solvout->brownbr*sizeof(pastix_int_t));

    solvblok = solvout->bloktab;
    for (solvcblk = solvout->cblktab; solvcblk  < solvout->cblktab + solvout->cblknbr; solvcblk++) {
        pastix_int_t bloknbr = (solvcblk+1)->fblokptr - solvcblk->fblokptr;
        solvcblk->fblokptr = solvblok;
        solvblok+= bloknbr;

        if ( flttype == -1 ) {
            solvcblk->lcoeftab = NULL;
            solvcblk->ucoeftab = NULL;
        }
        else {
            void *lcoeftab = solvcblk->lcoeftab;
            void *ucoeftab = solvcblk->ucoeftab;
            size_t size = cblk_colnbr( solvcblk ) * solvcblk->stride
                * pastix_size_of( flttype );

            if ( lcoeftab ) {
                MALLOC_INTERN( solvcblk->lcoeftab, size, char );
                memcpy(solvcblk->lcoeftab, lcoeftab, size );
            }
            else {
                solvcblk->lcoeftab = NULL;
            }

            if ( ucoeftab ) {
                MALLOC_INTERN( solvcblk->ucoeftab, size, char );
                memcpy(solvcblk->ucoeftab, ucoeftab, size );
            }
            else {
                solvcblk->ucoeftab = NULL;
            }
        }
    }
    solvcblk->fblokptr = solvblok;

    /* Copy ftgttab */
    if (solvout->ftgtnbr != 0)
    {
        MALLOC_INTERN(solvout->ftgttab, solvout->ftgtnbr, solver_ftgt_t);
        memcpy(solvout->ftgttab, solvin->ftgttab,
               solvout->ftgtnbr*sizeof(solver_ftgt_t));
    }
    /* Copy infotab of fan intarget */
    /*
     for(i=0;i<solvin->ftgtnbr;i++)
       memcpy(solvout->ftgttab[i].infotab, solvin->ftgttab[i].infotab, FTGT_MAXINFO*sizeof(pastix_int_t));
     */

    /* Copy indtab */
    MALLOC_INTERN(solvout->indtab, solvout->indnbr, pastix_int_t);
    memcpy(solvout->indtab, solvin->indtab, solvout->indnbr*sizeof(pastix_int_t));


    /* Copy ttsktab & ttsknbr */
    if (solvout->bublnbr>0)
    {
        MALLOC_INTERN(solvout->ttsknbr, solvout->bublnbr, pastix_int_t);
        memcpy(solvout->ttsknbr, solvin->ttsknbr, solvout->bublnbr*sizeof(pastix_int_t));
        MALLOC_INTERN(solvout->ttsktab, solvout->bublnbr, pastix_int_t*);

        for (i=0;i<solvout->bublnbr;i++)
        {
            solvout->ttsktab[i] = NULL;
            MALLOC_INTERN(solvout->ttsktab[i], solvout->ttsknbr[i], pastix_int_t);
            memcpy(solvout->ttsktab[i], solvin->ttsktab[i],
                   solvout->ttsknbr[i]*sizeof(pastix_int_t));
        }
    }
    else
    {
        solvout->ttsknbr = NULL;
        solvout->ttsktab = NULL;
    }

    MALLOC_INTERN(solvout->proc2clust, solvout->procnbr, pastix_int_t);
    memcpy(solvout->proc2clust, solvin->proc2clust,
           solvout->procnbr * sizeof(pastix_int_t));
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Generate a copy of a solver matrix structure.
 *
 * Every data is copied, event the coefficient if they are allocated and
 * initialized.
 * @warning This function is not able to copy a solver matrix with low rank
 * blocks yet.
 *
 *******************************************************************************
 *
 * @param[in] solvin
 *          The solver matrix structure to duplicate.
 *
 * @param[in] flttype
 *          The floating point arithmetic sued in the input solver matrix to
 *          know the size of the memory space to duplicate for the coefficients.
 *
 *******************************************************************************
 *
 * @return The pointer to the solver matrix internally allocated and that is a
 *         copy of the input solver. This pointer is NULL if the copy failed.
 *
 *******************************************************************************/
SolverMatrix *
solverCopy( const SolverMatrix *solvin,
            int                 flttype )
{
    SolverMatrix *solvout;

    MALLOC_INTERN(solvout, 1, SolverMatrix);
    memcpy(solvout, solvin, sizeof(SolverMatrix));

    solver_copy( solvin, solvout, flttype );

    return solvout;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Realloc in a contiguous way a given solver structure.
 *
 * All internal data of the solver structure are reallocated in a contiguous
 * manner to avoid the possible fragmentation from the initialization at
 * runtime.
 * @warning This function is not able to copy a solver matrix with low rank
 * blocks yet.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          On entry, the solver matrix to reallocate.
 *          On exit, the solver matrix with all internal data reallocated.
 *
 *******************************************************************************/
void
solverRealloc( SolverMatrix *solvmtx )
{
    SolverMatrix *tmp;

    MALLOC_INTERN(tmp, 1, SolverMatrix);
    /** copy general info **/
    memcpy(tmp, solvmtx, sizeof(SolverMatrix));

    solver_copy( tmp, solvmtx, -1 );

    /** Free the former solver matrix **/
    solverExit(tmp);
    memFree_null(tmp);
}
