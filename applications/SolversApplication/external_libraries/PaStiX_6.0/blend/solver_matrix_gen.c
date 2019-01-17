/**
 *
 * @file solver_matrix_gen.c
 *
 * PaStiX solver structure generation function.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>

#include "common.h"
#include "pastix/order.h"
#include "cost.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"
#include "extendVector.h"
#include "elimintree.h"
#include "cand.h"
#include "blendctrl.h"
#include "simu.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure
 *
 * This function takes all the global preprocessing steps: the symbol matrix,
 * and the resul of the simulation step to generate the solver matrix that hold
 * only local information to each PaStiX process.
 *
 *******************************************************************************
 *
 * @param[in] clustnum
 *          The index of the local PaStiX process.
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
int
solverMatrixGen( pastix_int_t           clustnum,
                 SolverMatrix          *solvmtx,
                 const symbol_matrix_t *symbmtx,
                 const pastix_order_t  *ordeptr,
                 const SimuCtrl        *simuctrl,
                 const BlendCtrl       *ctrl )
{
    pastix_int_t  p, c;
    pastix_int_t  cursor;
    pastix_int_t  i, j, k;
    pastix_int_t  ftgtnum      = 0;
    pastix_int_t  coefnbr      = 0;
    pastix_int_t  nodenbr      = 0;
    pastix_int_t  odb_nbr      = 0;
    pastix_int_t  cblknum      = 0;
    pastix_int_t  tasknum      = 0;
    pastix_int_t  brownum      = 0;
    pastix_int_t  indnbr       = 0;
    pastix_int_t *cblklocalnum = NULL;
    pastix_int_t *bloklocalnum = NULL;
    pastix_int_t *tasklocalnum = NULL;
    pastix_int_t *ftgtlocalnum = NULL;
    pastix_int_t  flaglocal    = 0;
    pastix_int_t  dof = symbmtx->dof;
    (void)ordeptr;

    assert( dof == 1 );
    assert( symbmtx->baseval == 0 );

    solverInit(solvmtx);

#ifdef PASTIX_DYNSCHED
    solvmtx->btree    = ctrl->btree;
#endif
    solvmtx->clustnum = ctrl->clustnum;
    solvmtx->clustnbr = ctrl->clustnbr;
    solvmtx->procnbr  = ctrl->total_nbcores;
    solvmtx->thrdnbr  = ctrl->local_nbthrds;
    solvmtx->bublnbr  = ctrl->local_nbctxts;
    solvmtx->ftgtcnt  = simuctrl->ftgtcnt;

    /* Copy the vector used to get a cluster number from a processor number */
    MALLOC_INTERN(solvmtx->proc2clust, solvmtx->procnbr, pastix_int_t);
    memcpy(solvmtx->proc2clust, ctrl->core2clust, sizeof(pastix_int_t)*solvmtx->procnbr);

    /*
     * Compute local indexes to compress the symbol information into solver
     */
    {
        pastix_int_t *localindex;

        MALLOC_INTERN(localindex, ctrl->clustnbr, pastix_int_t);
        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );

        /* Compute local number of tasks on each cluster */
        MALLOC_INTERN(tasklocalnum, simuctrl->tasknbr, pastix_int_t);
        for(i=0; i<simuctrl->tasknbr; i++) {
            c = simuctrl->bloktab[simuctrl->tasktab[i].bloknum].ownerclust;

            tasklocalnum[i] = localindex[c];
            localindex[c]++;
        }
        solvmtx->tasknbr = localindex[clustnum];

        /* Compute the local numbering of the bloks and cblks on each cluster */
        MALLOC_INTERN(bloklocalnum, symbmtx->bloknbr, pastix_int_t);
        MALLOC_INTERN(cblklocalnum, symbmtx->cblknbr, pastix_int_t);

        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );
        cblknum = 0;
        for(i=0; i<symbmtx->cblknbr; i++)
        {
            flaglocal = 0;
            for(j = symbmtx->cblktab[i].bloknum; j<symbmtx->cblktab[i+1].bloknum; j++)
            {
                c = simuctrl->bloktab[j].ownerclust;
                bloklocalnum[j] = localindex[c];
                localindex[c]++;

                if (c == clustnum) {
                    flaglocal = 1;
                }
            }

            if(flaglocal) {
                cblklocalnum[i] = cblknum;
                cblknum++;
                brownum += symbmtx->cblktab[i+1].brownum - symbmtx->cblktab[i].brownum;
                assert( brownum <= symbmtx->cblktab[ symbmtx->cblknbr ].brownum );
            }
            else {
                cblklocalnum[i] = -i-1;
            }
        }
        solvmtx->cblknbr = cblknum;
        solvmtx->bloknbr = localindex[clustnum];
        solvmtx->brownbr = brownum;

        memFree_null(localindex);
    }

    /***************************************************************************
     * Fill in the local bloktab and cblktab
     */

    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk  );
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr,   SolverBlok  );
    MALLOC_INTERN(solvmtx->browtab, solvmtx->brownbr,   pastix_int_t);

    {
        SolverCblk    *solvcblk = solvmtx->cblktab;
        SolverBlok    *solvblok = solvmtx->bloktab;
        symbol_cblk_t *symbcblk = symbmtx->cblktab;
        symbol_blok_t *symbblok = symbmtx->bloktab;
        SimuBlok      *simublok = simuctrl->bloktab;
        Cand          *candcblk = ctrl->candtab;
#if defined(PASTIX_SUPERNODE_STATS)
        pastix_int_t  *sndetab = ordeptr->sndetab;
#endif
        pastix_int_t   blokamax = 0; /* Maximum area of a block in the global matrix */
        pastix_int_t   nbcblk2d = 0;
        pastix_int_t   nbblok2d = 0;

        solvmtx->cblkmax1d  = -1;
        solvmtx->cblkmin2d  = solvmtx->cblknbr;
        solvmtx->cblkmaxblk = 1;
        solvmtx->cblkschur  = solvmtx->cblknbr;
        cblknum = 0;
        brownum = 0;
        nodenbr = 0;
        coefnbr = 0;
        for(i=0; i<symbmtx->cblknbr; i++, symbcblk++, candcblk++)
        {
            SolverBlok  *fblokptr = solvblok;
            pastix_int_t fbloknum = symbcblk[0].bloknum;
            pastix_int_t lbloknum = symbcblk[1].bloknum;
            pastix_int_t stride   = 0;
            pastix_int_t fcolnum, lcolnum, nbcols;
            pastix_int_t layout2D = candcblk->cblktype & CBLK_LAYOUT_2D;
            pastix_int_t tasks2D  = candcblk->cblktype & CBLK_TASKS_2D;

            if ( symbmtx->dof < 0 ) {
                fcolnum = symbmtx->dofs[ symbcblk->fcolnum    ];
                lcolnum = symbmtx->dofs[ symbcblk->lcolnum + 1] - 1;
            }
            else {
                fcolnum = symbmtx->dof *  symbcblk->fcolnum;
                lcolnum = symbmtx->dof * (symbcblk->lcolnum+1) - 1;
            }
            nbcols = lcolnum - fcolnum + 1;
            flaglocal = 0;

            for( j=fbloknum; j<lbloknum; j++, symbblok++, simublok++ ) {
                pastix_int_t frownum, lrownum, nbrows;

                if ( symbmtx->dof < 0 ) {
                    frownum = symbmtx->dofs[ symbblok->frownum     ];
                    lrownum = symbmtx->dofs[ symbblok->lrownum + 1 ] - 1;
                }
                else {
                    frownum = symbmtx->dof *  symbblok->frownum;
                    lrownum = symbmtx->dof * (symbblok->lrownum+1) - 1;
                }
                nbrows = lrownum - frownum + 1;
                blokamax = pastix_imax( blokamax, nbrows * nbcols );

                if(simublok->ownerclust == clustnum)
                {
                    flaglocal = 1;

                    /* Init the blok */
                    solvblok->handler[0] = NULL;
                    solvblok->handler[1] = NULL;
                    solvblok->fcblknm = cblklocalnum[symbblok->fcblknm];
                    solvblok->lcblknm = cblklocalnum[symbblok->lcblknm];
                    solvblok->frownum = frownum;
                    solvblok->lrownum = lrownum;
                    solvblok->coefind = layout2D ? stride * nbcols : stride;
                    solvblok->browind = -1;
                    solvblok->gpuid   = GPUID_UNDEFINED;
                    solvblok->LRblock = NULL;

                    stride += nbrows;
                    solvblok++;
                }
            }
            if(flaglocal)
            {
                pastix_int_t brownbr;

                /*
                 * 2D tasks: Compute the number of cblk split in 2D tasks, and
                 * the smallest id
                 */
                if (tasks2D) {
                    if (cblknum < solvmtx->cblkmin2d) {
                        solvmtx->cblkmin2d = cblknum;
                    }
                    nbcblk2d++;
                    nbblok2d += solvblok - fblokptr;
                }
                else {
                    if (cblknum > solvmtx->cblkmax1d) {
                        solvmtx->cblkmax1d = cblknum;
                    }
                }

                /*
                 * Store first local cblk in Schur
                 */
                if (candcblk->cblktype & CBLK_IN_SCHUR) {
                    if (cblknum < solvmtx->cblkschur) {
                        solvmtx->cblkschur = cblknum;
                    }
                }

                /*
                 * Compute the maximum number of block per cblk for data
                 * structure in PaRSEC/StarPU
                 */
                if ((cblknum >= solvmtx->cblkmin2d) &&
                    ((solvblok - fblokptr) > solvmtx->cblkmaxblk) )
                {
                    solvmtx->cblkmaxblk = solvblok - fblokptr;
                }

                /* Init the cblk */
                solvcblk->lock     = PASTIX_ATOMIC_UNLOCKED;
                solvcblk->ctrbcnt  = -1;
                solvcblk->cblktype = candcblk->cblktype;
                solvcblk->gpuid    = GPUID_UNDEFINED;
                solvcblk->fcolnum  = fcolnum;
                solvcblk->lcolnum  = lcolnum;
                solvcblk->fblokptr = fblokptr;
                solvcblk->stride   = stride;
                solvcblk->lcolidx  = nodenbr;
                solvcblk->brownum  = brownum;
                solvcblk->lcoeftab = NULL;
                solvcblk->ucoeftab = NULL;
                solvcblk->gcblknum = i;

                assert( nodenbr == fcolnum );

#if defined(PASTIX_SUPERNODE_STATS)
                {
                    while ( sndetab[1] <= solvcblk->lcolnum ) {
                        sndetab++;
                        assert( (sndetab - ordeptr->sndetab) < ordeptr->sndenbr );
                    }
                    assert( (sndetab[0] <= solvcblk->fcolnum) &&
                            (sndetab[1] >  solvcblk->lcolnum) );
                    solvcblk->sndeidx = sndetab - ordeptr->sndetab;
                }
#else
                solvcblk->sndeidx = i;
#endif
                solvcblk->handler[0] = NULL;
                solvcblk->handler[1] = NULL;

                /*
                 * Copy browtab information
                 * In case of 2D tasks, we reorder the browtab to first store
                 * the 1D contributions, and then the 2D updates.
                 * This might also be used for low rank compression, to first
                 * accumulate small dense contributions, and then, switch to a
                 * low rank - low rank update scheme.
                 */
                brownbr = symbmtx->cblktab[i+1].brownum
                    -     symbmtx->cblktab[i].brownum;
                solvcblk->brown2d = solvcblk->brownum + brownbr;
                if (brownbr)
                {
                    if ( tasks2D ) {
                        SolverBlok *blok;
                        SolverCblk *cblk;
                        pastix_int_t j2d, j1d, j, *b;

                        j2d = -1;
                        /* First pass to copy 1D updates */
                        for( j=0, j1d=0; j<brownbr; j++ ) {
                            b = symbmtx->browtab + symbmtx->cblktab[i].brownum + j;
                            blok = solvmtx->bloktab + (*b);
                            cblk = solvmtx->cblktab + blok->lcblknm;
                            if (! (cblk->cblktype & CBLK_TASKS_2D) ) {
                                solvmtx->browtab[brownum + j1d] = (*b);
                                *b = -1;
                                j1d++;
                            }
                            else {
                                /* Store the first non 1D index to not rediscover the begining */
                                if (j2d == -1) {
                                    j2d = j;
                                }
                            }
                        }

                        /* Store the index of the first 2D contribution in the array */
                        assert(j1d <= brownbr);
                        solvcblk->brown2d = solvcblk->brownum + j1d;

                        /* Second pass to copy 2D updates to the end */
                        if (j2d != -1) {
                            for( j=j2d; j<brownbr; j++ ) {
                                b = symbmtx->browtab + symbmtx->cblktab[i].brownum + j;
                                if (*b != -1) {
                                    solvmtx->browtab[brownum + j1d] = *b;
                                    j1d++;
                                }
                            }
                        }
                        assert(j1d == brownbr);
                    }
                    else {
                        memcpy( solvmtx->browtab + brownum,
                                symbmtx->browtab + symbmtx->cblktab[i].brownum,
                                brownbr * sizeof(pastix_int_t) );
                    }
                }
                brownum += brownbr;
                assert( brownum <= solvmtx->brownbr );

                assert(solvcblk->brown2d >= solvcblk->brownum);
                assert(solvcblk->brown2d <= solvcblk->brownum + brownbr);

                /* Extra statistic informations */
                nodenbr += nbcols;
                coefnbr += stride * nbcols;

                cblknum++; solvcblk++;
            }
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if (cblknum > 0)
        {
            solvcblk->lock     = PASTIX_ATOMIC_UNLOCKED;
            solvcblk->ctrbcnt  = -1;
            solvcblk->cblktype = 0;
            solvcblk->gpuid    = GPUID_NONE;
            solvcblk->fcolnum  = solvcblk[-1].lcolnum + 1;
            solvcblk->lcolnum  = solvcblk[-1].lcolnum + 1;
            solvcblk->fblokptr = solvblok;
            solvcblk->stride   = 0;
            solvcblk->lcolidx  = nodenbr;
            solvcblk->brownum  = symbcblk->brownum;
            solvcblk->brown2d  = symbcblk->brownum;
            solvcblk->gcblknum = -1;
            solvcblk->lcoeftab = NULL;
            solvcblk->ucoeftab = NULL;
            solvcblk->handler[0] = NULL;
            solvcblk->handler[1] = NULL;
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;
        solvmtx->arftmax = blokamax;

        solvmtx->nb2dcblk = nbcblk2d;
        solvmtx->nb2dblok = nbblok2d;

        assert( solvmtx->cblkmax1d+1 >= solvmtx->cblkmin2d );
        assert( solvmtx->cblknbr == cblknum );
        assert( solvmtx->bloknbr == solvblok - solvmtx->bloktab );
    }

    /*
     * Update browind fields
     */
    for(i=0; i<solvmtx->brownbr; i++)
    {
        pastix_int_t bloknum = solvmtx->browtab[i];
        if ( simuctrl->bloktab[bloknum].ownerclust == clustnum ) {
            solvmtx->bloktab[ bloknum ].browind = i;
        }
    }

    /*
     * Fill in tasktab
     */
    MALLOC_INTERN(solvmtx->tasktab, solvmtx->tasknbr+1, Task);
    {
        SimuTask    *simutask = simuctrl->tasktab;
        Task        *solvtask = solvmtx->tasktab;
        pastix_int_t nbftmax  = 0;

        tasknum = 0;
        ftgtnum = 0;
        indnbr  = 0;

        for(i=0; i<simuctrl->tasknbr; i++, simutask++)
        {
            nbftmax = pastix_imax( nbftmax, simutask->ftgtcnt );
            if( simuctrl->bloktab[ simutask->bloknum ].ownerclust == clustnum )
            {
                assert( tasknum == tasklocalnum[i] );

                solvtask->taskid  = COMP_1D;
                solvtask->prionum = simutask->prionum;
                solvtask->cblknum = cblklocalnum[ simutask->cblknum ];
                solvtask->bloknum = bloklocalnum[ simutask->bloknum ];
                solvtask->ftgtcnt = simutask->ftgtcnt;
                solvtask->ctrbcnt = simutask->ctrbcnt;
                solvtask->indnum  = indnbr;

                /*
                 * Count number of index needed in indtab:
                 *  => number of off-diagonal block below the block (including the block itself)
                 */
                odb_nbr = symbmtx->cblktab[ simutask->cblknum + 1 ].bloknum - simutask->bloknum - 1;

                indnbr += (odb_nbr*(odb_nbr+1))/2;
                tasknum++; solvtask++;
            }
        }
        assert(tasknum == solvmtx->tasknbr);

        /* One more task to avoid side effect */
        solvtask->taskid  = -1;
        solvtask->prionum = -1;
        solvtask->cblknum = solvmtx->cblknbr+1;
        solvtask->bloknum = solvmtx->bloknbr+1;
        solvtask->ftgtcnt = 0;
        solvtask->ctrbcnt = 0;
        solvtask->indnum  = indnbr;

        /* Store the final indnbr */
        solvmtx->indnbr  = indnbr;
        solvmtx->nbftmax = nbftmax;
    }

    /*
     * Fill in the ttsktab arrays (one per thread)
     *
     * TODO: This would definitely be better if each thread was initializing
     * it's own list on its own memory node.
     */
    {
        SimuProc *simuproc = &(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum]);

        /* Number of processor in this cluster */
        k = solvmtx->bublnbr;
        MALLOC_INTERN(solvmtx->ttsknbr, k, pastix_int_t  );
        MALLOC_INTERN(solvmtx->ttsktab, k, pastix_int_t* );

        for(p=0; p<k; p++, simuproc++)
        {
            pastix_int_t priomin = PASTIX_INT_MAX;
            pastix_int_t priomax = 0;
            pastix_int_t ttsknbr = extendint_Size( simuproc->tasktab );
            pastix_int_t j, jloc;

            solvmtx->ttsknbr[p] = ttsknbr;

            if(ttsknbr > 0) {
                MALLOC_INTERN(solvmtx->ttsktab[p], ttsknbr, pastix_int_t);
            }
            else {
                solvmtx->ttsktab[p] = NULL;
            }

            for(i=0; i<ttsknbr; i++)
            {
                j    = extendint_Read(simuproc->tasktab, i);
                jloc = tasklocalnum[j];
                solvmtx->ttsktab[p][i] = jloc;

#if defined(PASTIX_DYNSCHED)
                solvmtx->tasktab[jloc].threadid = p;
#endif
                priomax = pastix_imax( solvmtx->tasktab[jloc].prionum, priomax );
                priomin = pastix_imin( solvmtx->tasktab[jloc].prionum, priomin );
            }

#if defined(PASTIX_DYNSCHED)
            solvmtx->btree->nodetab[p].priomin = priomin;
            solvmtx->btree->nodetab[p].priomax = priomax;
#endif
        }
    }

    /*
     * Fill in ftgttab
     */
    {
        solvmtx->ftgtnbr = 0;
        solvmtx->ftgttab = NULL;

        /* Compute local number of outgoing contributions */
        for(c=0; c<ctrl->clustnbr; c++) {
            if(c == clustnum) {
                assert( extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c])) == 0 );
                continue;
            }
            solvmtx->ftgtnbr += extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));
        }

        if(solvmtx->ftgtnbr > 0) {
            SimuCluster   *simuclust = &(simuctrl->clustab[clustnum]);
            solver_ftgt_t *solvftgt;
            pastix_int_t   ftgtnbr;

            MALLOC_INTERN(solvmtx->ftgttab, solvmtx->ftgtnbr, solver_ftgt_t);

            /* Allocate array to store local indexes */
            ftgtnbr = simuctrl->bloktab[symbmtx->bloknbr].ftgtnum;
            MALLOC_INTERN(ftgtlocalnum, ftgtnbr, pastix_int_t);
            memset(ftgtlocalnum, -1, ftgtnbr * sizeof(pastix_int_t));

            cursor = 0;
            solvftgt = solvmtx->ftgttab;

            for(c=0; c<ctrl->clustnbr; c++)
            {
                ftgtnbr = extendint_Size(&(simuclust->ftgtsend[c]));
                for(i=0; i<ftgtnbr; i++)
                {
                    ftgtnum = extendint_Read(&(simuclust->ftgtsend[c]), i);
                    ftgtlocalnum[ftgtnum] = cursor;

                    /* Copy information computed during simulation */
                    memcpy(solvftgt->infotab, simuctrl->ftgttab[ftgtnum].ftgt.infotab, FTGT_MAXINFO*sizeof(pastix_int_t));

                    /* Update with Degre of freedoms */
                    solvftgt->infotab[FTGT_FCOLNUM] *= dof;
                    solvftgt->infotab[FTGT_LCOLNUM] *= dof;
                    solvftgt->infotab[FTGT_LCOLNUM] += dof - 1;
                    solvftgt->infotab[FTGT_FROWNUM] *= dof;
                    solvftgt->infotab[FTGT_LROWNUM] *= dof;
                    solvftgt->infotab[FTGT_LROWNUM] += dof - 1;

                    /* Convert to local numbering */
                    solvftgt->infotab[FTGT_TASKDST] = tasklocalnum[solvftgt->infotab[FTGT_TASKDST]];
                    solvftgt->infotab[FTGT_BLOKDST] = bloklocalnum[solvftgt->infotab[FTGT_BLOKDST]];

                    /* Restore ctrbcnt (modified durind simulation) */
                    solvftgt->infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR];
                    solvftgt->coeftab = NULL;

                    cursor++; solvftgt++;
                }
            }
        }
    }

    /*
     * Fill in indtab
     */
    {
        solvmtx->indtab = NULL;
        if (solvmtx->indnbr) {
            MALLOC_INTERN(solvmtx->indtab, solvmtx->indnbr, pastix_int_t);
        }

        indnbr = 0;
        for(i=0; i<simuctrl->tasknbr; i++)
        {
            pastix_int_t bloknum = simuctrl->tasktab[i].bloknum;
            pastix_int_t cblknum = simuctrl->tasktab[i].cblknum;

            if(simuctrl->bloktab[bloknum].ownerclust != clustnum)
                continue;

            assert(indnbr == solvmtx->tasktab[tasklocalnum[i]].indnum);
            assert(bloklocalnum[simuctrl->tasktab[i].bloknum] == solvmtx->tasktab[tasklocalnum[i]].bloknum);
            assert(cblklocalnum[simuctrl->tasktab[i].cblknum] == solvmtx->tasktab[tasklocalnum[i]].cblknum);

            {
                pastix_int_t fbloknum = symbmtx->cblktab[cblknum  ].bloknum+1;
                pastix_int_t lbloknum = symbmtx->cblktab[cblknum+1].bloknum;

                /*
                 * For each couple (bloknum,j) \ j>=bloknum of off-diagonal
                 * block, check where goes the contribution
                 */
                for(bloknum=fbloknum; bloknum<lbloknum; bloknum++)
                {
                    pastix_int_t firstbloknum = 0;
                    pastix_int_t facebloknum  = 0;

                    for(j=bloknum; j<lbloknum; j++)
                    {
                        facebloknum = pastixSymbolGetFacingBloknum(symbmtx, bloknum, j, firstbloknum, ctrl->ricar);

                        if(facebloknum >= 0) {
                            firstbloknum = facebloknum;

                            if(simuctrl->bloktab[facebloknum].ownerclust != clustnum)
                            {
                                solvmtx->indtab[indnbr] = ftgtlocalnum[CLUST2INDEX(facebloknum, clustnum)];
                                assert(solvmtx->indtab[indnbr] < solvmtx->ftgtnbr);
                            }
                            else
                            {
                                solvmtx->indtab[indnbr] = - tasklocalnum[simuctrl->bloktab[facebloknum].tasknum];
                                assert( (- solvmtx->indtab[indnbr]) < solvmtx->tasknbr );
                            }
                        }
                        else {
                            /* The facing block does not exist */
                            solvmtx->indtab[indnbr] =  solvmtx->ftgtnbr+1;
                        }
                        indnbr++;
                    }
                }
            }
        }
        assert(indnbr == solvmtx->indnbr);
    }

    /*
     * Compute the maximum area of the temporary buffer used during computation
     *
     * During this loop, we compute the maximum area that will be used as
     * temporary buffers, and statistics:
     *    - diagmax: Only for hetrf/sytrf factorization, this the maximum size
     *               of a panel of MAXSIZEOFBLOCKS width in a diagonal block
     *    - gemmmax: For all, this is the maximum area used to compute the
     *               compacted gemm on a CPU.
     *
     * Rk: This loop is not merged with the main block loop, since strides have
     * to be peviously computed.
     */
    {
        SolverCblk  *solvcblk = solvmtx->cblktab;
        SolverBlok  *solvblok = solvmtx->bloktab;
        pastix_int_t gemmmax = 0;
        pastix_int_t offdmax = 0;
        pastix_int_t blokmax = 0;
        pastix_int_t gemmarea, offdarea, cblk_m, acc_m;

        /* Let's keep the block dimensions to print statistics informations */
        pastix_int_t maxg_m = 0;
        pastix_int_t maxg_n = 0;
        pastix_int_t maxo_m = 0;
        pastix_int_t maxo_n = 0;

        for(i=0;i<solvmtx->cblknbr;i++, solvcblk++)
        {
            SolverBlok *lblok = solvcblk[1].fblokptr;
            pastix_int_t m = solvcblk->stride;
            pastix_int_t n = cblk_colnbr( solvcblk );
            pastix_int_t k = blok_rownbr( solvblok );

            m -= n;

            /*
             * Compute the surface of the off-diagonal block in a panel for
             * LDL^[th] factorizations
             */
            offdarea = m * n;
            if ( offdarea > offdmax ) {
                offdmax = offdarea;
                maxo_m = m;
                maxo_n = n;
            }

            /*
             * Compute the maximum area for 1d temporary workspace in GEMM
             */
            solvblok++;
            cblk_m = -1;
            acc_m  = 0;
            for( ; solvblok<lblok; solvblok++ ) {
                k = blok_rownbr( solvblok );

                /*
                 * Temporary workspace for GEMM
                 * m+1 to store the diagonal in case of GEMDM
                 */
                if ( !(solvcblk->cblktype & CBLK_LAYOUT_2D) ) {
                    gemmarea = (m+1) * k;
                    if ( gemmarea > gemmmax ) {
                        gemmmax = gemmarea;
                        maxg_m = m;
                        maxg_n = k;
                    }
                }

                /*
                 * Max size for off-diagonal blocks for 2-terms version of the
                 * 2D LDL
                 */
                if ( solvcblk->cblktype & (CBLK_TASKS_2D | CBLK_COMPRESSED) ) {
                    if ( solvblok->fcblknm == cblk_m ) {
                        acc_m += k;
                    }
                    else {
                        cblk_m = solvblok->fcblknm;
                        acc_m = k;
                    }
                    blokmax = pastix_imax( n * acc_m, blokmax );
                }
                m -= k;
            }
        }

        solvmtx->offdmax = offdmax;
        solvmtx->gemmmax = gemmmax;
        solvmtx->blokmax = blokmax;
#if defined(PASTIX_DEBUG_BLEND)
        if (ctrl->iparm[IPARM_VERBOSE]>PastixVerboseNo) {
            pastix_print(clustnum, 0,
                         "Coefmax: off-diagonal %6ld ((%4ld+1) x %3ld)\n"
                         "         update       %6ld ( %4ld    x %3ld)\n",
                         (long)offdmax, (long)maxo_m, (long)maxo_n,
                         (long)gemmmax, (long)maxg_m, (long)maxg_n );
        }
#endif
        (void)maxo_m;
        (void)maxo_n;
        (void)maxg_m;
        (void)maxg_n;
    }

    memFree_null(cblklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);

    if(ftgtlocalnum != NULL)
        memFree_null(ftgtlocalnum);

    return PASTIX_SUCCESS;
}
