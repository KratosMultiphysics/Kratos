/**
 *
 * @file solver_check.c
 *
 * PaStiX check function fo rthe solver structure.
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
#include <stdio.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "simu.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Checks the consistency of the given solver matrix structure.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure to check.
 *
 *******************************************************************************
 *
 * @retval 0 if the structure is correct
 * @retval 1 if incorrect
 *
 *******************************************************************************/
int
solverCheck(const SolverMatrix *solvmtx)
{
    pastix_int_t i, j, k = 0;
    pastix_int_t cblknum, bloknum, ftgtnum;
    pastix_int_t indnum, tasknum;
    pastix_int_t total;

    /** Check the task **/
    for(i=0;i<solvmtx->tasknbr;i++)
    {
        cblknum = solvmtx->tasktab[i].cblknum;
        bloknum = solvmtx->tasktab[i].bloknum;
        indnum  = solvmtx->tasktab[i].indnum;
        assert(cblknum < solvmtx->cblknbr);
        assert(bloknum < solvmtx->bloknbr);
        if(indnum >= solvmtx->indnbr) {
            printf("tasknbr %ld Tasknum %ld type %ld indnum %ld indnbr %ld \n", (long)solvmtx->tasknbr, (long)i, (long)solvmtx->tasktab[i].taskid, (long)indnum, (long)solvmtx->indnbr);
        }
        /** OIMBE ce test foire (ligne precedente) si il y a du 1D jusqu'au bout !! mais a priori on s'en fout **/
        /*if(solvmtx->tasktab[i].taskid != 0 && solvmtx->cblktab[solvmtx->tasktab[i].cblknum]
         assert(indnum < solvmtx->indnbr);*/
        if(indnum >= solvmtx->indnbr) {
            printf("cblknbr %ld cblknum %ld indnum %ld indnbr %ld \n", (long)solvmtx->cblknbr, (long)solvmtx->tasktab[i].cblknum, (long)indnum, (long)solvmtx->indnbr);
        }
        assert(solvmtx->tasktab[i].taskid >= 0);
        assert(solvmtx->tasktab[i].prionum >= 0);
        switch(solvmtx->tasktab[i].taskid)
        {
        case COMP_1D:
            assert(bloknum  == (solvmtx->cblktab[cblknum].fblokptr-solvmtx->bloktab));
            for(j = bloknum+1;j<solvmtx->cblktab[cblknum+1].fblokptr-solvmtx->bloktab;j++)
            {
                for(k=j;k<solvmtx->cblktab[cblknum+1].fblokptr-solvmtx->bloktab;k++)
                {
                    /*
                     * Test for ILU(k) to restore when we will do it
                     */
                    /* if(solvmtx->indtab[indnum] > solvmtx->ftgtnbr) /\** No ftgt **\/ */
                    /* { */
                    /*     indnum++; */
                    /*     continue; */
                    /* } */
                    if(solvmtx->indtab[indnum] < 0)
                    {
                        tasknum = -solvmtx->indtab[indnum];
                        switch(solvmtx->tasktab[tasknum].taskid)
                        {
                        case COMP_1D:
                        {
                            pastix_int_t facebloknum, facecblknum;
                            facecblknum = solvmtx->bloktab[j].fcblknm;
                            assert(facecblknum >= 0);
                            assert(facecblknum == solvmtx->tasktab[tasknum].cblknum);
                            facebloknum = solvmtx->cblktab[facecblknum].fblokptr-solvmtx->bloktab;

                            while ( !(( (solvmtx->bloktab[k].frownum >= solvmtx->bloktab[facebloknum].frownum) &&
                                        (solvmtx->bloktab[k].frownum <= solvmtx->bloktab[facebloknum].lrownum) ) ||
                                      ( (solvmtx->bloktab[k].lrownum >= solvmtx->bloktab[facebloknum].frownum) &&
                                        (solvmtx->bloktab[k].lrownum <= solvmtx->bloktab[facebloknum].lrownum) ) ||
                                      ( (solvmtx->bloktab[k].frownum <= solvmtx->bloktab[facebloknum].frownum) &&
                                        (solvmtx->bloktab[k].lrownum >= solvmtx->bloktab[facebloknum].lrownum) )) )
                            {
                                facebloknum++;
                            }

                            assert(solvmtx->bloktab[k].frownum >= solvmtx->bloktab[facebloknum].frownum);
                            assert(solvmtx->bloktab[k].lrownum <= solvmtx->bloktab[facebloknum].lrownum);
                            assert(solvmtx->bloktab[j].frownum >= solvmtx->cblktab[facecblknum].fcolnum);
                            assert(solvmtx->bloktab[j].lrownum <= solvmtx->cblktab[facecblknum].lcolnum);
                        }
                        break;
                        default:
                            errorPrint("tasknum %ld tasknbr %ld taskid %ld",
                                       (long)tasknum, (long)solvmtx->tasknbr,
                                       (long)solvmtx->tasktab[tasknum].taskid);
                            EXIT(MOD_BLEND,INTERNAL_ERR);
                        }
                    }
                    else
                    {
                        ftgtnum = solvmtx->indtab[indnum];


                        assert(solvmtx->bloktab[j].frownum >= solvmtx->ftgttab[ftgtnum].infotab[FTGT_FCOLNUM]);
                        assert(solvmtx->bloktab[j].lrownum <= solvmtx->ftgttab[ftgtnum].infotab[FTGT_LCOLNUM]);
                        assert(solvmtx->bloktab[k].frownum >= solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM]);
                        assert(solvmtx->bloktab[k].lrownum <= solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM]);

                        solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT]--;

#if defined(PASTIX_DEBUG_BLEND)
                        if( solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT] == 0) {
                            if(solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM] != solvmtx->tasktab[i].prionum) {
                                fprintf( stdout, "Task1D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n",
                                         (long)i, (long)ftgtnum, (long)solvmtx->tasktab[i].prionum,
                                         (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM] );
                            }
                        }
#endif
                    }
                    indnum++;
                }
            }
            break;
        default:
            fprintf(stderr, "solver_check: The task %ld has no type \n", (long)i);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    }
    for(i=0;i<solvmtx->ftgtnbr;i++)
    {
        assert(solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]>0);
        assert(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0);
    }

    /** Reset the ftgt ctrbcnt **/
    for(i=0;i<solvmtx->ftgtnbr;i++)
    {
        assert(solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]>0);
        assert(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0);
        solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];
    }

    /** Test the task partition on the thread of the cluster **/
    total = 0;
    for(i=0;i<solvmtx->bublnbr;i++){
        printf("i = %d, ttsknbr = %d\n", (int)i, (int)(solvmtx->ttsknbr[i]));
        total += solvmtx->ttsknbr[i];
    }
    if(total != solvmtx->tasknbr) {
        fprintf(stderr, " total %ld tasknbr %ld \n", (long)total, (long)solvmtx->tasknbr);
    }

    assert(total == solvmtx->tasknbr);

#if defined(PASTIX_DEBUG_BLEND)
    {
        pastix_int_t * flag;
        MALLOC_INTERN(flag, solvmtx->tasknbr, pastix_int_t);
        bzero(flag, sizeof(pastix_int_t)*solvmtx->tasknbr);

        for(i=0;i<solvmtx->bublnbr;i++) {
            for(j=0;j<solvmtx->ttsknbr[i];j++)
            {
                if(flag[solvmtx->ttsktab[i][j]] != 0) {
                    fprintf(stderr, "flag %ld thread %ld task %ld already on another thread \n", (long)flag[solvmtx->ttsktab[i][j]], (long)i, (long)solvmtx->ttsktab[i][j]);
                }
                flag[solvmtx->ttsktab[i][j]]++;
            }
        }

        for(i=0;i<solvmtx->tasknbr;i++) {
            assert(flag[i] == 1);
        }
        memFree(flag);
    }
#endif

#if defined(PASTIX_DYNSCHED)
    {
        int k, bubnbr = solvmtx->bublnbr;

        for(k = 0; k<bubnbr; k++)
        {
            int father = BFATHER( solvmtx->btree, k );
            if ( (father != -1) &&
                 (  solvmtx->btree->nodetab[k].priomax >  solvmtx->btree->nodetab[father].priomin ) )
            {
                fprintf(stderr, "We have a problem of task distribution\n"
                        " Bubble[%d] priorities (%ld,%ld) intersect with bubble[%d](%ld,%ld)\n",
                        k,      (long)solvmtx->btree->nodetab[k].priomin,      (long)solvmtx->btree->nodetab[k].priomax,
                        father, (long)solvmtx->btree->nodetab[father].priomin, (long)solvmtx->btree->nodetab[father].priomax) ;
            }
        }
    }
#endif

    return 0;
}
