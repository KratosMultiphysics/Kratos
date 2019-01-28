/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
#include "slu_mt_ddefs.h"

void
pxgstrf_scheduler(const int_t pnum, const int_t n, const int_t *etree, 
		  int_t *cur_pan, int_t *bcol, pxgstrf_shared_t *pxgstrf_shared)
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * =======
 *
 * pxgstrf_scheduler() gets a panel for the processor to work on. 
 * It schedules a panel in decreasing order of priority:
 *   (1) the current panel's parent, if it can be done without pipelining
 *   (2) any other panel in the queue that can be done without pipelining
 *       ("CANGO" status)
 *   (3) any other panel in the queue that can be done with pipelining
 *       ("CANPIPE" status)
 *
 * Arguments
 * =========
 * pnum    (input) int_t
 *         Processor number.
 *
 * n       (input) int_t
 *         Column dimension of the matrix.
 *
 * etree   (input) int_t*
 *         Elimination tree of A'*A, size n.
 *         Note: etree is a vector of parent pointers for a forest whose
 *         vertices are the integers 0 to n-1; etree[root] = n.
 *
 * cur_pan (input/output) int_t*
 *         On entry, the current panel just finished by this processor;
 *         On exit, [0, n-1]: the new panel to work on;
 *                  EMPTY:    failed to get any work, will try later;
 *                  n:        all panels are taken; ready to terminate.
 *
 * taskq   (input/output) queue_t*
 *         Global work queue.
 *
 * fb_cols (input/output) int_t*
 *         The farthest busy descendant of each (leading column of the) panel.
 *
 * bcol    (output) int_t*
 *         The most distant busy descendant of cur_pan in the *linear*
 *         pipeline of busy descendants. If all proper descendants of
 *         cur_pan are done, bcol is returned equal to cur_pan.
 *
 * Defining terms
 * ==============
 *   o parent(panel) = etree(last column in the panel)
 *   o the kids of a panel = collective kids of all columns in the panel
 *     kids[REP] = SUM_{j in panel} ( kids[j] ) 
 *   o linear pipeline - what does it mean in the panel context?
 *       if ukids[REP] = 0, then the panel becomes a leaf (CANGO)
 *       if ukids[REP] = 1 && ukids[firstcol] = 1, then the panel can
 *                       be taken with pipelining (CANPIPE)
 *
 * NOTES
 * =====
 *   o When a "busy" panel finishes, if its parent has only one remaining
 *     undone child there is no check to see if the parent should change
 *     from "unready" to "canpipe". Thus a few potential pipelinings will
 *     be lost, but checking out this pipeline opportunity may be costly.
 *
 */

    register int_t dad, dad_ukids, jcol, w, j;
    int_t *fb_cols = pxgstrf_shared->fb_cols;
    queue_t *taskq = &pxgstrf_shared->taskq;
    Gstat_t *Gstat = pxgstrf_shared->Gstat;
#ifdef PROFILE
    double t;
#endif

    jcol = *cur_pan;
    if ( jcol != EMPTY ) {
#ifdef DOMAINS
	if ( in_domain[jcol] == TREE_DOMAIN )
	    dad = etree[jcol];
	else
#endif
	    dad = DADPANEL (jcol);
    }

    /* w_top = sp_ienv(1)/2;
       if ( w_top == 0 ) w_top = 1;*/

#ifdef PROFILE
    TIC(t);
#endif
#if ( MACH==SUN )
    mutex_lock( &pxgstrf_shared->lu_locks[SCHED_LOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_lock( &pxgstrf_shared->lu_locks[SCHED_LOCK] );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock(pxgstrf_shared->lu_locks[SCHED_LOCK])
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard SCHED_LOCK
#elif ( MACH==OPENMP )
#pragma omp critical ( SCHED_LOCK )
#endif    

{   /* ---- START CRITICAL SECTION ---- */
    
    /* Update the status of the current panel and its parent, so that
     * the other processors waiting on it can proceed.
     * If all siblings are done, and dad is not busy, then take dad.
     */
    if ( jcol != EMPTY ) { /* jcol was just finished by this processor */    
	dad_ukids = --pxgstrf_shared->pan_status[dad].ukids;
	
#if ( DEBUGlevel>=1 )
	printf("(%d) DONE %d in Scheduler(), dad %d, STATE %d, dad_ukids %d\n",
	       pnum, jcol, dad, STATE(dad), dad_ukids);
#endif	

	if ( dad_ukids == 0 && STATE( dad ) > BUSY ) { /* dad not started */
	    jcol = dad;
#if ( DEBUGlevel>=1 )
	    printf("(%d) Scheduler[1] Got dad %d, STATE %d\n",
		   pnum, jcol, STATE(dad));
#endif
#ifdef PROFILE
	    ++(Gstat->panhows[DADPAN]);
#endif	    
	} else {
	    /* Try to get a panel from the task Q. */
	    while ( 1 ) {
		/*>>if ( (j = Dequeue(taskq, &item)) == EMPTY ) {*/
		if ( taskq->count <= 0 ) {
		    jcol = EMPTY;
		    break;
		} else {
		    jcol = taskq->queue[taskq->head++];
		    --taskq->count;
		    if ( STATE( jcol ) >= CANGO ) { /* CANGO or CANPIPE */
#if ( DEBUGlevel>=1 )
			printf("(%d) Dequeue[1] Got %d, STATE %d, Qcount %d\n",
			       pnum, jcol, STATE(jcol), j);
#endif
#ifdef PROFILE
			if (STATE( jcol ) == CANGO) ++(Gstat->panhows[NOPIPE]);
			else ++(Gstat->panhows[PIPE]);
#endif			
		        break;
		    }
		}
	    } /* while */
	}
    } else {
	/*
	 * jcol was EMPTY; Try to get a panel from the task Q.
	 */
    	while ( 1 ) {
    	    /*>>if ( (j = Dequeue(taskq, &item)) == EMPTY ) {*/
	    if ( taskq->count <= 0 ) {
		jcol = EMPTY;
		break;
	    } else {
		jcol = taskq->queue[taskq->head++];
		--taskq->count;
		if ( STATE( jcol ) >= CANGO ) { /* CANGO or CANPIPE */
#if ( DEBUGlevel>=1 )
		    printf("(%d) Dequeue[2] Got %d, STATE %d, Qcount %d\n",
			   pnum, jcol, STATE(jcol), j);
#endif
#ifdef PROFILE
		    if (STATE( jcol ) == CANGO) ++(Gstat->panhows[NOPIPE]);
		    else ++(Gstat->panhows[PIPE]);
#endif			
		    break;
		}
	    }
	} /* while */
    }
    
    /*
     * Update the status of the new panel "jcol" and its parent "dad".
     */
    if ( jcol != EMPTY ) {
	    --pxgstrf_shared->tasks_remain;
#ifdef DOMAINS
	if ( in_domain[jcol] == TREE_DOMAIN ) {
	    /* Dequeue the first descendant of this domain */
	    *bcol = taskq->queue[taskq->head++];
	    --taskq->count;
	} else
#endif
	{
	    STATE( jcol ) = BUSY;
	    w = pxgstrf_shared->pan_status[jcol].size;

	    for (j = jcol; j < jcol+w; ++j) pxgstrf_shared->spin_locks[j] = 1;
	    dad = DADPANEL (jcol);
	    if ( dad < n && pxgstrf_shared->pan_status[dad].ukids == 1 ) {
		STATE( dad ) = CANPIPE;
		/*>> j = Enqueue(taskq, dad);*/
		taskq->queue[taskq->tail++] = dad;
		++taskq->count;
#if ( DEBUGlevel>=1 )
		printf("(%d) Enqueue() %d's dad %d ->CANPIPE, Qcount %d\n",
		       pnum, jcol, dad, j);
#endif
	    }

#ifdef PROFILE
	    Gstat->procstat[pnum].panels++;
#endif
	
	    /* Find the farthest busy descendant of the new panel
	       and its parent.*/
	    *bcol = fb_cols[jcol];
#if ( DEBUGlevel>=1 )
	    printf("(%d) Scheduler[2] fb_cols[%d]=%d, STATE %d\n",
		   pnum, jcol, *bcol, STATE( *bcol ));
#endif
	    while ( STATE( *bcol ) == DONE ) *bcol = DADPANEL (*bcol);
	    fb_cols[dad] = *bcol;
	
	} /* else regular_panel */

    } /* if jcol != empty */

    *cur_pan = jcol;

#if ( DEBUGlevel>=1 )
    printf("(%d) Exit C.S. tasks_remain %d, cur_pan %d\n", 
	   pnum, pxgstrf_shared->tasks_remain, jcol);
#endif

} /* ---- END CRITICAL SECTION ---- */
    
#if ( MACH==SUN )
    /* Exit C.S. */
    mutex_unlock( &pxgstrf_shared->lu_locks[SCHED_LOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_unlock( &pxgstrf_shared->lu_locks[SCHED_LOCK] );
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard SCHED_LOCK
#endif    

#ifdef PROFILE
    Gstat->procstat[pnum].cs_time += SuperLU_timer_() - t;
#endif

    return;
}


/* @@@@@@@@@@@@@@ not called @@@@@@@@@@@@@@@@@@ */
#if 0
/* Fix the order of the panels to be taken. */
void
Preorder(const int_t pnum, const int_t n, const int_t *etree, int_t *cur_pan,
         queue_t *taskq, int_t *fb_cols, int_t *bcol,
	 pxgstrf_shared_t *pxgstrf_shared)
{
    register int_t w, dad, dad_ukids;

#undef POSTORDER
#ifdef POSTORDER
    if ( *cur_pan == EMPTY ) {
	*cur_pan = 0;
    } else {
	w = pxgstrf_shared->pan_status[*cur_pan].size;
	*cur_pan += w;
    }
#else /* Breadth-first bottom up */
    if ( *cur_pan != EMPTY ) {
	dad = DADPANEL (*cur_pan);
	dad_ukids = --pxgstrf_shared->pan_status[dad].ukids;
	if ( dad_ukids == 0 ) {
	    taskq->queue[taskq->tail++] = dad;
	    ++taskq->count;
	}
    }
    *cur_pan = taskq->queue[taskq->head++];
    --taskq->count;
#endif
    --pxgstrf_shared->tasks_remain;
    *bcol = *cur_pan;
}
#endif
/* @@@@@@@@@@@@@@ not called @@@@@@@@@@@@@@@@@@ */
