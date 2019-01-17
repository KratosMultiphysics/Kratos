/**
 *
 * @file order_find_supernodes.c
 *
 * PaStiX order function to find supernodes out of a given permutation.
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
#include "graph.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief Computes the size of each subtree.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of nodes.
 *
 * @param[in] father
 *          Array of size n.
 *          List of father to each node. If node i is a root then father[i] = i.
 *
 * @param[in] iperm
 *          Array of size n. The inverse permutation vector.
 *
 * @param[out] T
 *          Array of size n.
 *          On exit, size of each subtree.
 *
 *******************************************************************************/
static inline void
compute_subtree_size(      pastix_int_t  n,
                     const pastix_int_t *father,
                     const pastix_int_t *iperm,
                           pastix_int_t *T)
{
    /********************************************/
    /*  Compute the size of each subtree given  */
    /*  the number of the father of each node   */
    /********************************************/
    pastix_int_t k, i;

    /* TODO pas la peine d'utiliser un tas; il suffit de parcourir iperm pour assurer
     de toujours traiter un fils avant son pere */

    bzero(T, sizeof(pastix_int_t)*n);

    for(k=0;k<n;k++)
    {
        i = iperm[k];
        T[i]++;
        if(i!=father[i])
            T[father[i]] += T[i];
    }

#if defined(PASTIX_DEBUG_ORDERING)
    {
        pastix_int_t sum = 0;

        for(i=0;i<n;i++) {
            if(father[i] == i)
                sum += T[i];
        }
        if(sum != n)
            errorPrint("compute_subtree_size: sum of the subtree = %ld n = %ld", (long)sum, (long)n);
        assert(n == sum);
    }
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief Computes the post order of the elimination tree given on
 * entry.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of nodes.
 *
 * @param[in] father
 *          Array of size n.
 *          List of father to each node. If node i is a root then father[i] = i.
 *
 * @param[inout] perm
 *          Array of size n. The permutation vector.
 *          On exit, the postorder permutation vector.
 *
 * @param[inout] invp
 *          Array of size n. The inverse permutation vector.
 *          On exit, the postorder inverse permutation vector.
 *
 * @param[inout] T
 *          Workspace of size n.
 *
 *******************************************************************************/
static inline void
compute_post_order(      pastix_int_t n,
                   const pastix_int_t *father,
                         pastix_int_t *perm,
                         pastix_int_t *invp,
                         pastix_int_t *T)
{
    pastix_int_t i;
    pastix_int_t j, k, t;

    /*
     * First compute the number of node in the subtree rooted in node i
     */
    compute_subtree_size(n, father, invp, T);

    /*
     * When multiple roots are present we have to compute the start index of
     * each root
     */
    t = 0;
    for(k=0;k<n;k++) {
        i = invp[k];
        if(father[i] == i) {
            /* This is a root */
            j = T[i];
            T[i] += t;
            t += j;
        }
    }

#if defined(PASTIX_DEBUG_ORDERING)
    for(i=0;i<n;i++)
        assert(T[i] <= T[father[i]]);
#endif

    for(k=n-1;k>=0;k--)
    {
        i = invp[k];
        perm[i] = T[father[i]]; /* We MUST HAVE father[i] == i for a root ! */
        T[father[i]] -= T[i];
        T[i] = perm[i]-1;
        assert(perm[father[i]] >= perm[i]);
    }

    /*
     * We need to retrieve 1 for the C numbering compliance
     */
    for(i=0;i<n;i++)
        perm[i]--;

#if defined(PASTIX_DEBUG_ORDERING)
    /* Check the permutation vector */
    for(i=0;i<n;i++)
    {
        assert(perm[i] >= 0);
        assert(perm[i] <  n);
    }

    bzero(invp, sizeof(pastix_int_t)*n);
    for(i=0;i<n;i++)
        invp[perm[i]]++;

    k = 0;
    for(i=0;i<n;i++)
        if(invp[i] != 1)
            k++;
    if(k>0)
        errorPrint("Number of errors in perm vector in postorder %ld", (long)k);
    assert(k==0);
#endif

    /*
     * Compute the invp vector
     */
    for(i=0; i<n; i++)
        invp[perm[i]] = i;
}

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief Compute the elimination tree of a matrix A
 * (without computing the symbolic factorization) associated with a reordering
 * of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of vertices in the given CSC
 *
 * @param[in] ia
 *          Array of size n+1.
 *          The column pointer tabular for the given CSC. (In C numbering)
 *
 * @param[in] ja
 *          Array of size ia[n].
 *          The list of edges in the given CSC (In C numbering)
 *
 * @param[in] perm
 *          Array of size n. The permutation vector.
 *
 * @param[in] invp
 *          Array of size n. The inverse permutation vector.
 *
 * @param[inout] father
 *          Array of size n.
 *          On entry, an allocated array of size n.
 *          On exit, father[i] = father of ith node on the eliminination
 *          tree. If node i is a root then father[i] = i.
 *
 *******************************************************************************/
static inline void
compute_elimination_tree(      pastix_int_t n,
                         const pastix_int_t *ia,
                         const pastix_int_t *ja,
                         const pastix_int_t *perm,
                         const pastix_int_t *invp,
                               pastix_int_t *father)
{
    pastix_int_t i, j, k;
    pastix_int_t node;
    pastix_int_t vroot;

    /* Optim */
    pastix_int_t flag, ind;
    pastix_int_t *jrev = NULL;
    pastix_int_t *jf   = NULL;

    MALLOC_INTERN(jrev, n, pastix_int_t);
    for(i=0;i<n;i++)
        jrev[i] = -1;

    MALLOC_INTERN(jf, n, pastix_int_t);
    bzero(jf, sizeof(pastix_int_t)*n);

    for(i=0;i<n;i++)
        father[i] = -1;

    for(i=0;i<n;i++)
    {
        ind = 0;
        node = invp[i];
        for(j=ia[node];j<ia[node+1];j++)
        {

            k = ja[j];
            if(perm[k] < perm[node])
            {
                flag = 1;
                vroot = k;
                while(father[vroot] != -1 && father[vroot] != node)
                {
                    if(jrev[vroot] >= 0)
                    {
                        flag = 0;
                        break;
                    }
                    jrev[vroot] = ind;
                    jf[ind] = vroot;
                    ind++;

                    vroot = father[vroot];
                }
                if(flag == 1)
                    father[vroot] = node;
            }
        }
        /* reinit jrev */
        for(j=0;j<ind;j++)
            jrev[jf[j]]=-1;
    }

    memFree_null(jrev);
    memFree_null(jf);

    for(i=0;i<n;i++)
        if(father[i] == -1)
            father[i] = i;

#if defined(PASTIX_DEBUG_ORDERING)
    /* Check to see if a father has a lower rank in the permutation array than one of its sons */
    for(i=0;i<n;i++)
    {
        if(perm[i] > perm[father[i]])
        {
            fprintf(stderr, "Node %ld perm=%ld Father %ld perm=%ld \n",
                    (long)i, (long)perm[i], (long)father[i], (long)perm[father[i]]);
            assert(perm[i] <= perm[father[i]]);
        }
    }
#endif
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Computes the set of supernodes for a given permutation.
 *
 * The permutation of the matrix given on entry is modified to obtain a
 * postorder of the elimination tree: this does not affect the fill-in properties
 * of the initial ordering.
 *
 * WARNING: The matrix pattern is assumed to be symmetric.
 *
 * NOTE: This function can take on entry the lower triangular part or the whole
 * matrix A.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          The graph associated with the order structuree on which we need to
 *          find the supernodes.
 *
 * @param[inout] ordeptr
 *          Pointer to a pastix_order_t structure, that will be further initialized by
 *          the routine.
 *
 *          On entry:
 *            - orderptr->permtab: the original permutation vector for the
 *                elimination tree.
 *            - orderptr->permtab: the original inverse permutation vector for the
 *                elimination tree.
 *
 *          On exit:
 *            - orderptr->cblknbr: The number of supernodes found.
 *            - orderptr->rangtab: Contains the first element of each
 *                supernode. rangtab[i] is the first element of the ith
 *                supernode.
 *            - orderptr->permtab: the permutation vector for the
 *                postorder of the nodes in the elimination tree.
 *            - orderptr->permtab: the inverse permutation vector for the
 *                postorder of the nodes in the elimination tree.
 *            - orderptr->treetab: treetab[i] is the number of the father of
 *                supernode i in the supernodal elimination tree.
 *
 *******************************************************************************/
void
pastixOrderFindSupernodes( const pastix_graph_t *graph,
                           pastix_order_t * const ordeptr )
{
    pastix_int_t *father     = NULL; /* father[i] is the father of node i in he elimination tree of A */
    pastix_int_t *T;                 /* T[j] is the number of node in the subtree rooted in node j in
                                      the elimination tree of A */
    pastix_int_t *S          = NULL; /* S[i] is the number of sons for node i in the elimination tree */
    pastix_int_t *isleaf     = NULL;
    pastix_int_t *prev_rownz = NULL;
    pastix_int_t *treetab    = NULL;
    pastix_int_t i, j, k;
    pastix_int_t pi, pj;

    pastix_int_t  snodenbr;
    pastix_int_t *perm = ordeptr->permtab;
    pastix_int_t *invp = ordeptr->peritab;
    pastix_int_t  n  = graph->n;
    pastix_int_t *ia = graph->colptr;
    pastix_int_t *ja = graph->rows;

    assert( graph->colptr[0] == 0 );

    /* Free the rangtab/treetab array if existing */
    if ( ordeptr->rangtab != NULL ) {
        memFree_null( ordeptr->rangtab );
    }
    if ( ordeptr->treetab != NULL ) {
        memFree_null( ordeptr->treetab );
    }

    MALLOC_INTERN(S,      n, pastix_int_t);
    MALLOC_INTERN(father, n, pastix_int_t);

#if defined(PASTIX_DEBUG_ORDERING)
    /* Check that the permutation vector is 0 based */
    for(i=0;i<n;i++) {
        assert(perm[i] >= 0);
        assert(perm[i] < n);
        assert(invp[i] >= 0);
        assert(invp[i] < n);
    }

    bzero(S, sizeof(pastix_int_t)*n);
    for(i=0;i<n;i++)
        S[perm[i]]++;

    k = 0;
    for(i=0;i<n;i++) {
        if(S[i] != 1) {
            k++;
        }
    }

    if(k>0) {
        errorPrint("perm array is not valid, number of error =  %ld", (long)k);
    }
    assert(k==0);
#endif

    /*
     * Compute the elimination tree of A
     */
    compute_elimination_tree(n, ia, ja, perm, invp, father);

    /*
     * Compute the postorder of the elimination tree
     * Warning: This operation modifies perm and invp
     */
    MALLOC_INTERN(T, n+1, pastix_int_t);
    compute_post_order(n, father, perm, invp, T);

    /*
     * Compute the number of descendant of each node i in the elimination tree
     */
    compute_subtree_size(n, father, invp, T);

    MALLOC_INTERN(isleaf,     n, pastix_int_t);
    MALLOC_INTERN(prev_rownz, n, pastix_int_t);
    bzero(isleaf,     sizeof(pastix_int_t)*n);
    bzero(prev_rownz, sizeof(pastix_int_t)*n);

    for(j=0;j<n;j++)
    {
        pj = invp[j];
        for(i=ia[pj];i<ia[pj+1];i++)
        {
            pi = perm[ja[i]];
            if(pi > j)
            {
                k = prev_rownz[pi];
                if(k < j - T[pj]+1 ) {
                    isleaf[j] = 1;
                }

                prev_rownz[pi] = j;
            }
        }
    }
    memFree(prev_rownz);

    /*
     * Compute the number of sons of each node in the elimination tree
     * The snodetab/rangtab is computed in the workspace T.
     */
    bzero(S, sizeof(pastix_int_t)*n);
    for(i=0;i<n;i++) {
        if(father[i] != i) {
            S[father[i]]++;
        }
    }

    for(i=0;i<n;i++) {
        if(S[i] != 1) {
            isleaf[perm[i]] = 1;
        }
    }

    snodenbr = 0;
    for(i=0;i<n;i++) {
        if(isleaf[i] == 1)
        {
            T[snodenbr] = i;
            snodenbr++;
        }
    }
    T[snodenbr] = n;

    memFree(isleaf);

    /*
     * If the treetab is required, we compute it before to free the local data.
     */
    assert( snodenbr > 0 );
    MALLOC_INTERN( ordeptr->treetab, snodenbr, pastix_int_t );
    treetab = ordeptr->treetab;
    {
        pastix_int_t dad;

        /* Node to supernode conversion vector */
        for(i=0;i<snodenbr;i++) {
            for(j=T[i];j<T[i+1];j++) {
                S[j] = i;
            }
        }

        /* Fill the treetab info */
        for(i=0;i<snodenbr;i++)
        {
            k=snodenbr;
            for(j=T[i];j<T[i+1];j++)
            {
                dad = S[perm[father[invp[j]]]];
                if( dad < k && dad > i) {
                    k = dad;
                }
            }
            treetab[i] = k;
            if(k == snodenbr)
            {
                treetab[i] = -1; /* This is a root */
            }
            assert((treetab[i] == -1) || (treetab[i] >= i));
        }
    }

    memFree(father);
    memFree(S);

    /*
     * Save result to ordeptr (switch former rangtab and new one allocated to
     * restricted array)
     */
    ordeptr->cblknbr = snodenbr;
    MALLOC_INTERN( ordeptr->rangtab, snodenbr+1, pastix_int_t );
    memcpy( ordeptr->rangtab, T, (snodenbr+1)*sizeof(pastix_int_t) );
    memFree(T);

    /* Check that the order is 0 based */
    assert(ordeptr->rangtab[0] == 0);
}
