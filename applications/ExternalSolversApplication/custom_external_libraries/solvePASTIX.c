/*-----------------------------------------------------------------*
 * main test driver for the ARMS2 preconditioner for
 * Matrices in the COO/Harwell Boeing format
 *-----------------------------------------------------------------*
 * Yousef Saad - Aug. 2005.                                        *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu                 *
 *-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pastix.h"
#include "omp.h"

#define MPI_COMM_WORLD 0

void CompRow_to_CompCol(int m, int n, int nnz,
                        double *a, size_t *colind, size_t *rowptr,
                        pastix_float_t **at, pastix_int_t **rowind, pastix_int_t **colptr)
{
	register int i, j, col, relpos;
	pastix_int_t *marker = (pastix_int_t *)malloc(n*sizeof(pastix_int_t));
	for(i = 0; i<n; i++)
		marker[i] = 0;

	/* Allocate storage for another copy of the matrix. */
//    *at = (double *) (pastix_float_t *)malloc(nnz);
//    *rowind = (int *) (pastix_int_t *)malloc(nnz);
//    *colptr = (int *) (pastix_int_t *)malloc(n+1);

//printf("111\n");

	/* Get counts of each column of A, and set up column pointers */
	for(i = 0; i < m; ++i)
		for(j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
	(*colptr)[0] = 0;
	for(j = 0; j < n; ++j)
	{
		(*colptr)[j+1] = (*colptr)[j] + marker[j];
		marker[j] = (*colptr)[j];
	}
//printf("111\n");

	/* Transfer the matrix into the compressed column storage. */
	for(i = 0; i < m; ++i)
	{
		for(j = rowptr[i]; j < rowptr[i+1]; ++j)
		{
			col = colind[j];
			relpos = marker[col];
//printf("col %d \n",col);
//printf("relpos %d \n",relpos);
//printf("j %d \n",j);
//printf("i %d \n",i);
			(*rowind)[relpos] = i;
//printf("(*rowind)[relpos]  %d \n",(*rowind)[relpos] );
			(*at)[relpos] = a[j];
			++marker[col];
		}
	}
//printf("111\n");

	free(marker);
}


int solvePASTIX(int echo_level,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double *x, double* b)
{

	printf("000\n");
	pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
	pastix_int_t    ncol = mat_size;               /* Size of the matrix                                        */
	pastix_int_t   *rows      = (pastix_int_t *)malloc(sizeof(pastix_int_t)*nnz);  /* Indexes of first element of each column in row and values */
	pastix_int_t   *colptr        = (pastix_int_t *)malloc(sizeof(pastix_int_t)*(mat_size+1));   /* Row of each element of the matrix                         */
	pastix_float_t *values      = (pastix_float_t *)malloc(sizeof(pastix_float_t)*nnz);   /* Value of each element of the matrix                       */
	pastix_float_t *rhs         = (pastix_float_t *)malloc(sizeof(pastix_float_t)*mat_size);  /* right hand side                                           */
//	pastix_float_t *rhssaved    = NULL; /* right hand side (save)                                    */
//	pastix_float_t *ax          = NULL; /* A times X product                                         */
	pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
	double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
	pastix_int_t   *perm        = (pastix_int_t *)malloc(ncol*sizeof(pastix_int_t));; /* Permutation tabular                                       */
	pastix_int_t   *invp        = (pastix_int_t *)malloc(ncol*sizeof(pastix_int_t));; /* Reverse permutation tabular                               */
//		char           *type        = NULL; /* type of the matrix                                        */
//		char           *rhstype     = NULL; /* type of the right hand side                               */
//	pastix_int_t             mpid = 0;
//		driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
//	pastix_int_t             nbmatrices = 1;         /* Number of matrices given by user                          */
	pastix_int_t             nbthread = omp_get_max_threads();           /* Number of thread wanted by user                           */
	pastix_int_t             verbosemode = 2;        /* Level of verbose mode (0, 1, 2)                           */
	int             ordering = API_ORDER_SCOTCH;           /* Ordering to use                                           */
	pastix_int_t             nbrhs = 1;
	int             incomplete = 0;         /* Indicate if we want to use incomplete factorisation       */
	int             level_of_fill = 4;      /* Level of fill for incomplete factorisation                */
	int             amalgamation = 5;       /* Level of amalgamation for Kass                            */
	int             ooc = 2000;                /* OOC limit (Mo/percent depending on compilation options)   */
	pastix_int_t    mat_type = API_SYM_NO;
//		long            i;
//		double norme1, norme2;
	int i;
	printf("aaa\n");
	//compute the transpose
	CompRow_to_CompCol(mat_size, mat_size, nnz,AA, JA, IA,&values, &rows, &colptr);

//	/*
//	 * Matrix needs :
//	 *    - to be in fortran numbering
//	 *    - to have only the lower triangular part in symmetric case
//	 *    - to have a graph with a symmetric structure in unsymmetric case
//	 */
	if(NO_ERR != pastix_checkMatrix(0, verbosemode,
	                                mat_type,
	                                API_YES,
	                                ncol, &colptr, &rows, &values, NULL, 1))
		return 1;

	//copy b to the solution. It will be overwritten
	for(i = 0; i < mat_size; i++)
		x[i] = b[i];
	printf("bbb\n")	;
	/*******************************************/
	/* Initialize parameters to default values */
	/*******************************************/
	iparm[IPARM_MODIFY_PARAMETER] = API_NO;
	pastix(&pastix_data, MPI_COMM_WORLD,
	       ncol, colptr, rows, values,
	       perm, invp, x, 1, iparm, dparm);
	printf("ccc\n")	;

	/*******************************************/
	/*       Customize some parameters         */
	/*******************************************/
	iparm[IPARM_THREAD_NBR] = nbthread;
	iparm[IPARM_SYM] = mat_type;
	switch(mat_type)
	{
	case API_SYM_YES:
		iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
		break;
	case API_SYM_HER:
		iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
		break;
	default:
		iparm[IPARM_FACTORIZATION] = API_FACT_LU;
	}
	iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
	iparm[IPARM_VERBOSE]             = verbosemode;
	iparm[IPARM_ORDERING]            = ordering;
	iparm[IPARM_INCOMPLETE]          = incomplete;
	iparm[IPARM_OOC_LIMIT]           = ooc;

	if(incomplete == 1)
	{
		dparm[DPARM_EPSILON_REFINEMENT] = 1e-7;
	}
	iparm[IPARM_LEVEL_OF_FILL]       = level_of_fill;
	iparm[IPARM_AMALGAMATION_LEVEL]  = amalgamation;
	iparm[IPARM_RHS_MAKING]          = API_RHS_B;
	iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
	iparm[IPARM_END_TASK]            = API_TASK_CLEAN;




	/*******************************************/
	/*           Call pastix                   */
	/*******************************************/
	//note that we pass "x" instead of the rhs
	pastix(&pastix_data, 0,
	       ncol, colptr, rows, values,
	       perm, invp, x, nbrhs, iparm, dparm);

//		PRINT_RHS("SOL", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);
//		CHECK_SOL(rhs, rhssaved, ncol, mpid);




	free(colptr);
	free(rows);
	free(values);
	free(perm);
	free(invp);
	free(rhs);
//		free(rhssaved);
//		free(type);
//		free(rhstype);

	return EXIT_SUCCESS;
}
#undef MPI_COMM_WORLD
