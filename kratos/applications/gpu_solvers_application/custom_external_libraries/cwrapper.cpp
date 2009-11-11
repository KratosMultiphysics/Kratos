/*
==============================================================================
KratosGPUApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2009
Pooyan Dadvand, Riccardo Rossi, Isaac Gallego, Farshid Mossaiby 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
isaac.gallego.pla@gmail.com
mossaiby@yahoo.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#include "cwrapper.h"
#include "AMGpreconditioner.h"
#include "linear_solvers.h"
#include "Diagonalpreconditioner.h"
#include "mmio.h"
#include <stdio.h>

extern "C" {

static GPUPreconditioner* preconditioner;
static bool preconditionerSet = false;
static size_t* preSweeps;
static size_t* postSweeps;
static bool sweepsSet = false;


void AMGInitialize(double _W, size_t _numLevelsRoh, BOOLEAN _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t minimumSizeAllowed,
        size_t* _preSweeps, size_t* _postSweeps,
        size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, BOOLEAN dataIsChanged, BOOLEAN structureIsChanged, BOOLEAN actAsPreconditioner){
	if(preconditionerSet) AMGClean();
	preconditioner = new AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep == TRUE, _numMaxHierarchyLevels, minimumSizeAllowed, _preSweeps, _postSweeps, actAsPreconditioner == TRUE);
	preconditioner->initialize(ptr_cpu, indices_cpu, values_cpu,
		ptr_gpu, indices_gpu, values_gpu, numRows, numCols, numNNZ, dataIsChanged, structureIsChanged);
	preconditionerSet = true;
}
size_t AMGSolve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double precision, size_t maxIters){
	if(preconditionerSet)
    		return preconditioner->solve(b_gpu, b_cpu, x_gpu, x_cpu, precision, maxIters);
	return 0;
}
void AMGSingleStep(double* b_gpu, double* x_gpu){
	if(preconditionerSet)
    		preconditioner->singleStep(b_gpu, x_gpu);
}

void AMGClean(){
	preconditioner->cleanPreconditioner();
	delete preconditioner;
	preconditionerSet = false;
}

/** Solving strategies **/
void Solve(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double *mBNorm, double *mResidualNorm, size_t *mIterationsNumber, GPULinearSolverType a, GPUPreconditionerType b){
	if(preconditionerSet) AMGClean();
	
	if(b == AMG){
		size_t* preSweeps = new size_t[10];
		size_t* postSweeps = new size_t[10];
		for(size_t i = 0; i < 10; i++){
			preSweeps[i] = postSweeps[i] = 1;
		}
		preconditioner = new AMGpreconditioner(4.0/3.0, 2, true, 10, 100, preSweeps, postSweeps, true);
		preconditionerSet = true;
		sweepsSet = true;
	}else if(b == DIAGONAL){
		preconditioner = new Diagonalpreconditioner(true);
		preconditionerSet = true;
	}else if(b == NOPRECOND){
		preconditioner = NULL;
		preconditionerSet = false;
	}
	
	if(a == GPU_CG){
		CG_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, *mBNorm, *mResidualNorm, *mIterationsNumber, *preconditioner);
	}else if(a == GPU_BICGSTAB){
		BICGSTAB_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, *mBNorm, *mResidualNorm, *mIterationsNumber, *preconditioner);
	}
	//deletings
	if(sweepsSet){
		delete[] preSweeps;
		delete[] postSweeps;
		sweepsSet = false;
	}
	if(preconditionerSet){
		delete preconditioner;
		preconditionerSet = false;
	}
}

void Solve_AMGpreconditioner(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double *mBNorm, double *mResidualNorm, size_t *mIterationsNumber, GPULinearSolverType a, 
	double _W, size_t _numLevelsRoh, BOOLEAN _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps,
	BOOLEAN dataIsChanged, BOOLEAN structureIsChanged){
	if(preconditionerSet) AMGClean();
	preconditioner = new AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep == TRUE, _numMaxHierarchyLevels, _minimumSizeAllowed, _preSweeps, _postSweeps, true);
	preconditionerSet = true;
	if(a == GPU_CG){
		CG_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, *mBNorm, *mResidualNorm, *mIterationsNumber, *preconditioner);
	}else if(a == GPU_BICGSTAB){
		BICGSTAB_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, *mBNorm, *mResidualNorm, *mIterationsNumber, *preconditioner);
	}

	//deletings
	if(preconditionerSet){
		AMGClean();
	}
}

void loadDataFromFile(char* name, size_t ** indices_cpu, size_t **ptr_cpu, double **values_cpu,
	size_t *numRows, size_t *numCols, size_t *numNNZ){
    int ret_code;
    MM_typecode matcode;
    FILE *f;



    if ((f = fopen(name, "r")) == NULL)
        exit(1);

	int a;
    if ((a = mm_read_banner(f, &matcode)) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
	printf("%d\n", a);
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    int m, n, nnz;
    if ((ret_code = mm_read_mtx_crd_size(f, &n, &m, &nnz)) !=0)
        exit(1);

	

    *numCols = m;
    *numRows = n;
    *numNNZ = nnz;


    /* reseve memory for matrices */

    *indices_cpu = new size_t[*numNNZ];
    *values_cpu = new double[*numNNZ];
    *ptr_cpu = new size_t[*numRows+1];


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    size_t counter = 0;
    size_t currentRow, nextRow;
    currentRow = 1;
    nextRow = 1;
    (*ptr_cpu)[0] = 0;


    for (size_t i=0; i < *numNNZ; i++)
    {
        int position;
        double zz = 0.0;

        int err = fscanf(f, "%u %u %lg\n", &nextRow, &(*indices_cpu)[i], &zz);
	/*if(err != 1){
		printf("error reading mm\n");
		printf("%u, %u, %lg\n", nextRow, (*indices_cpu)[i], zz);
		exit(1);
	}*/
        position = (*indices_cpu)[i];
        (*values_cpu)[i] = zz;
        if(nextRow > currentRow){
            (*ptr_cpu)[currentRow] = counter + (*ptr_cpu)[currentRow-1];
            counter = 0;
            currentRow = nextRow;
        }
        (*indices_cpu)[i]--;  /* adjust from 1-based to 0-based */
        counter++;
    }
    (*ptr_cpu)[currentRow] = counter + (*ptr_cpu)[currentRow-1];

    fclose(f);
}

void loadDataFromFile_vector(char* name, double **values_cpu,
	size_t *numRows, size_t *numCols){
    int ret_code;
    MM_typecode matcode;
    FILE *f;



    if ((f = fopen(name, "r")) == NULL)
        exit(1);

	int a;
    if ((a = mm_read_banner(f, &matcode)) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
	printf("%d\n", a);
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    int m, n;
/*    if ((ret_code = mm_read_mtx_crd_size(f, &n, &m, &nnz)) !=0)
        exit(1);*/
	int err = fscanf(f, "%u %u \n", &n, &m);

	

    *numCols = m;
    *numRows = n;


    /* reseve memory for matrices */

    *values_cpu = new double[*numRows];


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    size_t counter = 0;
    size_t currentRow, nextRow;
    currentRow = 1;
    nextRow = 1;


    for (size_t i=0; i < *numRows; i++)
    {

        int err = fscanf(f, "%lg\n", &(*values_cpu)[i]);
	
    }

    fclose(f);
}

double* allocate_double(size_t size){
	return new double[size];
}

size_t* allocate_size_t(size_t size){
	return new size_t[size];
}

void deallocate_double(double* a){
	delete[] a;
}

void deallocate_size_t(size_t* a){
	delete[] a;
}

}
