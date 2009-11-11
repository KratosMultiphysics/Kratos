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

#ifndef _CWRAPPER_H
#define	_CWRAPPER_H

#include <stdlib.h>

typedef int BOOLEAN;

#define FALSE 0
#define TRUE 1



typedef enum GPULinearSolverType{
	GPU_CG,
	GPU_BICGSTAB
}GPULinearSolverType;
typedef enum GPUPreconditionerType{
	AMG,
	DIAGONAL,
	NOPRECOND
}GPUPreconditionerType;

#ifdef __cplusplus
extern "C" {
#endif

/** For AMG use **/
void AMGInitialize(double _W, size_t _numLevelsRoh, BOOLEAN _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t minimumSizeAllowed,
        size_t* _preSweeps, size_t* _postSweeps,
        size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, BOOLEAN dataIsChanged, BOOLEAN structureIsChanged, BOOLEAN actAsPreconditioner);
size_t AMGSolve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double precision, size_t maxIters);
void AMGSingleStep(double* b_gpu, double* x_gpu);
void AMGClean();

/** Solving strategies **/
void Solve(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double *mBNorm, double *mResidualNorm, size_t *mIterationsNumber, GPULinearSolverType , GPUPreconditionerType);

void Solve_AMGpreconditioner(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double *mBNorm, double *mResidualNorm, size_t *mIterationsNumber, GPULinearSolverType, 
	double _W, size_t _numLevelsRoh, BOOLEAN _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps,
	BOOLEAN dataIsChanged, BOOLEAN structureIsChanged);

/** Loading matrix market **/
void loadDataFromFile(char* name, size_t ** indices_cpu, size_t **ptr_cpu, double **values_cpu,
	size_t *numRows, size_t *numCols, size_t *numNNZ);
void loadDataFromFile_vector(char* name, double **values_cpu,
	size_t *numRows, size_t *numCols);

double* allocate_double(size_t size);
size_t* allocate_size_t(size_t size);
void deallocate_double(double* a);
void deallocate_size_t(size_t* a);

#ifdef __cplusplus
}
#endif

#endif	/* _CWRAPPER_H */

