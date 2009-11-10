/* 
 * File:   cwrapper.h
 * Author: isaac
 *
 * Created on 20 / octubre / 2009, 09:03
 */

#ifndef _CWRAPPER_H
#define	_CWRAPPER_H

#include <cstdlib>

enum GPULinearSolverType{
	GPU_CG,
	GPU_BICGSTAB,
};
enum GPUPreconditionerType{
	AMG,
	DIAGONAL,
	NOPRECOND
};

/** For AMG use **/
void AMGInitialize(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t minimumSizeAllowed,
        size_t* _preSweeps, size_t* _postSweeps,
        size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged);
size_t AMGSolve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double precision, size_t maxIters);
void AMGSingleStep(double* b_gpu, double* x_gpu);
void AMGClean();

/** Solving strategies **/
void Solve(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPULinearSolverType, GPUPreconditionerType);

void Solve_AMGpreconditioner(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPULinearSolverType, 
	double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps,
	bool dataIsChanged, bool structureIsChanged);

#endif	/* _CWRAPPER_H */

