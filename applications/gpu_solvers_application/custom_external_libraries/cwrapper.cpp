#include "cwrapper.h"
#include "AMGpreconditioner.h"
#include "linear_solvers.h"
#include "Diagonalpreconditioner.h"

GPUPreconditioner* preconditioner;
bool preconditionerSet = false;
size_t* preSweeps;
size_t* postSweeps;
bool sweepsSet = false;

void AMGInitialize(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t minimumSizeAllowed,
        size_t* _preSweeps, size_t* _postSweeps,
        size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged){
	if(preconditionerSet) AMGClean();
	preconditioner = new AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep, _numMaxHierarchyLevels, minimumSizeAllowed, _preSweeps, _postSweeps);
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
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPULinearSolverType a, GPUPreconditionerType b){
	if(preconditionerSet) AMGClean();
	
	if(b == AMG){
		size_t* preSweeps = new size_t[10];
		size_t* postSweeps = new size_t[10];
		preconditioner = new AMGpreconditioner(4.0/3.0, 2, true, 10, 1000, preSweeps, postSweeps);
		preconditionerSet = true;
		sweepsSet = true;
	}else if(b == DIAGONAL){
		preconditioner = new Diagonalpreconditioner();
		preconditionerSet = true;
	}else if(b == NOPRECOND){
		preconditioner = NULL;
	}
	
	if(a == GPU_CG){
		CG_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, mBNorm, mResidualNorm, mIterationsNumber, *preconditioner);
	}else if(a == GPU_BICGSTAB){
		BICGSTAB_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, mBNorm, mResidualNorm, mIterationsNumber, *preconditioner);
	}

	//deletings
	if(sweepsSet){
		delete[] preSweeps;
		delete[] postSweeps;
	}
	if(preconditionerSet){
		AMGClean();
	}
}

void Solve_AMGpreconditioner(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPULinearSolverType a, 
	double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps,
	bool dataIsChanged, bool structureIsChanged){
	if(preconditionerSet) AMGClean();
	preconditioner = new AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep, _numMaxHierarchyLevels, _minimumSizeAllowed, _preSweeps, _postSweeps);
	preconditionerSet = true;
	if(a == GPU_CG){
		CG_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, mBNorm, mResidualNorm, mIterationsNumber, *preconditioner);
	}else if(a == GPU_BICGSTAB){
		BICGSTAB_GPU(A_Size1, A_size2, A_NNZ, A_values, A_indices, A_ptr, vectorSize, X_values, B_values, tol, maxIterations, mBNorm, mResidualNorm, mIterationsNumber, *preconditioner);
	}

	//deletings
	if(preconditionerSet){
		AMGClean();
	}
}
