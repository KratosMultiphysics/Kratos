/* 
 * File:   AMGpreconditioner.h
 * Author: isaac
 *
 * Created on 1 / octubre / 2009, 09:39
 */

#ifndef _AMGPRECONDITIONER_H
#define	_AMGPRECONDITIONER_H

#include "GPUPreconditioner.h"
#include "gpu_sparse.h"

using namespace Kratos::GPUSparse;


class AMGpreconditioner : public GPUPreconditioner {
public:
    AMGpreconditioner();
    AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed);
    AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps);
    virtual ~AMGpreconditioner();
    void initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged);
    size_t solve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double precision, size_t maxIters);
    void singleStep(double* b_gpu, double* x_gpu);
    void cleanPreconditioner();
private:
    _Matrix* Matrices;
    _Matrix *P, *R, *G;
    _Vector b;
    size_t numFinalLevels;
    double W;
    double threshold;
    size_t numLevelsRoh;

    size_t* preSweeps;
    size_t* postSweeps;

    bool assumeZerosForEachStep;//this variable choose between solver (false) or singlestep (true)
    size_t numMaxHierarchyLevels;
    size_t minimumSizeAllowed;

    size_t counter;

protected:

	void setPreSweeps(size_t* _preSweeps);
	void setPostSweeps(size_t* _postSweeps);
};

#endif	/* _AMGPRECONDITIONER_H */

