/* 
 * File:   Diagonalpreconditioner.h
 * Author: isaac
 *
 * Created on 4 / novembre / 2009, 11:26
 */

#ifndef _DIAGONALPRECONDITIONER_H
#define	_DIAGONALPRECONDITIONER_H

#include "GPUPreconditioner.h"
#include "gpu_sparse.h"

using namespace Kratos::GPUSparse;

class Diagonalpreconditioner : public GPUPreconditioner {
public:
    Diagonalpreconditioner();
    virtual ~Diagonalpreconditioner();

    void initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged);
    size_t solve(double* b_cpu, double* b_gpu, double* x_cpu, double* x_gpu, double precision, size_t iterations);
    void singleStep(double* b_gpu, double* x_gpu);
    void cleanPreconditioner() ;
private:
    _Vector G;
    bool allocated;
    double threshold;
    _Matrix A;
};

#endif	/* _DIAGONALPRECONDITIONER_H */

