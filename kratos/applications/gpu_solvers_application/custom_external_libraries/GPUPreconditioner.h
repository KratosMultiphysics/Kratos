/* 
 * File:   preconditioner.h
 * Author: isaac
 *
 * Created on 1 / octubre / 2009, 09:31
 */

#ifndef _GPUPRECONDITIONER_H
#define	_GPUPRECONDITIONER_H

#include <cstdlib>

class GPUPreconditioner {
public:
    GPUPreconditioner() {};
    virtual ~GPUPreconditioner() {};

    virtual void initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged){};
    virtual size_t solve(double* b_cpu, double* b_gpu, double* x_cpu, double* x_gpu, double precision, size_t maxIters){return 0;}
    virtual void singleStep(double* b_gpu, double* x_gpu){};
	virtual void cleanPreconditioner() {};

private:

};

#endif	/* _PRECONDITIONER_H */

