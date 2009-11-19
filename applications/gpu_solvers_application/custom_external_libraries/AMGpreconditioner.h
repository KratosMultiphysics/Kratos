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

#ifndef _AMGPRECONDITIONER_H
#define	_AMGPRECONDITIONER_H

#include "GPUPreconditioner.h"
#include "gpu_sparse.h"

using namespace Kratos::GPUSparse;


class AMGpreconditioner : public GPUPreconditioner {
public:
    AMGpreconditioner();
    AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, bool actAsPreconditioner);
    AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps, bool actAsPreconditioner);
    virtual ~AMGpreconditioner();
    void initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged);
    size_t solve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double precision, size_t maxIters);
    void singleStep(double* b_gpu, double* x_gpu);
    void cleanPreconditioner();
private:
    GPUCSRMatrix **Matrices;
    GPUCSRMatrix **P, **R, **G;
    GPUVector *b;
    size_t numFinalLevels;
    double W;
    double threshold;
    size_t numLevelsRoh;

    bool isPreconditioner;

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

