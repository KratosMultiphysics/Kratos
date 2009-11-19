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
    GPUVector *G;
    bool allocated;
    double threshold;
    GPUCSRMatrix *A;
    bool isPreconditioner;
};

#endif	/* _DIAGONALPRECONDITIONER_H */

