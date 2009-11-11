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

#include "Diagonalpreconditioner.h"
#include <cstdio>


Diagonalpreconditioner::Diagonalpreconditioner(bool actAsPreconditioner) {
	isPreconditioner = actAsPreconditioner;
}

Diagonalpreconditioner::~Diagonalpreconditioner() {
    if(allocated)
        cleanPreconditioner();
}

void Diagonalpreconditioner::initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu, size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu, size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged){
    //Sacar la diagonal y guardarla en cpu-gpu
    allocated = true;
    //aqui gpu
    A.numRows = numRows;
    A.numCols = numCols;
    A.numNNZ = numNNZ;
    A.ptr_cpu = ptr_cpu;
    A.indices_cpu = indices_cpu;
    A.values_cpu = values_cpu;
    A.ptr_gpu = ptr_gpu;
    A.indices_gpu = indices_gpu;
    A.values_gpu = values_gpu;
    createDiagonal(A, G);
}

void Diagonalpreconditioner::cleanPreconditioner(){
    //borrar los dos vectores de la diagonal
    deletingStuff(G.values_gpu);
    delete[] G.values_cpu;
    allocated = false;
}

void Diagonalpreconditioner::singleStep(double* b_gpu, double* x_gpu){
    //aplicar una multiplicación vector vector, de nuestra diagonal a su vector
	if(isPreconditioner)
		copyMem(b_gpu, x_gpu, G.numElems, 2);
    GPU_VectorMultiply(G.values_gpu, x_gpu, G.numElems);
}

size_t Diagonalpreconditioner::solve(double* b_cpu, double* b_gpu, double* x_cpu, double* x_gpu, double _precision, size_t maxIters){
    threshold = _precision;
    _Vector u;
    double residual;
    u.numElems = G.numElems;
    u.values_cpu = x_cpu;
    u.values_gpu = x_gpu;
    _Vector b;
    b.numElems = G.numElems;
    b.values_cpu = b_cpu;
    b.values_gpu = b_gpu;

    residual = checkResidual(u, b, A);

    bool done = false;
    size_t iterations = 0;

    while(!done && iterations < maxIters){
        iterations++;
        //multilevel(Matrices, P, R, G, b, u, 0, numFinalLevels);
        singleStep(b_gpu, x_gpu);
        copyMem(u.values_gpu, u.values_cpu, u.numElems, 1);
        done = checkConvergence(u, b, A, residual, threshold);
    }
    if(iterations == maxIters)
        printf("It haven't converged!!\n");
    else
        printf("Convergence achieved int %u iterations\n", iterations);

    return iterations;
}
