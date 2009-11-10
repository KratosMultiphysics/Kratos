/* 
 * File:   Diagonalpreconditioner.cpp
 * Author: isaac
 * 
 * Created on 4 / novembre / 2009, 11:26
 */

#include "Diagonalpreconditioner.h"
#include <cstdio>


Diagonalpreconditioner::Diagonalpreconditioner() {
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
