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

#include "AMGpreconditioner.h"
#include <cstdio>

AMGpreconditioner::AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, size_t* _preSweeps, size_t* _postSweeps, bool actAsPreconditioner) {
        W = _W;
        numLevelsRoh = _numLevelsRoh;
	assumeZerosForEachStep = _assumeZerosForEachStep;
	numMaxHierarchyLevels = _numMaxHierarchyLevels;
	minimumSizeAllowed = _minimumSizeAllowed;

	preSweeps = _preSweeps;
	postSweeps = _postSweeps;
	isPreconditioner = actAsPreconditioner;
	numFinalLevels = 0;
        //printf("Minimum size allowed set in constructor: %lu\n", minimumSizeAllowed);
}

AMGpreconditioner::AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, bool actAsPreconditioner){
	W = _W;
        numLevelsRoh = _numLevelsRoh;
	assumeZerosForEachStep = _assumeZerosForEachStep;
	numMaxHierarchyLevels = _numMaxHierarchyLevels;
	minimumSizeAllowed = _minimumSizeAllowed;
	isPreconditioner = actAsPreconditioner;
	numFinalLevels = 0;
}

AMGpreconditioner::AMGpreconditioner(){
	numFinalLevels = 0;
}

AMGpreconditioner::~AMGpreconditioner() {

}

void AMGpreconditioner::cleanPreconditioner(){
	if(numFinalLevels > 0){
		//P
		delete[] P[0].indices_cpu;
		delete[] P[0].values_cpu;
		delete[] P[0].ptr_cpu;
		deletingStuff(P[0].indices_gpu);
		deletingStuff(P[0].values_gpu);
		deletingStuff(P[0].ptr_gpu);
		//R
		delete[] R[0].indices_cpu;
		delete[] R[0].values_cpu;
		delete[] R[0].ptr_cpu;
		deletingStuff(R[0].indices_gpu);
		deletingStuff(R[0].values_gpu);
		deletingStuff(R[0].ptr_gpu);
		//G
		delete[] G[0].indices_cpu;
		delete[] G[0].values_cpu;
		delete[] G[0].ptr_cpu;
		deletingStuff(G[0].indices_gpu);
		deletingStuff(G[0].values_gpu);
		deletingStuff(G[0].ptr_gpu);
		for(size_t i = 1; i < numFinalLevels; i++){
		    //P
		    delete[] P[i].indices_cpu;
		    delete[] P[i].values_cpu;
		    delete[] P[i].ptr_cpu;
		    deletingStuff(P[i].indices_gpu);
		    deletingStuff(P[i].values_gpu);
		    deletingStuff(P[i].ptr_gpu);
		    //R
		    delete[] R[i].indices_cpu;
		    delete[] R[i].values_cpu;
		    delete[] R[i].ptr_cpu;
		    deletingStuff(R[i].indices_gpu);
		    deletingStuff(R[i].values_gpu);
		    deletingStuff(R[i].ptr_gpu);
		    //G
		    delete[] G[i].indices_cpu;
		    delete[] G[i].values_cpu;
		    delete[] G[i].ptr_cpu;
		    deletingStuff(G[i].indices_gpu);
		    deletingStuff(G[i].values_gpu);
		    deletingStuff(G[i].ptr_gpu);
		    //A
		    delete[] Matrices[i].indices_cpu;
		    delete[] Matrices[i].values_cpu;
		    delete[] Matrices[i].ptr_cpu;
		    deletingStuff(Matrices[i].indices_gpu);
		    deletingStuff(Matrices[i].values_gpu);
		    deletingStuff(Matrices[i].ptr_gpu);
		}
    	}
	delete[] Matrices[numFinalLevels].matAuxValues;
	delete[] P;
	delete[] R;
	delete[] G;
	delete[] Matrices;
}

void AMGpreconditioner::initialize(size_t* ptr_cpu, size_t* indices_cpu, double* values_cpu,
        size_t* ptr_gpu, size_t* indices_gpu, double* values_gpu,
        size_t numRows, size_t numCols, size_t numNNZ, bool dataIsChanged, bool structureIsChanged){

	/*printf("PRINTING from AMGpreconditioner initialize, variable values:\n W = %f, Roh = %u, Zeros = %s, HierarchyLevels = %u, minimumSize = %u, firstPre = %u, secondPre = %u\n", W, numLevelsRoh, (assumeZerosForEachStep)?"true":"false", numMaxHierarchyLevels, minimumSizeAllowed, preSweeps[0], preSweeps[1]);*/

	//printf("El valor de minimumSizeAllowed es: %u\n", minimumSizeAllowed);

    Matrices = new _Matrix[numMaxHierarchyLevels];
    P = new _Matrix[numMaxHierarchyLevels];
    R = new _Matrix[numMaxHierarchyLevels];
    G = new _Matrix[numMaxHierarchyLevels];
    b.numElems = numCols;

    Matrices[0].numRows = numRows;
    Matrices[0].numCols = numCols;
    Matrices[0].numNNZ = numNNZ;
    Matrices[0].indices_cpu = indices_cpu;
    Matrices[0].ptr_cpu = ptr_cpu;
    Matrices[0].values_cpu = values_cpu;

    Matrices[0].indices_gpu = indices_gpu;
    Matrices[0].ptr_gpu = ptr_gpu;
    Matrices[0].values_gpu = values_gpu;

        /** Generating hierarchy **/
    numFinalLevels = generateHierarchy(Matrices, P, R, G, W, numLevelsRoh, numMaxHierarchyLevels, minimumSizeAllowed);
	printf("Initialize finalized, with numFinalLevels %u\n", numFinalLevels);
}

void AMGpreconditioner::singleStep(double* b_gpu, double* x_gpu){

	_Vector u;
	u.numElems = Matrices[0].numCols;
	u.values_cpu = new double[u.numElems];
	if(isPreconditioner){
		GPU_fillWithZeros(u.numElems, x_gpu);
	}

	u.values_gpu = x_gpu;
	b.values_gpu = b_gpu;
	b.values_cpu = new double[u.numElems];
	multilevel(Matrices, P, R, G, b, u, 0, numFinalLevels, preSweeps, postSweeps, assumeZerosForEachStep);
	delete[] u.values_cpu;
	delete[] b.values_cpu;
}

size_t AMGpreconditioner::solve(double* b_gpu, double* b_cpu, double* x_gpu, double* x_cpu, double _precision, size_t maxIters){
    threshold = _precision;
    _Vector u;
    double residual;
    u.numElems = Matrices[0].numCols;
    u.values_cpu = x_cpu;
    u.values_gpu = x_gpu;
    b.values_cpu = b_cpu;
    b.values_gpu = b_gpu;

    residual = checkResidual(u, b, Matrices[0]);

    bool done = false;
    size_t iterations = 0;

    while(!done && iterations < maxIters){
        iterations++;
        //multilevel(Matrices, P, R, G, b, u, 0, numFinalLevels);
        singleStep(b_gpu, x_gpu);        
        done = checkConvergence(u, b, Matrices[0], residual, threshold);
    }
    copyMem(u.values_gpu, u.values_cpu, u.numElems, 1);
    if(iterations == maxIters)
        printf("It haven't converged!!\n");
    else
        printf("Convergence achieved int %u iterations\n", iterations);

    return iterations;
}

void AMGpreconditioner::setPreSweeps( size_t* _preSweeps) {
	preSweeps = _preSweeps;
}

void AMGpreconditioner::setPostSweeps(size_t* _postSweeps) {
	postSweeps = _postSweeps;

}

