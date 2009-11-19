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

#if !defined(LINEAR_GPU_SOLVERS_H_INCLUDED )
#define  LINEAR_GPU_SOLVERS_H_INCLUDED


#include <cmath>
#include <iostream>
#include "gpu_sparse.h"
#include "GPUPreconditioner.h"

using namespace Kratos::GPUSparse;

void BICGSTAB_GPU(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPUPreconditioner &preconditioner)
      {
	bool havePreconditioner = (&preconditioner != NULL);


	// Inputs
	GPUCSRMatrix gA(A_NNZ, A_Size1, A_size2, A_indices, A_ptr, A_values, false);
	GPU_CHECK(gA.GPU_Allocate());
	GPU_CHECK(gA.Copy(CPU_GPU, false));
	
	GPUVector gX(vectorSize, X_values);
	GPU_CHECK(gX.GPU_Allocate());
	GPU_CHECK(gX.Copy(CPU_GPU));
	
	GPUVector gB(vectorSize, B_values);
	GPU_CHECK(gB.GPU_Allocate());
	GPU_CHECK(gB.Copy(CPU_GPU));

	if(havePreconditioner){
		preconditioner.initialize(A_ptr, A_indices, A_values,
			gA.GPU_RowIndices, gA.GPU_Columns, gA.GPU_Values,
			A_Size1, A_size2, A_NNZ, true, true);
	}

	double resid;
	double roh_1 = 0.0, roh_2 = 0.0, alpha = 0.0, beta = 0.0, omega = 0.0;


	const int size = vectorSize;
	mBNorm = 0.0;
	double norms = 0.0;
	mResidualNorm = 0.0;

	GPUVector p(size);
	GPU_CHECK(p.GPU_Allocate());
	GPUVector phat(size);
	GPU_CHECK(phat.GPU_Allocate());
	GPUVector s(size);
	GPU_CHECK(s.GPU_Allocate());
	GPUVector shat(size);
	GPU_CHECK(shat.GPU_Allocate());
	GPUVector t(size);
	GPU_CHECK(t.GPU_Allocate());
	GPUVector v(size);
	GPU_CHECK(v.GPU_Allocate());


		// Real normb = norm(b);
	GPU_CHECK(GPUGPUVectorNorm2(gB, mBNorm));

		//  Vector r = b - A * x;
	//VectorType r(size);
	GPUVector r(size);
	GPU_CHECK(r.GPU_Allocate());
	//SparseSpaceType::Mult(rA,rX,r); // r = rA*rX
	GPU_CHECK(GPUGPUCSRMatrixVectorMultiply(gA, gX, r));

	//SparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r); // r = rB - r
	GPU_CHECK(GPUGPUVectorScaleAndAdd(1.00, gB, -1.00, r));
		// rtilde
	//VectorType rtilde(size);
	GPUVector rtilde(size);
	GPU_CHECK(rtilde.GPU_Allocate());
	GPU_CHECK(rtilde.CopyFromGPU(r));
	
	if(mBNorm == 0.0)
		mBNorm = 1.0;

	mIterationsNumber = 0;
	GPU_CHECK(GPUGPUVectorNorm2(r, mResidualNorm));
	if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
		if(havePreconditioner)
			preconditioner.cleanPreconditioner();
		return;
	}
	size_t i = 1;

	do{
		//roh_1 = dotProduct(rtilde, r);
		GPU_CHECK(GPUGPUVectorVectorMultiply(rtilde, r, roh_1));

		if(roh_1 == 0){
			GPU_CHECK(GPUGPUVectorNorm2(r, mResidualNorm));
			break;
		}

		if(i == 1){
			//p = r;
			GPU_CHECK(p.CopyFromGPU(r));
		}else{
			beta = (roh_1/roh_2) * (alpha/omega);
      			//p = r + beta * (p - omega * v);
			GPU_CHECK(GPUGPUVectorScaleAndAdd(-omega, v, 1.0, p));
			GPU_CHECK(GPUGPUVectorScaleAndAdd(1.0, r, beta, p));
		}
		if(havePreconditioner){
			//phat = M.solve(p);
			//GPU_fillWithZeros(size, phat.GPU_Values);
			preconditioner.singleStep(p.GPU_Values, phat.GPU_Values);
		}else{
			GPU_CHECK(phat.CopyFromGPU(p));
		}
		
		//v = A * phat;
		GPU_CHECK(GPUGPUCSRMatrixVectorMultiply(gA, phat, v));
		//dot(rtilde, v)
		double dotRtildeV = 0.0;
		GPU_CHECK(GPUGPUVectorVectorMultiply(rtilde, v, dotRtildeV));

		alpha = roh_1 / dotRtildeV;

		//s = r - alpha(0) * v;
		GPU_CHECK(GPUGPUVectorScaleAndAdd(1.0, r, -alpha, v, s));

		//norm(s)
		GPU_CHECK(GPUGPUVectorNorm2(s, norms));

		if ((resid = norms/mBNorm) < tol) {
			//x += alpha(0) * phat;
			//tol = resid;
			//return 0;
			break;
		}
		if(havePreconditioner){
			//shat = M.solve(s);
			//GPU_fillWithZeros(size, shat.GPU_Values);
			preconditioner.singleStep(s.GPU_Values, shat.GPU_Values);
		}else{
			GPU_CHECK(shat.CopyFromGPU(s));
		}

		//t = A * shat;
		GPU_CHECK(GPUGPUCSRMatrixVectorMultiply(gA, shat, t));

		//omega = dot(t,s) / dot(t,t);
		double dotTS = 0.0, dotTT = 0.0;
		GPU_CHECK(GPUGPUVectorVectorMultiply(t, s, dotTS));
		GPU_CHECK(GPUGPUVectorVectorMultiply(t, t, dotTT));
		omega = dotTS / dotTT;
		
		//x += alpha(0) * phat + omega(0) * shat;
		GPU_CHECK(GPUGPUVectorScaleAndAdd_addingVersion(alpha, phat, omega, shat, gX));
		//r = s - omega(0) * t;
		GPU_CHECK(GPUGPUVectorScaleAndAdd(1.0, s, -omega, t, r));

		roh_2 = roh_1;
		//norm(r)
		GPU_CHECK(GPUGPUVectorNorm2(r, mResidualNorm));
		mIterationsNumber++;
		if ((resid = mResidualNorm / mBNorm) < tol) {
			//tol = resid;
			//max_iter = i;
			//return 0;
			break;
		}
		if (omega == 0) {
			//tol = norm(r) / normb;
			//return 3;
			break;
    		}
		i++;
//		std::cout << "end iteration" << std::endl;
		//std::cout << mResidualNorm << std::endl;
//		std::cout << mBNorm << std::endl;
		//std::cout << mResidualNorm/mBNorm << std::endl;
	}while((mIterationsNumber < maxIterations) && (mResidualNorm > tol * mBNorm));
	GPU_CHECK(gX.Copy(GPU_CPU));
	if(havePreconditioner)
		preconditioner.cleanPreconditioner();
}

void CG_GPU(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPUPreconditioner& preconditioner)
      {

	bool havePreconditioner = (&preconditioner != 0);
	const int size = vectorSize;

	//Allocating matrix A
	GPUCSRMatrix gpuA(A_NNZ, A_Size1, A_size2, A_indices, A_ptr, A_values, false);
	GPU_CHECK(gpuA.GPU_Allocate());
	GPU_CHECK(gpuA.Copy(CPU_GPU, false));

	//Allocating vector b
	GPUVector gpuB(vectorSize, B_values);
	GPU_CHECK(gpuB.GPU_Allocate());
	GPU_CHECK(gpuB.Copy(CPU_GPU));

	//Allocating vector x
	GPUVector gpuX(vectorSize, X_values);
	GPU_CHECK(gpuX.GPU_Allocate());
	GPU_CHECK(gpuX.Copy(CPU_GPU));


	double alpha, beta, roh, roh_1 = 1.0;

	//Ini preconditioner
	//clock_t s1 = clock();
	if(havePreconditioner){
		preconditioner.initialize(A_ptr, A_indices, A_values,
				gpuA.GPU_RowIndices, gpuA.GPU_Columns, gpuA.GPU_Values,
				gpuA.Size1, gpuA.Size2, A_NNZ, true, true);
	}
	//clock_t s2 = clock();
	//std::cout << "Time to create hierarchy" << double(s2-s1) / CLOCKS_PER_SEC << "s" << std::endl;
	//Norm(b)
	GPU_CHECK(GPUGPUVectorNorm2(gpuB, mBNorm));

	//r = b - A*x
	GPUVector gpuR(size);
	GPU_CHECK(gpuR.GPU_Allocate());
	GPU_CHECK(GPUGPUCSRMatrixVectorMultiply(gpuA, gpuX, gpuR));
	GPU_CHECK(GPUGPUVectorScaleAndAdd(1.00, gpuB, -1.00, gpuR));

	if(mBNorm == 0.0)
		mBNorm = 1.0;

	mIterationsNumber = 0;
	//Norm(r)
	double resid;
	GPU_CHECK(GPUGPUVectorNorm2(gpuR, mResidualNorm));
	if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
		if(havePreconditioner)
			preconditioner.cleanPreconditioner();
		return;
	}
	GPUVector gpuP(size), gpuZ(size), gpuQ(size);
	GPU_CHECK(gpuP.GPU_Allocate());
	GPU_CHECK(gpuZ.GPU_Allocate());
	GPU_CHECK(gpuQ.GPU_Allocate());
	
	//clock_t s3 = 0.0;
	size_t i = 1;
	do{
		//s1 = clock();
		if(havePreconditioner){	
			preconditioner.singleStep(gpuR.GPU_Values, gpuZ.GPU_Values);	
		}else{
			gpuZ.CopyFromGPU(gpuR);
		}
		//s2 = clock();
		//s3 += s2-s1;
		//roh = GPU_dotProduct(gpuR.Size, gpuR.GPU_Values, 1, gpuZ.GPU_Values, 1);
		GPU_CHECK(GPUGPUVectorVectorMultiply(gpuR, gpuZ, roh));

		if(i == 1)
			//gpuP = gpuZ
			gpuP.CopyFromGPU(gpuZ);
		else{
			beta = roh / roh_1; 
			GPU_CHECK(GPUGPUVectorScaleAndAdd(1.00, gpuZ, beta, gpuP));
		}

		GPU_CHECK(GPUGPUCSRMatrixVectorMultiply(gpuA, gpuP, gpuQ));
		//roh_1 = dotProduct(p, q);
		double dotPQ;
		GPU_CHECK(GPUGPUVectorVectorMultiply(gpuP, gpuQ, dotPQ));
		alpha = roh / dotPQ;
		GPU_CHECK(GPUGPUVectorScaleAndAdd(alpha, gpuP, 1.00, gpuX));
		GPU_CHECK(GPUGPUVectorScaleAndAdd(-alpha, gpuQ, 1.00, gpuR));
	
		GPU_CHECK(GPUGPUVectorNorm2(gpuR, mResidualNorm));
		//std::cout << mResidualNorm << std::endl;
		//std::cout << mResidualNorm/mBNorm << std::endl;
		if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
			break;
		}
		roh_1 = roh;
		mIterationsNumber++;			
		i++;
	}while((mIterationsNumber < maxIterations) && (mResidualNorm > tol * mBNorm));
	//std::cout << "Average time for single step" << (double(s3)/(i-1)) / CLOCKS_PER_SEC << "s" << std::endl;
	//std::cout << "the final RESIDUAL NORM is " << mResidualNorm << std::endl;
	GPU_CHECK(gpuX.Copy(GPU_CPU));
	if(havePreconditioner)
		preconditioner.cleanPreconditioner();
}

//to add in gpu_sparse.cu

#endif
