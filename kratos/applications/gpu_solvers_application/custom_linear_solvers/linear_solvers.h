#if !defined(LINEAR_GPU_SOLVERS_H_INCLUDED )
#define  LINEAR_GPU_SOLVERS_H_INCLUDED


#include <cmath>
#include <iostream>
#include "gpu_sparse.h"
#include "GPUPreconditioner.h"

#define KRATOS_GPU_CHECK(call)		if (!(call)) printf(#call " in file %s line %d failed.\n", __FILE__, __LINE__);
using namespace Kratos::GPUSparse;

void BICGSTAB_GPU(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPUPreconditioner &preconditioner)
      {
	bool havePreconditioner = (&preconditioner != NULL);

	// Inputs
	GPUCSRMatrix gA(A_NNZ, A_Size1, A_size2, A_indices, A_ptr, A_values);
	KRATOS_GPU_CHECK(gA.GPU_Allocate());
	KRATOS_GPU_CHECK(gA.Copy(CPU_GPU, false));
	
	GPUVector gX(vectorSize, X_values);
	KRATOS_GPU_CHECK(gX.GPU_Allocate());
	KRATOS_GPU_CHECK(gX.Copy(CPU_GPU));
	
	GPUVector gB(vectorSize, B_values);
	KRATOS_GPU_CHECK(gB.GPU_Allocate());
	KRATOS_GPU_CHECK(gB.Copy(CPU_GPU));

	if(havePreconditioner){
		preconditioner.initialize(gA.CPU_RowIndices, gA.CPU_Columns, gA.CPU_Values,
			gA.GPU_RowIndices, gA.GPU_Columns, gA.GPU_Values,
			gA.Size1, gA.Size2, gA.NNZ, true, true);
	}

	double resid;
	double roh_1 = 0.0, roh_2 = 0.0, alpha = 0.0, beta = 0.0, omega = 0.0;


	const int size = vectorSize;
	mBNorm = 0.0;
	double norms = 0.0;
	mResidualNorm = 0.0;

	GPUVector p(size);
	KRATOS_GPU_CHECK(p.GPU_Allocate());
	GPUVector phat(size);
	KRATOS_GPU_CHECK(phat.GPU_Allocate());
	GPUVector s(size);
	KRATOS_GPU_CHECK(s.GPU_Allocate());
	GPUVector shat(size);
	KRATOS_GPU_CHECK(shat.GPU_Allocate());
	GPUVector t(size);
	KRATOS_GPU_CHECK(t.GPU_Allocate());
	GPUVector v(size);
	KRATOS_GPU_CHECK(v.GPU_Allocate());


		// Real normb = norm(b);
	KRATOS_GPU_CHECK(GPU_VectorNorm2(gB, mBNorm));

		//  Vector r = b - A * x;
	//VectorType r(size);
	GPUVector r(size);
	KRATOS_GPU_CHECK(r.GPU_Allocate());
	//SparseSpaceType::Mult(rA,rX,r); // r = rA*rX
	KRATOS_GPU_CHECK(GPU_MatrixVectorMultiply(gA, gX, r));

	//SparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r); // r = rB - r
	KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gB, -1.00, r));
		// rtilde
	//VectorType rtilde(size);
	GPUVector rtilde(size);
	KRATOS_GPU_CHECK(rtilde.GPU_Allocate());
	KRATOS_GPU_CHECK(rtilde.CopyFromGPU(r));
	
	if(mBNorm == 0.0)
		mBNorm = 1.0;

	mIterationsNumber = 0;
	KRATOS_GPU_CHECK(GPU_VectorNorm2(r, mResidualNorm));
	if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
		if(havePreconditioner)
			preconditioner.cleanPreconditioner();
		return;
	}
	size_t i = 1;

	do{
		//roh_1 = dotProduct(rtilde, r);
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(rtilde, r, roh_1));

		if(roh_1 == 0){
			KRATOS_GPU_CHECK(GPU_VectorNorm2(r, mResidualNorm));
			break;
		}

		if(i == 1){
			//p = r;
			KRATOS_GPU_CHECK(p.CopyFromGPU(r));
		}else{
			beta = (roh_1/roh_2) * (alpha/omega);
      			//p = r + beta * (p - omega * v);
			KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(-omega, v, 1.0, p));
			KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.0, r, beta, p));
		}
		if(havePreconditioner){
			//phat = M.solve(p);
			//GPU_fillWithZeros(size, phat.GPU_Values);
			KRATOS_GPU_CHECK(phat.CopyFromGPU(p));
			preconditioner.singleStep(p.GPU_Values, phat.GPU_Values);
		}else{
			KRATOS_GPU_CHECK(phat.CopyFromGPU(p));
		}
		
		//v = A * phat;
		KRATOS_GPU_CHECK(GPU_MatrixVectorMultiply(gA, phat, v));
		//dot(rtilde, v)
		double dotRtildeV = 0.0;
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(rtilde, v, dotRtildeV));

		alpha = roh_1 / dotRtildeV;

		//s = r - alpha(0) * v;
		KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.0, r, -alpha, v, s));

		//norm(s)
		KRATOS_GPU_CHECK(GPU_VectorNorm2(s, norms));

		if ((resid = norms/mBNorm) < tol) {
			//x += alpha(0) * phat;
			//tol = resid;
			//return 0;
			break;
		}
		if(havePreconditioner){
			//shat = M.solve(s);
			//GPU_fillWithZeros(size, shat.GPU_Values);
			KRATOS_GPU_CHECK(shat.CopyFromGPU(s));
			preconditioner.singleStep(s.GPU_Values, shat.GPU_Values);
		}else{
			KRATOS_GPU_CHECK(shat.CopyFromGPU(s));
		}

		//t = A * shat;
		KRATOS_GPU_CHECK(GPU_MatrixVectorMultiply(gA, shat, t));

		//omega = dot(t,s) / dot(t,t);
		double dotTS = 0.0, dotTT = 0.0;
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(t, s, dotTS));
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(t, t, dotTT));
		omega = dotTS / dotTT;
		
		//x += alpha(0) * phat + omega(0) * shat;
		KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd_addingVersion(alpha, phat, omega, shat, gX));
		//r = s - omega(0) * t;
		KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.0, s, -omega, t, r));

		roh_2 = roh_1;
		//norm(r)
		KRATOS_GPU_CHECK(GPU_VectorNorm2(r, mResidualNorm));
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
//		std::cout << mResidualNorm << std::endl;
//		std::cout << mBNorm << std::endl;
//		std::cout << mResidualNorm/mBNorm << std::endl;
	}while((mIterationsNumber < maxIterations) && (mResidualNorm > tol * mBNorm));
	KRATOS_GPU_CHECK(gX.Copy(GPU_CPU));
	if(havePreconditioner)
		preconditioner.cleanPreconditioner();
}

void CG_GPU(size_t A_Size1, size_t A_size2, size_t A_NNZ, double *A_values, size_t *A_indices, size_t *A_ptr, size_t vectorSize, double *X_values, double *B_values, 
	double tol, size_t maxIterations, double &mBNorm, double &mResidualNorm, size_t &mIterationsNumber, GPUPreconditioner& preconditioner)
      {

	bool havePreconditioner = (&preconditioner != 0);
	const int size = vectorSize;

	//Allocating matrix A
	GPUCSRMatrix gpuA(A_NNZ, A_Size1, A_size2, A_indices, A_ptr, A_values, true);
	KRATOS_GPU_CHECK(gpuA.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuA.Copy(CPU_GPU, false));

	//Allocating vector b
	GPUVector gpuB(vectorSize, B_values);
	KRATOS_GPU_CHECK(gpuB.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuB.Copy(CPU_GPU));

	//Allocating vector x
	GPUVector gpuX(vectorSize, X_values);
	KRATOS_GPU_CHECK(gpuX.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuX.Copy(CPU_GPU));

	double alpha, beta, roh, roh_1 = 1.0;
		
	//Ini preconditioner
	//clock_t s1 = clock();
	if(havePreconditioner){
		preconditioner.initialize(gpuA.CPU_RowIndices, gpuA.CPU_Columns, gpuA.CPU_Values,
				gpuA.GPU_RowIndices, gpuA.GPU_Columns, gpuA.GPU_Values,
				gpuA.Size1, gpuA.Size2, gpuA.NNZ, true, true);
	}
	//clock_t s2 = clock();
	//std::cout << "Time to create hierarchy" << double(s2-s1) / CLOCKS_PER_SEC << "s" << std::endl;

	//Norm(b)
	KRATOS_GPU_CHECK(GPU_VectorNorm2(gpuB, mBNorm));

	//r = b - A*x
	GPUVector gpuR(size);
	KRATOS_GPU_CHECK(gpuR.GPU_Allocate());
	KRATOS_GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuX, gpuR));
	KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuB, -1.00, gpuR));

	if(mBNorm == 0.0)
		mBNorm = 1.0;

	mIterationsNumber = 0;
	//Norm(r)
	double resid;
	KRATOS_GPU_CHECK(GPU_VectorNorm2(gpuR, mResidualNorm));
	if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
		if(havePreconditioner)
			preconditioner.cleanPreconditioner();
		std::cout << "leaving 1" << std::endl;
		return;
	}
	GPUVector gpuP(size), gpuZ(size), gpuQ(size);
	KRATOS_GPU_CHECK(gpuP.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuZ.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuQ.GPU_Allocate());
	
	//clock_t s3 = 0.0;
	size_t i = 1;
	do{
		//s1 = clock();
		if(havePreconditioner){
			GPU_fillWithZeros(size, gpuZ.GPU_Values);
			//gpuZ.CopyFromGPU(gpuR);		
			preconditioner.singleStep(gpuR.GPU_Values, gpuZ.GPU_Values);	
		}else{
			gpuZ.CopyFromGPU(gpuR);
		}
		//s2 = clock();
		//s3 += s2-s1;
		//roh = GPU_dotProduct(gpuR.Size, gpuR.GPU_Values, 1, gpuZ.GPU_Values, 1);
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(gpuR, gpuZ, roh));

		if(i == 1)
			//gpuP = gpuZ
			gpuP.CopyFromGPU(gpuZ);
		else{
			beta = roh / roh_1; 
			KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuZ, beta, gpuP));
		}

		KRATOS_GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuP, gpuQ));
		//roh_1 = dotProduct(p, q);
		double dotPQ;
		KRATOS_GPU_CHECK(GPU_VectorVectorMultiply(gpuP, gpuQ, dotPQ));
		alpha = roh / dotPQ;
		KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(alpha, gpuP, 1.00, gpuX));
		KRATOS_GPU_CHECK(GPU_VectorScaleAndAdd(-alpha, gpuQ, 1.00, gpuR));
	
		KRATOS_GPU_CHECK(GPU_VectorNorm2(gpuR, mResidualNorm));

		if((resid = mResidualNorm/mBNorm) <= 1.0e-30){
			break;
					std::cout << "leaving 2" << std::endl;
		}
		roh_1 = roh;
		mIterationsNumber++;			
		i++;
	}while((mIterationsNumber < maxIterations) && (mResidualNorm > tol * mBNorm));
	//std::cout << "Average time for single step" << (double(s3)/(i-1)) / CLOCKS_PER_SEC << "s" << std::endl;
			std::cout << "leaving 3" << std::endl;
	KRATOS_GPU_CHECK(gpuX.Copy(GPU_CPU));
	if(havePreconditioner)
		preconditioner.cleanPreconditioner();
}

//to add in gpu_sparse.cu

#endif
