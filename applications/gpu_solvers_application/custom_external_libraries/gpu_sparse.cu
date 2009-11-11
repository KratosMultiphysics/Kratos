/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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

//
// Sparse matrix and vector operations on GPU

//
// Notes!
//
//   * In case of VectorVectorMultiply and VectorNorm2, cuBlas library has been used, as it is much faster than a code without proper optimizations.
//
//   * In case of VectorScaleAndAdd, as there is no direct way of using cuBlas for this, we had to write our own kernel. It was 2 or 3 times faster
//     than using cuBlas indirectly.
//
//   * For compilation command, simply issue a 'make' command.
//
//   * Removed templates to be able to link to Kratos.
//
//   * Modified GPUCSRMatrix to keep no. of non-zeros per row constant and equal to HALF_WARP_SIZE (16) and used same no. of threads to multiply a row

// More notes!
//
//	* Checks for error in Bind_X() and Unbind_X() have been removed for more consistent error checking after kernel calls; can they fail?
//
//	* Added an optional parameter to GPUCSRMatrix constructor, so that user can optionally avoid making non-zeros in a row a multiple of HALF_WARP_SIZE (16)
//
//	* Added a make file; in the command line use emu=1 for emulation mode and dbg=1 for a debug version

// Includes, system

#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <cstdio>

// Includes, project

#include "gpu_sparse.h"
#include "gpu_sparse_utils.h"
#include "gpu_sparse_kernels.h"

//Includes, preconditioner
#include <lapackd.h>
#include <laslv.h>
#include <gmd.h>
#include <gmf.h>
#include <vector>

namespace Kratos
{

namespace GPUSparse
{

// GPUVector class definition

GPUVector::GPUVector(size_t _Size, double *_CPU_Values): Size(_Size), CPU_Values(_CPU_Values), GPU_Values(0), Allocated(false)
{
	// Nothing to do!
}

GPUVector::GPUVector(size_t _Size): Size(_Size), CPU_Values(0), GPU_Values(0), Allocated(false)
{
	// Nothing to do!
}

GPUVector::~GPUVector()
{
	if (Allocated)
		GPU_Free();
}

bool GPUVector::GPU_Allocate()
{
	if (Allocated)
		return false;

	Allocated = true;

	return CUDA_Success(cudaMalloc(reinterpret_cast <void **> (&GPU_Values), Size * sizeof(double)));
}

bool GPUVector::GPU_Free()
{
	if (!Allocated)
		return false;

	Allocated = false;

	return CUDA_Success(cudaFree(GPU_Values));
}

bool GPUVector::Copy(CopyDirection Direction)
{
	if (!Allocated)
		return false;

	switch (Direction)
	{
		case CPU_GPU:

			return CUDA_Success(cudaMemcpy(GPU_Values, CPU_Values, Size * sizeof(double), cudaMemcpyHostToDevice));

		case GPU_CPU:

			return CUDA_Success(cudaMemcpy(CPU_Values, GPU_Values, Size * sizeof(double), cudaMemcpyDeviceToHost));
	}

	// We should never get here!
	return false;
}

bool GPUVector::CopyFromGPU(GPUVector &V)
{
	if (V.Size != Size || !Allocated || !V.Allocated)
		return false;

	return CUDA_Success(cudaMemcpy(GPU_Values, V.GPU_Values, Size * sizeof(double), cudaMemcpyDeviceToDevice));
}

// GPUCSRMatrix class definition

GPUCSRMatrix::GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2, size_t *_CPU_Columns, size_t *_CPU_RowIndices, double *_CPU_Values, bool _NZMultiple16): NNZ(_NNZ), Size1(_Size1), Size2(_Size2), CPU_Columns(0), CPU_RowIndices(0), CPU_Values(0), GPU_Columns(0), GPU_RowIndices(0), GPU_Values(0), Allocated(false)

{

	if (_NZMultiple16)
	{

		NNZ = 0;
		// Temporary RowIndices vector
		size_t *Temp_CPU_RowIndices = new size_t[Size1 + 1];

		Temp_CPU_RowIndices[0] = 0;

		// Find out how many non-zeros are needed to pad all rows to 16 while building the RowIndices
		for (size_t i = 0; i < Size1; i++)
		{
			size_t NZ = _CPU_RowIndices[i + 1] - _CPU_RowIndices[i];

			size_t R = NZ & HALF_WARP_SIZE_MASK;

			if (R != 0)
				NZ += HALF_WARP_SIZE - R;

			NNZ += NZ;
			Temp_CPU_RowIndices[i + 1] = Temp_CPU_RowIndices[i] + NZ;
		}

		// Allocate CPU memory for CSR structure using only one chunk of page-locked memory to speed up data transfer between CPU and GPU
		void *CSR_Data;


		if (!CUDA_Success(cudaMallocHost(&CSR_Data, NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t))))  // TODO: What should be done?!
			CSR_Data = 0;

		// We are sure that using this order, the memory alignment conditions will be satisfied as NNZ is a multiple of HALF_WARP_SIZE (16) and sizeof(double) = 8
		// TODO: Check this!
/*		CPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		CPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		CPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)*/
CPU_Values = (double*)CSR_Data;
CPU_Columns = (size_t*)((double*) CSR_Data + NNZ);
CPU_RowIndices = (size_t*)((size_t*)((double*) CSR_Data + NNZ) + NNZ);
		// Move temporary data

		memcpy(CPU_RowIndices, Temp_CPU_RowIndices, (Size1 + 1) * sizeof(size_t));


		delete[] Temp_CPU_RowIndices;

		// Build ECSR structure from given CSR
		for (size_t i = 0; i < Size1; i++)
		{
			size_t _Start = _CPU_RowIndices[i], Start = CPU_RowIndices[i];

			for (size_t j = 0; j < _CPU_RowIndices[i + 1] - _CPU_RowIndices[i]; j++)
			{
				CPU_Columns[Start + j] = _CPU_Columns[_Start + j];
				CPU_Values[Start + j] = _CPU_Values[_Start + j];
			}

			size_t LastCol = _CPU_Columns[_CPU_RowIndices[i + 1] - 1];

			for (size_t j = _CPU_RowIndices[i + 1] - _CPU_RowIndices[i]; j < CPU_RowIndices[i + 1] - CPU_RowIndices[i]; j++)
			{
				CPU_Columns[Start + j] = LastCol;  // To maintain coalescing as much as possible
				CPU_Values[Start + j] = 0.00;
			}
		}
	}

	else
	{
		// Allocate CPU memory for CSR structure using only one chunk of page-locked memory to speed up data transfer between CPU and GPU
		void *CSR_Data;

		if (!CUDA_Success(cudaMallocHost(&CSR_Data, NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t))))  // TODO: What should be done?!
			CSR_Data = 0;

		// We are sure that using this order, the memory alignment conditions will be satisfied as NNZ is a multiple of HALF_WARP_SIZE (16) and sizeof(double) = 8
		// TODO: Check this!
/*		CPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		CPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		CPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)*/
CPU_Values = (double*)CSR_Data;
CPU_Columns = (size_t*)((double*) CSR_Data + NNZ);
CPU_RowIndices = (size_t*)((size_t*)((double*) CSR_Data + NNZ) + NNZ);

		// Move temporary data
		memcpy(CPU_Values, _CPU_Values, NNZ * sizeof(double));
		memcpy(CPU_Columns, _CPU_Columns, NNZ * sizeof(size_t));
		memcpy(CPU_RowIndices, _CPU_RowIndices, (Size1 + 1) * sizeof(size_t));
	}
}

GPUCSRMatrix::GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2): NNZ(_NNZ), Size1(_Size1), Size2(_Size2), CPU_Columns(0), CPU_RowIndices(0), CPU_Values(0), GPU_Columns(0), GPU_RowIndices(0), GPU_Values(0), Allocated(false)
{
	// Nothing to do!
}

GPUCSRMatrix::~GPUCSRMatrix()
{
	// Free CSR data; as it is allocated in one chunk of memory, we need only to free the begining address
	cudaFreeHost(CPU_Values);

	if (Allocated)
		GPU_Free();
}

bool GPUCSRMatrix::GPU_Allocate()
{
	if (Allocated)
		return false;

	Allocated = true;

	// Allocate GPU memory for CSR structure using only one chunk of memory to speed up data transfer between CPU and GPU

	void *CSR_Data;


	if (CUDA_Success(cudaMalloc(&CSR_Data, NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t))))
	{
/*		GPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		GPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		GPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)*/
GPU_Values = (double*)CSR_Data;
		GPU_Columns = (size_t*)((double*) CSR_Data + NNZ);
		GPU_RowIndices = (size_t*)((size_t*)((double*) CSR_Data + NNZ) + NNZ);

		return true;
	}

	else
		return false;
}

bool GPUCSRMatrix::GPU_Free()
{
	if (!Allocated)
		return false;

	Allocated = false;
	
	return CUDA_Success(cudaFree(GPU_Values));
}

bool GPUCSRMatrix::Copy(CopyDirection Direction, bool CopyValuesOnly)
{
	if (!Allocated)
		return false;

	switch (Direction)
	{
		case CPU_GPU:
			if (CopyValuesOnly)
				return CUDA_Success(cudaMemcpy(GPU_Values, CPU_Values, NNZ * sizeof(double), cudaMemcpyHostToDevice));
			else
				return CUDA_Success(cudaMemcpy(GPU_Values, CPU_Values, NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t), cudaMemcpyHostToDevice));

		case GPU_CPU: 
			if (CopyValuesOnly)
				return CUDA_Success(cudaMemcpy(CPU_Values, GPU_Values, NNZ * sizeof(double), cudaMemcpyDeviceToHost));
			else
				return CUDA_Success(cudaMemcpy(CPU_Values, GPU_Values, NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t), cudaMemcpyDeviceToHost));

	}

	// We should never get here!
	return false;
}

bool GPUCSRMatrix::CopyFromGPU(GPUCSRMatrix &M, bool CopyStructure, bool CopyValues)
{
	if (M.Size1 != Size1 || M.Size2 != Size2 || M.NNZ != NNZ || !Allocated || !M.Allocated)
		return false;

	size_t CopyLength;
	void *CopyFrom, *CopyTo;

	if (CopyStructure && CopyValues)
	{
		CopyTo = GPU_Values;
		CopyFrom = M.GPU_Values;
		CopyLength = NNZ * (sizeof(double) + sizeof(size_t)) + (Size1 + 1) * sizeof(size_t);
	}
	else if (CopyStructure && !CopyValues)
	{
		CopyTo = reinterpret_cast <void *> (reinterpret_cast <size_t> (GPU_Values) + NNZ * sizeof(double));
		CopyFrom = reinterpret_cast <void *> (reinterpret_cast <size_t> (M.GPU_Values) + NNZ * sizeof(double));
		CopyLength = (NNZ + Size1 + 1) * sizeof(size_t);
	}
	else if (!CopyStructure && CopyValues)
	{
		CopyTo = GPU_Values;
		CopyFrom = M.GPU_Values;
		CopyLength = NNZ * sizeof(double);
	}
	else if (!CopyStructure && !CopyValues)
	{
		CopyTo = 0;
		CopyFrom = 0;
		CopyLength = 0;
	}

	if (CopyLength != 0)
		return CUDA_Success(cudaMemcpy(CopyTo, CopyFrom, CopyLength, cudaMemcpyDeviceToDevice));
	else
		return true;
}

// Operations defined on GPUCSRMatrix and GPUVector

//
// CPU_MatrixVectorMultiply
// Matrix-Vector multiply on CPU

bool CPU_MatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y)
{
	// Primary checks
	if (A.Size2 != X.Size || X.Size != Y.Size)
		return false;

	for (size_t i = 0; i < A.Size1; i++)
	{
		double YI = static_cast <double> (0);

		for (size_t j = A.CPU_RowIndices[i]; j < A.CPU_RowIndices[i + 1]; j++)
			YI += A.CPU_Values[j] * X.CPU_Values[A.CPU_Columns[j]];

		Y.CPU_Values[i] = YI;
	}

	return true;
}

//
// GPU_MatrixVectorMultiply
// Matrix-Vector multiply on GPU

bool GPU_MatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y)
{
	// Primary checks
	if (A.Size2 != X.Size || X.Size != Y.Size || !X.Allocated || !Y.Allocated)
		return false;

#ifdef USE_TEXTURE_CACHING

	// Bind the texture memory to X
	Bind_X(X.GPU_Values);

#endif


	bool UseVectorizedVersion = (A.NNZ / A.Size2) > 10;	// From nVidia forum

	if (UseVectorizedVersion)
	{
		dim3 Grid = Build_Grid(A.Size1 *  HALF_WARP_SIZE, BLOCK_SIZE);
		GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values, Y.GPU_Values);

		if (!GPUSparse::CUDA_Success(cudaGetLastError()))
			return false;
	}

	else

	{
		dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
		GPU_MatrixVectorMultiply_CSR_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values, Y.GPU_Values);

		if (!GPUSparse::CUDA_Success(cudaGetLastError()))
			return false;
	}

#ifdef USE_TEXTURE_CACHING

	// Unbind the texture memory
	Unbind_X();
	}

#endif

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPU_MatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on CPU

bool CPU_MatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X)
{
	// Primary checks
	if (A.Size1 != A.Size2 || A.Size2 != X.Size)
		return false;

	for (size_t i = 0; i < A.Size1; i++)
	{
		X.CPU_Values[i] = static_cast <double> (0);

		for (size_t j = A.CPU_RowIndices[i]; j < A.CPU_RowIndices[i + 1]; j++)
			if (A.CPU_Columns[j] == i)
				X.CPU_Values[i] = A.CPU_Values[j];
	}

	return true;
}

//
// GPU_MatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on GPU

bool GPU_MatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X)
{
	// Primary checks
	if (A.Size1 != A.Size2 || A.Size2 != X.Size || !A.Allocated || !X.Allocated)
		return false;

	dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPU_MatrixGetDiagonals_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPU_MatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on CPU

bool CPU_MatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A)
{
	// Primary checks
	if (X.Size != A.Size1)
		return false;

	for (size_t i = 0; i < X.Size; i++)
	{
		double t = X.CPU_Values[i];

		for (size_t j = A.CPU_RowIndices[i]; j < A.CPU_RowIndices[i + 1]; j++)
				A.CPU_Values[j] *= t;
	}

	return true;
}

//
// GPU_MatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on GPU

bool GPU_MatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A)
{
	// Primary checks
	if (X.Size != A.Size1 || !X.Allocated || !A.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);
	GPU_MatrixMatrixDiagonalMultiply_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPU_VectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on CPU

bool CPU_VectorPrepareDiagonalPreconditionerValues(GPUVector &X)
{
	for (size_t i = 0; i < X.Size; i++)
		if (X.CPU_Values[i] == 0.00)
			X.CPU_Values[i] = 1.00;
		else
			X.CPU_Values[i] = 1.00 / X.CPU_Values[i];
//			X.CPU_Values[i] = 1.00 / sqrt(abs(X.CPU_Values[i]));

	return true;
}

//
// GPU_VectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on GPU

bool GPU_VectorPrepareDiagonalPreconditionerValues(GPUVector &X)
{
	// Primary check
	if (!X.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);
	GPU_VectorPrepareDiagonalPreconditionerValues_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// GPU_PrepareSPAIPreconditioner
// Prepare SPAI preconditioner on GPU

bool GPU_PrepareSPAIPreconditioner(GPUCSRMatrix &A, GPUCSRMatrix &M)
{
	// Primary checks
	if (A.Size1 != M.Size1 || A.Size2 != M.Size2 || A.NNZ != M.NNZ || !A.Allocated || !M.Allocated)
		return false;

	dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPU_SPAIPreconditioner_CSR_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, M.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPU_VectorVectorMultiply
// Vector-Vector multiply on CPU

bool CPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result)
{
	// Primary check
	if (X.Size != Y.Size)
		return false;

	Result = static_cast <double> (0);

	for (size_t i = 0; i < X.Size; i++)
		Result += X.CPU_Values[i] * Y.CPU_Values[i];

	return true;
}

//
// GPU_VectorVectorMultiply
// Vector-Vector multiply on GPU

bool GPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result)
{
	// Primary check
	if (X.Size != Y.Size || !X.Allocated || !Y.Allocated)
		return false;

	Result = cublasDdot(X.Size, X.GPU_Values, 1, Y.GPU_Values, 1);

	return CUBLAS_Success(cublasGetError());
}

//
// CPU_VectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on CPU

bool CPU_VectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y,  GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Z.CPU_Values[i] = X.CPU_Values[i] * Y.CPU_Values[i];

	return true;
}

//
// GPU_VectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on GPU

bool GPU_VectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size || !X.Allocated || !Y.Allocated || !Z.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	GPU_VectorVectorMultiplyElementWise_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values, Y.GPU_Values, Z.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPU_VectorNorm2
// Vector norm 2 on CPU

bool CPU_VectorNorm2(GPUVector &X, double &Result)
{
	Result = static_cast <double> (0);

	for (size_t i = 0; i < X.Size; i++)
		Result += X.CPU_Values[i] * X.CPU_Values[i];

	Result = sqrt(Result);

	return true;
}

//
// GPU_VectorNorm2
// Vector norm 2 on GPU

bool GPU_VectorNorm2(GPUVector &X, double &Result)
{
	// Primary check
	if (!X.Allocated)
		return false;

	Result = cublasDnrm2(X.Size, X.GPU_Values, 1);

	return CUBLAS_Success(cublasGetError());
}

//
// CPU_VectorScaleAndAdd
// Vector scale-and-add on CPU

// Variant 1: Z = A * X + B * Y

bool CPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Z.CPU_Values[i] = A * X.CPU_Values[i] + B * Y.CPU_Values[i];

	return true;
}

// Variant 2: Y = A * X + B * Y

bool CPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y)
{
	// Primary check
	if (X.Size != Y.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Y.CPU_Values[i] = A * X.CPU_Values[i] + B * Y.CPU_Values[i];

	return true;
}

//
// GPU_VectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool GPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size || !X.Allocated || !Y.Allocated || !Z.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	if (A == 1.00)
	{
	//printf("A = 1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_A_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_B_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPU_VectorScaleAndAdd_1_E_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else if (A == -1.00)
	{
	//printf("A = -1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_C_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_D_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPU_VectorScaleAndAdd_1_F_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else
	{
	//printf("B = 1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_G_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_H_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPU_VectorScaleAndAdd_1_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

// Variant 2: Y = A * X + B * Y

bool GPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y)
{
	// Primary check
	if (X.Size != Y.Size || !X.Allocated || !Y.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	if (A == 1.00)
	{
		if (B == 1.00)
			GPU_VectorScaleAndAdd_2_A_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_2_B_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPU_VectorScaleAndAdd_2_E_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
	}

	else if (A == -1.00)
	{
		if (B == 1.00)
			GPU_VectorScaleAndAdd_2_C_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_2_D_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPU_VectorScaleAndAdd_2_F_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
	}

	else
	{
		if (B == 1.00)
			GPU_VectorScaleAndAdd_2_G_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_2_H_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPU_VectorScaleAndAdd_2_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
	}

	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

/*double GPU_dotProduct(size_t numElems, const double *firstVec, int incFirstVec,
    const double *secondVec, int incSecondVec){
    return cublasDdot(numElems, firstVec, incFirstVec, secondVec, incSecondVec);
}*/


/** ADDED FUNCTIONS **/

void GPU_fillWithZeros(size_t numElems, double* gpuVec){
	dim3 grid = Build_Grid(numElems, BLOCK_SIZE);
	fillWithZeros <<< grid, BLOCK_SIZE >>> (gpuVec, numElems);
}

bool GPU_VectorScaleAndAdd_addingVersion(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size || !X.Allocated || !Y.Allocated || !Z.Allocated){
		//printf("Falla x la comprovació\n");		
		return false;
	}
		
	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	if (A == 1.00)
	{
	//printf("A = 1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_A_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

			
		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_B_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
				
		else
			GPU_VectorScaleAndAdd_1_E_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	
	else if (A == -1.00)
	{
	//printf("A = -1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_C_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_D_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else
			GPU_VectorScaleAndAdd_1_F_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	
	else
	{
	//printf("B = 1.0\n");
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_G_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_H_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else
			GPU_VectorScaleAndAdd_1_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
	
}

/** temp variables for LU decomposition **/
LaVectorLongInt ipiv;

/** functions from scipy for mat-mat calculation **/
template <class I>
void csr_matmat_pass1(const I n_row,
                      const I n_col,
                      const I Ap[],
                      const I Aj[],
                      const I Bp[],
                      const I Bj[],
                            I Cp[]){
    //std::vector<I> mask(n_col,-1);
    int* mask = new int[n_col];
    for(size_t aux = 0; aux < n_col; aux++){
	mask[aux] = -1;
    }
    Cp[0] = 0;

    I nnz = 0;
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];
            for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
                I k = Bj[kk];
                if(mask[k] != (int)i){
                    mask[k] = (int)i;
                    nnz++;
                }
            }
        }
        Cp[i+1] = nnz;
    }
    delete[] mask;
}

template <class I, class T>
void csr_matmat_pass2(const I n_row,
      	              const I n_col,
      	              const I Ap[],
      	              const I Aj[],
      	              const T Ax[],
      	              const I Bp[],
      	              const I Bj[],
      	              const T Bx[],
      	                    I Cp[],
      	                    I Cj[],
      	                    T Cx[])
{
//    std::vector<I> next(n_col,-1);
 //   std::vector<T> sums(n_col, 0);

    int* next = new int[n_col];
    T* sums = new T[n_col];
    for(size_t aux = 0; aux < n_col; aux++){
	next[aux] = -1;
	sums[aux] = 0;
    }

    I nnz = 0;

    Cp[0] = 0;

    for(I i = 0; i < n_row; i++){
        int head   = -2;
        I length =  0;

        I jj_start = Ap[i];
        I jj_end   = Ap[i+1];
        for(I jj = jj_start; jj < jj_end; jj++){
            I j = Aj[jj];
            T v = Ax[jj];

            I kk_start = Bp[j];
            I kk_end   = Bp[j+1];
            for(I kk = kk_start; kk < kk_end; kk++){
                I k = Bj[kk];

                sums[k] += v*Bx[kk];

                if(next[k] == -1){
                    next[k] = head;
                    head = (int)k;
                    length++;
                }
            }
        }

        for(I jj = 0; jj < length; jj++){

            if(sums[head] != 0){
                Cj[nnz] = (size_t)head;
                Cx[nnz] = sums[head];
                nnz++;
            }

            I temp = (size_t)head;
            head = (I)next[head];

            next[temp] = -1; //clear arrays
            sums[temp] =  0;
        }

        Cp[i+1] = nnz;
    }

	delete[] sums;
	delete[] next;
}

/** maxLevels define the maxLevels of that execution
	G defines the diagonals of each lvl of A, created on previous step**/
void multilevel(_Matrix*& A, _Matrix*& P, _Matrix*& R, _Matrix*& G, _Vector& b, _Vector& u,
			unsigned short lvl, unsigned short maxLevels, size_t* preSweeps, size_t* postSweeps, bool assumeZeros)
{
    bool vectorized = (A[lvl].numNNZ / A[lvl].numCols) > 10;
    _Vector r;
    //calculateInstantVector(u, b, A[lvl], G[lvl]);
    if(lvl < maxLevels){
	//clock_t t1 = clock();
        if(assumeZeros) //we receive from the upper level a zero start vector
        {
            if(preSweeps[lvl] != 0)
            {
                //first iteration (does not require computation of residual
                dim3 Grid = Build_Grid(A[lvl].numRows, BLOCK_SIZE);
                GPU_MatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BLOCK_SIZE >>>(G[lvl].numRows, G[lvl].indices_gpu,
		    G[lvl].ptr_gpu, G[lvl].values_gpu, b.values_gpu, u.values_gpu);
                if(!CUDA_Success(cudaThreadSynchronize())){
                    cout << "Error en linea 130" << endl;
                }
                //from the second sweel on we need to recompute the residual
                for(size_t i = 1; i < preSweeps[lvl]; i++){
		    if(!vectorized)
			calculateInstantVector(u, b, A[lvl], G[lvl]);
		    else
			calculateInstantVector_vectorized(u, b, A[lvl], G[lvl]);
		}

                if(!vectorized)
                    generateResidual(R[lvl], b, A[lvl], u, r);
                else
                    generateResidual_vectorized(R[lvl], b, A[lvl], u, r);


            }
            else //preSweeps[0] == 0 case
            {
                //inefficient! -- in this case we do not need to recompute the residual
                if(!vectorized)
                    generateResidual(R[lvl], b, A[lvl], u, r);
                else
                    generateResidual_vectorized(R[lvl], b, A[lvl], u, r);

            }
        }
        else
        {
          //from the second sweel on we need to recompute the residual
                for(size_t i = 0; i < preSweeps[lvl]; i++){
		    if(!vectorized)
			calculateInstantVector(u, b, A[lvl], G[lvl]);
		    else
			calculateInstantVector_vectorized(u, b, A[lvl], G[lvl]);
		}


                if(!vectorized)
                    generateResidual(R[lvl], b, A[lvl], u, r);
                else
                    generateResidual_vectorized(R[lvl], b, A[lvl], u, r);
        }
        _Vector v;
        v.numElems = r.numElems;
        v.values_cpu = new double[v.numElems];
        malloc_(v.values_gpu, v.numElems);
        dim3 Grid = Build_Grid(v.numElems, BLOCK_SIZE);
        fillWithZeros <<< Grid, BLOCK_SIZE >>>(v.values_gpu, v.numElems);
        if(!CUDA_Success(cudaThreadSynchronize())){
                cout << "Error en linea 160" << endl;
        }

        multilevel(A, P, R, G, r, v, lvl+1, maxLevels, preSweeps, postSweeps, assumeZeros);

        _Vector pv;
        pv.numElems = P[lvl].numRows;
        // malloc de pv.values_gpu
        malloc_(pv.values_gpu, pv.numElems);
        Grid = Build_Grid(pv.numElems, BLOCK_SIZE);
        // here product matrix P with vector v
        GPU_MatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(P[lvl].numRows,
                P[lvl].indices_gpu, P[lvl].ptr_gpu, P[lvl].values_gpu, v.values_gpu, pv.values_gpu);
        if(!CUDA_Success(cudaThreadSynchronize())){
            cout << "Error en linea 173" << endl;
        }
        // here the addition of pv to u
        Grid = Build_Grid(u.numElems, BLOCK_SIZE);
        sumVectorVector <<< Grid, BLOCK_SIZE >>> (pv.values_gpu, u.values_gpu, u.numElems);
        if(!CUDA_Success(cudaThreadSynchronize())){
            cout << "Error en linea 180" << endl;
        }
        //delete de pv, v i r
        CUDA_CHECK(cudaFree(pv.values_gpu));
        CUDA_CHECK(cudaFree(v.values_gpu));
        CUDA_CHECK(cudaFree(r.values_gpu));
        delete[] r.values_cpu;
        delete[] v.values_cpu;

        //double norm2 = checkResidual(u, b, A[lvl]);

	for(size_t i = 0; i < postSweeps[lvl]; i++){
	    if(!vectorized)
		calculateInstantVector(u, b, A[lvl], G[lvl]);
	    else
		calculateInstantVector_vectorized(u, b, A[lvl], G[lvl]);
	}

    }else{
	//clock_t t1 = clock();
        //here lapack direct solver

        copyMem(u.values_gpu, u.values_cpu, u.numElems, 1);
        copyMem(b.values_gpu, b.values_cpu, b.numElems, 1);


        LaGenMatDouble _A(A[lvl].matAuxValues, A[lvl].numRows, A[lvl].numCols);
        LaGenMatDouble _b(b.values_cpu, b.numElems, 1);
        LaGenMatDouble _x(u.values_cpu, u.numElems, 1);

    	_x.inject(_b);            // will throw exception if not conformant

	integer info = 0;
	int M = _A.size(0);
	integer Ml = M;
	integer lda = _A.inc(0) * _A.gdim(0);

	integer K = _x.size(1);
	integer ldx = _x.inc(0) * _x.gdim(0);
	F77NAME(dgetrs) ("No transpose", &Ml, &K, &_A(0,0), &lda, &ipiv(0), &_x(0,0), &ldx, &info);

	//int res = clapack_dgetrs(CblasRowMajor, CblasNoTrans, &Ml, &K, &_A(0,0), &lda, &ipiv(0), &_x(0,0), &ldx);
	//std::cout << "Problem with lapack, num " << res << std:endl;
	copyMem(u.values_cpu, u.values_gpu, u.numElems, 0);

	//clock_t t2 = clock();
	//cout << "Lower lvl timing " << double(t2-t1) / CLOCKS_PER_SEC << "s" << endl;
   }

}
/** This function is a wrapper for u += G ( A, b, u) **/
void calculateInstantVector(_Vector& u, const _Vector& b, const _Matrix& A, const _Matrix& G)
{
    /** Au **/
    _Vector auxAU;
    auxAU.numElems = A.numRows;
    malloc_(auxAU.values_gpu, auxAU.numElems);
    dim3 Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    GPU_MatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.numRows, A.indices_gpu,
            A.ptr_gpu, A.values_gpu, u.values_gpu, auxAU.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 238" << endl;
    }

    /** b - AU **/
    _Vector auxABU;
    auxABU.numElems = auxAU.numElems;
    malloc_(auxABU.values_gpu, auxABU.numElems);
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.values_gpu, auxAU.values_gpu, auxABU.values_gpu, auxABU.numElems);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 249" << endl;
    }
    /** u += G ( b - Au ) **/
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    GPU_MatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BLOCK_SIZE >>>(G.numRows, G.indices_gpu,
            G.ptr_gpu, G.values_gpu, auxABU.values_gpu, u.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 255" << endl;
    }
    //deleting structures
    CUDA_CHECK(cudaFree(auxABU.values_gpu));
    CUDA_CHECK(cudaFree(auxAU.values_gpu));

}
void calculateInstantVector_vectorized(_Vector& u, const _Vector& b, const _Matrix& A, const _Matrix& G)
{
    /** Au **/
    _Vector auxAU;
    auxAU.numElems = A.numRows;
    malloc_(auxAU.values_gpu, auxAU.numElems);
    dim3 Grid = Build_Grid(A.numRows * HWS, BS);
    GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel <<< Grid, BS >>>(A.numRows, A.indices_gpu,
            A.ptr_gpu, A.values_gpu, u.values_gpu, auxAU.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 272" << endl;
    }
    /** b - AU **/
    _Vector auxABU;
    auxABU.numElems = auxAU.numElems;
    malloc_(auxABU.values_gpu, auxABU.numElems);
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.values_gpu, auxAU.values_gpu, auxABU.values_gpu, auxABU.numElems);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 282" << endl;
    }
    /** u += G ( b - Au ) **/
    Grid = Build_Grid(A.numRows, BS);
    GPU_MatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BS >>>(G.numRows, G.indices_gpu,
            G.ptr_gpu, G.values_gpu, auxABU.values_gpu, u.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 288" << endl;
    }
    //deleting structures
    CUDA_CHECK(cudaFree(auxABU.values_gpu));
    CUDA_CHECK(cudaFree(auxAU.values_gpu));

}

void generateResidual(const _Matrix& R, const _Vector& b, const _Matrix& A, const _Vector& u, _Vector& r){
    /** Au **/
    _Vector auxAU;
    auxAU.numElems = A.numRows;
    malloc_(auxAU.values_gpu, auxAU.numElems);
    dim3 Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    GPU_MatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.numRows, A.indices_gpu,
            A.ptr_gpu, A.values_gpu, u.values_gpu, auxAU.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 305" << endl;
    }
    /** b - AU **/
    _Vector auxABU;
    auxABU.numElems = auxAU.numElems;
    malloc_(auxABU.values_gpu, auxABU.numElems);
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.values_gpu, auxAU.values_gpu,
            auxABU.values_gpu, auxABU.numElems);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 315" << endl;
    }
    /** r = R ( b - Au ) **/
    r.numElems = R.numRows;
    malloc_(r.values_gpu, r.numElems);
    r.values_cpu = new double[r.numElems];
    Grid = Build_Grid(R.numRows, BLOCK_SIZE);
    GPU_MatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(R.numRows, R.indices_gpu,
            R.ptr_gpu, R.values_gpu, auxABU.values_gpu, r.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 325" << endl;
    }
    //deleting structures
    CUDA_CHECK(cudaFree(auxABU.values_gpu));
    CUDA_CHECK(cudaFree(auxAU.values_gpu));
}
void generateResidual_vectorized(const _Matrix& R, const _Vector& b, const _Matrix& A, const _Vector& u, _Vector& r){
    /** Au **/
    _Vector auxAU;
    auxAU.numElems = A.numRows;
    malloc_(auxAU.values_gpu, auxAU.numElems);
    dim3 Grid = Build_Grid(A.numRows * HWS, BS);
    GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel <<< Grid, BS >>>(A.numRows, A.indices_gpu,
            A.ptr_gpu, A.values_gpu, u.values_gpu, auxAU.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 340" << endl;
    }
    /** b - AU **/
    _Vector auxABU;
    auxABU.numElems = auxAU.numElems;
    malloc_(auxABU.values_gpu, auxABU.numElems);
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.values_gpu, auxAU.values_gpu,
            auxABU.values_gpu, auxABU.numElems);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 350" << endl;
    }
    /** r = R ( b - Au ) **/
    r.numElems = R.numRows;
    malloc_(r.values_gpu, r.numElems);
    r.values_cpu = new double[r.numElems];
    Grid = Build_Grid(R.numRows *  HWS, BS);
    GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel <<< Grid, BS >>>(R.numRows, R.indices_gpu,
            R.ptr_gpu, R.values_gpu, auxABU.values_gpu, r.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 360" << endl;
    }
    //deleting structures
    CUDA_CHECK(cudaFree(auxABU.values_gpu));
    CUDA_CHECK(cudaFree(auxAU.values_gpu));
}

/** This function will return the hierarchy reconstructed of A and b;
    Additionally return the number of real hierarchy levels.
    All matrices and vectors are returned as GPU structures**/
size_t generateHierarchy(_Matrix*& Matrices, _Matrix*& Pmat, _Matrix*& Qmat,
        _Matrix*& Gmat, double W, size_t numLevelsRoh, size_t max_levels, size_t min_system_size)
{
    size_t i = 0;
    for(i = 0; i < max_levels; i++){
        /** This condition controls MAX_SYSTEM_SIZE for the last matrix in hierarchy **/
        if(Matrices[i].numRows < min_system_size || i == max_levels-1){
/*		cout << "BEFORE" << endl;
		for(size_t j = 0; j < Matrices[i].numNNZ; j++){
			cout << Matrices[i].values_cpu[j] << " ";
		}
		cout << endl << endl;*/
            computeDenseMatrix(Matrices[i], Matrices[i].matAuxValues);
            LaGenMatDouble A(Matrices[i].matAuxValues, Matrices[i].numRows, Matrices[i].numCols);

            int M = A.size(0);
            integer Ml = M;
            integer lda = A.inc(0) * A.gdim(0);

            LaVectorLongInt ipiv_( M);
            ipiv = ipiv_;
            integer info = 0;
            int max = M > 1 ? M : 1;
            if(M < 0){
                info = -1;
            }else if(M < 0){
                //este caso lo saltamos porque tratamos matrices cuadradas y por tanto se asume en el caso anterior
            }else if(lda < max){
                info  = -4;
            }
            F77NAME(dgetrf) (&Ml, &Ml, &A(0,0), &lda, &ipiv(0), &info);

		
	    if( i > 0 ){ //A matrix is not the first one in hierarchy
		    if(Matrices[i].ptr_cpu != NULL) delete[] Matrices[i].indices_cpu;
		    if(Matrices[i].values_cpu != NULL) delete[] Matrices[i].values_cpu;
		    if(Matrices[i].indices_cpu != NULL) delete[] Matrices[i].ptr_cpu;
	    }
            Matrices[i].numValuesDenseRep = A.inc(0) * A.gdim(0);
/*		cout << "AFTER" << endl;
		for(size_t j = 0; j < Matrices[i].numValuesDenseRep; j++){
			cout << Matrices[i].matAuxValues[j] << " ";
		}
		cout << endl << endl;
		exit(1);*/
            break;
        }
        _Matrix newDiag;
        _Vector diag;
        clock_t t1 = clock();
        createDiagonal_vCPU(Matrices[i], newDiag, diag);
        clock_t t2 = clock();
        Gmat[i] = newDiag;
        /** Generating P and Q for the current A level **/
        _Matrix P = generateP_vCPU(Matrices[i], diag, W, numLevelsRoh);
        clock_t t3 = clock();
        _Matrix Q = generateQ(P);
        clock_t t4 = clock();
        Pmat[i] = P;
        Qmat[i] = Q;

        //CUDA_CHECK(cudaFree(diag.values_gpu));
	delete[] diag.values_cpu;

        /** Allocating result matrix, and partialResult matrix **/
        _Matrix matResult, matPartialResult;

        matPartialResult.numRows = Q.numRows;
        matPartialResult.numCols = Matrices[i].numCols;
        matPartialResult.ptr_cpu = new size_t[matPartialResult.numRows+1];

        csr_matmat_pass1(Q.numRows,
                      Matrices[i].numCols,
                      Q.ptr_cpu,
                      Q.indices_cpu,
                      Matrices[i].ptr_cpu,
                      Matrices[i].indices_cpu,
                            matPartialResult.ptr_cpu);
        matPartialResult.numNNZ = matPartialResult.ptr_cpu[matPartialResult.numRows];
        matPartialResult.indices_cpu = new size_t[matPartialResult.numNNZ];
        matPartialResult.values_cpu = new double[matPartialResult.numNNZ];
        csr_matmat_pass2(Q.numRows,
      	              Matrices[i].numCols,
      	              Q.ptr_cpu,
                      Q.indices_cpu,
      	              Q.values_cpu,
      	              Matrices[i].ptr_cpu,
                      Matrices[i].indices_cpu,
      	              Matrices[i].values_cpu,
      	                    matPartialResult.ptr_cpu,
      	                    matPartialResult.indices_cpu,
                            matPartialResult.values_cpu);
        matPartialResult.numNNZ = matPartialResult.ptr_cpu[matPartialResult.numRows];

        matResult.numRows = matPartialResult.numRows;
        matResult.numCols = P.numCols;
        matResult.ptr_cpu = new size_t[matResult.numRows+1];
        csr_matmat_pass1(matPartialResult.numRows,
                      P.numCols,
                      matPartialResult.ptr_cpu,
                      matPartialResult.indices_cpu,
                      P.ptr_cpu,
                      P.indices_cpu,
                            matResult.ptr_cpu);
        matResult.numNNZ = matResult.ptr_cpu[matResult.numRows];
        matResult.indices_cpu = new size_t[matResult.numNNZ];
        matResult.values_cpu = new double[matResult.numNNZ];
        csr_matmat_pass2(matPartialResult.numRows,
                      P.numCols,
                      matPartialResult.ptr_cpu,
                      matPartialResult.indices_cpu,
      	              matPartialResult.values_cpu,
      	              P.ptr_cpu,
                      P.indices_cpu,
      	              P.values_cpu,
      	                    matResult.ptr_cpu,
      	                    matResult.indices_cpu,
      	                    matResult.values_cpu);
        matResult.numNNZ = matResult.ptr_cpu[matResult.numRows];

        sortMatrix(matResult, true);
        if(matResult.numRows >= min_system_size || i+1 == max_levels-1 ){
            mallocAndCopyMem(matResult.indices_cpu, matResult.indices_gpu, matResult.numNNZ);
            mallocAndCopyMem(matResult.ptr_cpu, matResult.ptr_gpu, matResult.numRows+1);
            mallocAndCopyMem(matResult.values_cpu, matResult.values_gpu, matResult.numNNZ);
        }
        clock_t t5 = clock();
        /** Store new matrix in the next lvl of hierarchy, i+1 **/
        Matrices[i+1] = matResult;

        /** Free memory from useless structures **/
        delete[] matPartialResult.indices_cpu;
        delete[] matPartialResult.ptr_cpu;
        delete[] matPartialResult.values_cpu;

/*        cout << "Level " << i << endl;
        cout << "   Time to create diagonal " << double(t2-t1) / CLOCKS_PER_SEC << "s" << endl;
        cout << "   Time to create P " << double(t3-t2) / CLOCKS_PER_SEC << "s" << endl;
        cout << "   Time to create Q " << double(t4-t3) / CLOCKS_PER_SEC << "s" << endl;
        cout << "   Time to create mat " << double(t5-t4) / CLOCKS_PER_SEC << "s" << endl;*/
    }
    return i;
}

_Matrix generateP_vCPU(const _Matrix& A, const _Vector& diag, double W, size_t numLevelsRoh){
    _Matrix Ptent;
    createPTent(A, Ptent);

    _Matrix P;
    /** Create P from Ptent **/
    //wDA
    _Matrix prodMat;
    prodMat.numNNZ = A.numNNZ;
    prodMat.numRows = A.numRows;
    prodMat.numCols = A.numCols;

    prodMat.indices_cpu = new size_t[prodMat.numNNZ];
    prodMat.values_cpu = new double[prodMat.numNNZ];
    prodMat.ptr_cpu = new size_t[prodMat.numRows+1];

    prodMat.ptr_cpu[0] = 0;
    size_t currentIndice = 0;
    for(size_t i = 0; i < A.numRows; i++){
        prodMat.ptr_cpu[i+1] = prodMat.ptr_cpu[i];
        if(diag.values_cpu[i] != 0.0){
            for(size_t r = A.ptr_cpu[i]; r < A.ptr_cpu[i+1]; r++){
                prodMat.ptr_cpu[i+1]++;
                prodMat.indices_cpu[currentIndice] = A.indices_cpu[r];
                prodMat.values_cpu[currentIndice] = A.values_cpu[r] * diag.values_cpu[i];
                currentIndice++;
            }
        }
    }

    prodMat.numNNZ = prodMat.ptr_cpu[prodMat.numRows];
    mallocAndCopyMem(prodMat.indices_cpu, prodMat.indices_gpu, prodMat.numNNZ);
    mallocAndCopyMem(prodMat.values_cpu, prodMat.values_gpu, prodMat.numNNZ);
    mallocAndCopyMem(prodMat.ptr_cpu, prodMat.ptr_gpu, prodMat.numRows+1);

    double roh_ = roh(prodMat, numLevelsRoh);
    double W_ = W / roh_;
    //cout << "La W queda aixi: " << W_ << ", i la W es: " << W << endl;

    for(size_t i = 0; i < prodMat.numNNZ; i++){
        prodMat.values_cpu[i] *= W_;
    }

    //printMatrix(prodMat);

    //I - wDA
    _Matrix subMat;
    subMat.numNNZ = prodMat.numRows + prodMat.numNNZ;
    subMat.numRows = prodMat.numRows;
    subMat.numCols = prodMat.numCols;

    subMat.ptr_cpu = new size_t[subMat.numRows+1];
    subMat.indices_cpu = new size_t[subMat.numNNZ];
    subMat.values_cpu = new double[subMat.numNNZ];
    subIdentityMatrix_cpu(prodMat, subMat);
    subMat.numNNZ = subMat.ptr_cpu[subMat.numRows];


    //(I - wDA) * PTent
    P.numRows = subMat.numRows;
    P.numCols = Ptent.numCols;
    P.ptr_cpu = new size_t[P.numRows+1];
    
    csr_matmat_pass1(subMat.numRows,
                      Ptent.numCols,
                      subMat.ptr_cpu,
                      subMat.indices_cpu,
                      Ptent.ptr_cpu,
                      Ptent.indices_cpu,
                            P.ptr_cpu);
    P.numNNZ = P.ptr_cpu[P.numRows];
    P.indices_cpu = new size_t[P.numNNZ];
    P.values_cpu = new double[P.numNNZ];
    csr_matmat_pass2(subMat.numRows,
                      Ptent.numCols,
                      subMat.ptr_cpu,
                      subMat.indices_cpu,
      	              subMat.values_cpu,
      	              Ptent.ptr_cpu,
                      Ptent.indices_cpu,
      	              Ptent.values_cpu,
      	                    P.ptr_cpu,
                            P.indices_cpu,
                            P.values_cpu);
    P.numNNZ = P.ptr_cpu[P.numRows];
    sortMatrix(P, true);


    mallocAndCopyMem(P.ptr_cpu, P.ptr_gpu, P.numRows+1);
    mallocAndCopyMem(P.indices_cpu, P.indices_gpu, P.numNNZ);
    mallocAndCopyMem(P.values_cpu, P.values_gpu, P.numNNZ);

    /** Free resources **/
    delete[] prodMat.ptr_cpu;
    delete[] prodMat.indices_cpu;
    delete[] prodMat.values_cpu;
    delete[] Ptent.indices_cpu;
    delete[] Ptent.ptr_cpu;
    delete[] Ptent.values_cpu;


    CUDA_CHECK(cudaFree(prodMat.ptr_gpu));
    CUDA_CHECK(cudaFree(prodMat.indices_gpu));
    CUDA_CHECK(cudaFree(prodMat.values_gpu));

    delete[] subMat.ptr_cpu;
    delete[] subMat.indices_cpu;
    delete[] subMat.values_cpu;

    /** Return P **/
    return P;
}

/** This is a simple function that transposes P assuming symmetric matrix
 *  Soon it will be need to implementate the non-symmetric construction of Q*/
_Matrix generateQ(const _Matrix& P){
    _Matrix Q;
    Q.numNNZ = P.numNNZ;
    Q.numCols = P.numRows;
    Q.numRows = P.numCols;
    Q.indices_cpu = new size_t[Q.numNNZ];
    Q.values_cpu = new double[Q.numNNZ];
    Q.ptr_cpu = new size_t[Q.numRows+1];
    malloc_(Q.indices_gpu, Q.numNNZ);
    malloc_(Q.values_gpu, Q.numNNZ);
    malloc_(Q.ptr_gpu, Q.numRows+1);
    //ens es suficient amb utilitzar la funcio csr_tocsc i prendre-ho com csr
    csr_tocsc(P.numRows,
	           P.numCols,
	           P.ptr_cpu,
	           P.indices_cpu,
	           P.values_cpu,
	                 Q.ptr_cpu,
	                 Q.indices_cpu,
	                 Q.values_cpu);
    //copy back
    copyMem(Q.values_cpu, Q.values_gpu, Q.numNNZ, 0);
    copyMem(Q.indices_cpu, Q.indices_gpu, Q.numNNZ, 0);
    copyMem(Q.ptr_cpu, Q.ptr_gpu, Q.numRows+1, 0);
    return Q;
}

/** This function will create G from diagonal of A **/
void createDiagonal(const _Matrix& A, _Vector& res){
    res.numElems = A.numRows;
    //malloc_(res.values_gpu, res.numElems);
    res.values_cpu = new double[res.numElems];
    malloc_(res.values_gpu, res.numElems);
    
    for(size_t i = 0; i < A.numRows; i++){
        res.values_cpu[i] = 0.0;
        for(size_t r = A.ptr_cpu[i]; r < A.ptr_cpu[i+1]; r++){
            if(A.indices_cpu[r] == i){
                res.values_cpu[i] =1.0/A.values_cpu[r];
                break;
            }
        }
    }
    //copia cpu->gpu
    copyMem(res.values_cpu, res.values_gpu, res.numElems, 0);
}

void createDiagonal_vCPU(const _Matrix& A, _Matrix&G, _Vector& res){
    res.numElems = A.numRows;
    //malloc_(res.values_gpu, res.numElems);
    res.values_cpu = new double[res.numElems];

    G.numRows = A.numRows;
    G.numCols = A.numCols;
    malloc_(G.values_gpu, G.numRows);
    malloc_(G.indices_gpu, G.numRows);
    malloc_(G.ptr_gpu, G.numRows+1);
    G.ptr_cpu = new size_t[G.numRows+1];
    G.indices_cpu = new size_t[G.numRows];
    G.values_cpu = new double[G.numRows];

    size_t currentIndice = 0;
    G.ptr_cpu[0] = 0;
    for(size_t i = 0; i < A.numRows; i++){
        G.ptr_cpu[i+1] = G.ptr_cpu[i];
        res.values_cpu[i] = 0.0;
        for(size_t r = A.ptr_cpu[i]; r < A.ptr_cpu[i+1]; r++){
            if(A.indices_cpu[r] == i){
                G.indices_cpu[currentIndice] = A.indices_cpu[r];
                res.values_cpu[i] = G.values_cpu[currentIndice] = 1.0/A.values_cpu[r];
                G.ptr_cpu[i+1]++;
                currentIndice++;
                break;
            }
        }
    }

    //copia cpu->gpu

    G.numNNZ = G.ptr_cpu[G.numRows];
    copyMem(G.ptr_cpu, G.ptr_gpu, G.numRows+1, 0);
    copyMem(G.values_cpu, G.values_gpu, G.numNNZ, 0);
    copyMem(G.indices_cpu, G.indices_gpu, G.numNNZ, 0);

    //copyMem(res.values_cpu, res.values_gpu, res.numElems, 0);
}
/** Aggregation of A to generate the colored graph **/
template <class I>
I standardAggregation(const I n_row,
                       const I Ap[],
                       const I Aj[],
                             I  x[])
{
    // Bj[n] == -1 means i-th node has not been aggregated
    std::fill(x, x + n_row, 0);

    I next_aggregate = 1; // number of aggregates + 1

    //Pass #1
    for(I i = 0; i < n_row; i++){
        if(x[i]){ continue; } //already marked

        const I row_start = Ap[i];
        const I row_end   = Ap[i+1];

        //Determine whether all neighbors of this node are free (not already aggregates)
        bool has_aggregated_neighbors = false;
        bool has_neighbors            = false;
        for(I jj = row_start; jj < row_end; jj++){
            const I j = Aj[jj];
            if( i != j ){
                has_neighbors = true;
                if( x[j] ){
                    has_aggregated_neighbors = true;
                    break;
                }
            }
        }

        if(!has_neighbors){
            //isolated node, do not aggregate
            x[i] = -n_row;
        }
        else if (!has_aggregated_neighbors){
            //Make an aggregate out of this node and its neighbors
            x[i] = next_aggregate;
            for(I jj = row_start; jj < row_end; jj++){
                x[Aj[jj]] = next_aggregate;
            }
            next_aggregate++;
        }
    }


    //Pass #2
    // Add unaggregated nodes to any neighboring aggregate
    for(I i = 0; i < n_row; i++){
        if(x[i]){ continue; } //already marked

        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];

            const I xj = x[j];
            if(xj > 0){
                x[i] = -xj;
                break;
            }
        }
    }

    next_aggregate--;

    //Pass #3
    for(I i = 0; i < n_row; i++){
        const I xi = x[i];

        if(xi != 0){
            // node i has been aggregated
            if(xi > 0)
                x[i] = xi - 1;
            else if(xi == -n_row)
                x[i] = -1;
            else
                x[i] = -xi - 1;
            continue;
        }

        // node i has not been aggregated
        const I row_start = Ap[i];
        const I row_end   = Ap[i+1];

        x[i] = next_aggregate;

        for(I jj = row_start; jj < row_end; jj++){
            const I j = Aj[jj];

            if(x[j] == 0){ //unmarked neighbors
                x[j] = next_aggregate;
            }
        }
        next_aggregate++;
    }

    return next_aggregate; //number of aggregates
}

/** This function will create Ptent from a given A **/
void createPTent(const _Matrix& A, _Matrix& P){
    long *x;
    long size = A.numRows;
    x = new long[size];
    size_t maxColumn = standardAggregation(size, (long*)A.ptr_cpu, (long*)A.indices_cpu, x);
    P.ptr_cpu = new size_t[size+1];
    P.indices_cpu = new size_t[size];
    P.values_cpu = new double[size];

    long* auxValues = new long[size];
    //ini auxValues to 0
    for(size_t i = 0; i < size; i++){
        auxValues[i] = 0;
    }
    //calculate each value group num
    for(size_t i = 0; i < size; i++){
        auxValues[x[i]]++;
    }
    for(size_t i = 0; i < size; i++){
        if(x[i] >= size)
            cout << "a l'index " << i << " tenim un valor de " << x[i]  << ", i el size es: " << size << endl;
    }
    //assign indices and right values
    P.ptr_cpu[0] = 0;
    for(size_t i = 0; i < size; i++){
        size_t j = x[i];
        P.ptr_cpu[i+1] = P.ptr_cpu[i] + 1;
        P.indices_cpu[i] = j;
        P.values_cpu[i] = 1.0/sqrt(auxValues[j]);
    }

    P.numCols = maxColumn;
    P.numRows = size;
    P.numNNZ = size;

    delete[] x;
    delete[] auxValues;
}

void csr_tocsc(const size_t n_row,
	           const size_t n_col,
	           const size_t Ap[],
	           const size_t Aj[],
	           const double Ax[],
	                 size_t*& Bp,
	                 size_t* Bi,
	                 double* Bx)
{
    const size_t nnz = Ap[n_row];
    Bp = new size_t[n_col+1];
    //compute number of non-zero entries per column of A
    std::fill(Bp, Bp + n_col, 0);

    for (size_t n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(size_t col = 0, cumsum = 0; col < n_col; col++){
        size_t temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz;

    for(size_t row = 0; row < n_row; row++){
        for(size_t jj = Ap[row]; jj < Ap[row+1]; jj++){
            size_t col  = Aj[jj];
            size_t dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    for(size_t col = 0, last = 0; col <= n_col; col++){
        size_t temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }
}

double calculateNorm(_Vector& b){
    double finalNum = 0.0;
    for(size_t i = 0; i < b.numElems; i++){
        finalNum += pow(b.values_cpu[i], 2);
    }
    return sqrt(finalNum);
}

double calculateNorm_GPU(_Vector&b){
    return cublasDnrm2(b.numElems, b.values_gpu, 1);
}

double checkResidual(const _Vector& u, const _Vector& b, const _Matrix& A){
    /** Au **/
    _Vector auxAU;
    auxAU.numElems = A.numRows;
    malloc_(auxAU.values_gpu, auxAU.numElems);
    dim3 Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    GPU_MatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.numRows, A.indices_gpu,
            A.ptr_gpu, A.values_gpu, u.values_gpu, auxAU.values_gpu);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 1208" << endl;
    }
    /** b - AU **/
    _Vector auxABU;
    auxABU.numElems = auxAU.numElems;
    auxABU.values_cpu = new double[auxABU.numElems];
    malloc_(auxABU.values_gpu, auxABU.numElems);
    Grid = Build_Grid(A.numRows, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.values_gpu, auxAU.values_gpu,
            auxABU.values_gpu, auxABU.numElems);
    if(!CUDA_Success(cudaThreadSynchronize())){
        cout << "Error en linea 1219" << endl;
    }
    copyMem(auxABU.values_gpu, auxABU.values_cpu, auxABU.numElems, 1);
    double finalNum = calculateNorm_GPU(auxABU);
    delete[] auxABU.values_cpu;
    CUDA_CHECK(cudaFree(auxAU.values_gpu));
    CUDA_CHECK(cudaFree(auxABU.values_gpu));

    return finalNum;
}

bool checkConvergence(const _Vector& u, const _Vector& b, const _Matrix& A, const double lastResidual, const double threshold){
    //double threshold = 1e-9;
    double newResidual = checkResidual(u, b, A);
    cout << "New norm: " << newResidual << endl;
    cout << "Current status: " << newResidual/lastResidual << endl;
/*    for(size_t i = 0; i < 5; i++){
        cout << "Element " << i << ": " << u.values_cpu[i] << endl;
    }*/
    cout << endl;
    if(newResidual/lastResidual < threshold){
        return true;
    }
    return false;

}

void computeDenseMatrix(const _Matrix& A, double *& vec){
    vec = new double[A.numCols * A.numRows];
    size_t pointer = 0;
    size_t currentIndice = 0;
    for(size_t i = 0; i < A.numRows; i++){
        size_t nonZeros = A.ptr_cpu[i+1] - A.ptr_cpu[i];
        size_t columnPointer = 0;
        for(size_t j = 0; j < A.numCols; j++){
            if(nonZeros > 0  && A.indices_cpu[currentIndice] == columnPointer){
                vec[pointer] = (double)A.values_cpu[currentIndice];
                currentIndice++;
                nonZeros--;
            }else{
                vec[pointer] = 0.0;
            }
            pointer++;
            columnPointer++;
        }
    }
}

void subIdentityMatrix_cpu(const _Matrix& A, _Matrix& sub){
    size_t currentIndex = 0;
    sub.ptr_cpu[0] = 0;
    for(size_t i = 0; i < A.numRows; i++){
        bool haveDiagonal = false;
        size_t numElems = 0;
        long lastIndex = -1;
        for(size_t r = A.ptr_cpu[i]; r < A.ptr_cpu[i+1]; r++){
            lastIndex = A.indices_cpu[r];
            if(A.indices_cpu[r] > i && !haveDiagonal){
                sub.indices_cpu[currentIndex] = i;
                sub.values_cpu[currentIndex] = 1.0;
                currentIndex++;
                haveDiagonal = true;
                numElems++;
            }
            if(A.indices_cpu[r] == i){
                haveDiagonal = true;
                if(A.values_cpu[r] != 1.0){
                    sub.indices_cpu[currentIndex] = A.indices_cpu[r];
                    sub.values_cpu[currentIndex] = 1.0 - A.values_cpu[r];
                    currentIndex++;
                    numElems++;
                }
            }else{
                sub.indices_cpu[currentIndex] = A.indices_cpu[r];
                sub.values_cpu[currentIndex] = -A.values_cpu[r];
                currentIndex++;
                numElems++;
            }
        }
        if(lastIndex < i){
            sub.indices_cpu[currentIndex] = i;
            sub.values_cpu[currentIndex] = 1.0;
            currentIndex++;
            numElems++;
        }
        sub.ptr_cpu[i+1] = sub.ptr_cpu[i] + numElems;
    }
}

double eigVals(_Vector& H, size_t finalIters){
    LaGenMatDouble A( H.values_cpu, finalIters, finalIters);
    LaVectorDouble eigvals_real(finalIters);
    LaVectorDouble eigvals_imag(finalIters);
    LaGenMatDouble VR(finalIters, finalIters);

    LaEigSolve(A, eigvals_real,
		eigvals_imag, VR);
    double *real = eigvals_real.addr();
    double *imag = eigvals_imag.addr();

    double max = real[0];
    if(max < 0)
        max = -max;
    for(size_t i = 1; i < finalIters; i++){
        if(real[i] < 0)
            real[i] = -real[i];
        if(real[i] > max)
            max = real[i];
    }

 /*   cout << "ROH VALUES " << endl;
    for(size_t i = 0; i < finalIters; i++)
        cout << real[i] << " ";
    cout << endl;*/
    return max;
}

double roh(const _Matrix& A, size_t iter){
    double threshold = 1e-10;

    size_t maxIter;
    maxIter = A.numCols < iter ?  A.numCols : iter;

    _Vector *V = new _Vector[maxIter+1];

    V[0].numElems = A.numCols;
    V[0].values_cpu = new double[V[0].numElems];
    srand(0);
    for(size_t i = 0; i < V[0].numElems; i++){
        V[0].values_cpu[i] = (double)(((int)rand())%100000000)/100000000.0;
        //cout << V[0].values_cpu[i]<< endl;
    }
    double v0Norm = calculateNorm(V[0]);
    for(size_t i = 0; i < V[0].numElems; i++){
        V[0].values_cpu[i] /= v0Norm;
    }
    mallocAndCopyMem(V[0].values_cpu, V[0].values_gpu, V[0].numElems);
    //delete[] V[0].values_cpu;

    _Vector H;
    H.numElems = (maxIter+1) * (maxIter+1);
    H.values_cpu = new double[H.numElems];
    for(size_t q = 0; q < H.numElems; q++){
        H.values_cpu[q] = 0.0;
    }
    size_t numCurrentV = 1;

    size_t j;

    for(j = 0; j < maxIter; j++){
        V[numCurrentV].numElems = A.numRows;
        malloc_(V[numCurrentV].values_gpu, V[numCurrentV].numElems);
        V[numCurrentV].values_cpu = new double[V[numCurrentV].numElems];
        dim3 grid = Build_Grid(V[numCurrentV].numElems, BLOCK_SIZE);
        GPU_MatrixVectorMultiply_CSR_Kernel <<< grid, BLOCK_SIZE >>> (A.numRows, A.indices_gpu, A.ptr_gpu,
            A.values_gpu, V[numCurrentV-1].values_gpu, V[numCurrentV].values_gpu);
        if(!CUDA_Success(cudaThreadSynchronize())){
            cout << "Error en linea 1379" << endl;
        }
        copyMem(V[numCurrentV].values_gpu, V[numCurrentV].values_cpu, V[numCurrentV].numElems, 1);

        _Vector auxVec;
        auxVec.numElems = V[numCurrentV].numElems;
        malloc_(auxVec.values_gpu, auxVec.numElems);

        grid = Build_Grid(V[numCurrentV-1].numElems, BLOCK_SIZE);
        for(size_t i = 0; i < numCurrentV; i++){
            size_t matrixIndice = (i*(maxIter+1))+j;
            double auxVal = H.values_cpu[matrixIndice] = cublasDdot(V[i].numElems, V[i].values_gpu, 1, V[numCurrentV].values_gpu, 1);
            subVectorConstantValue <<< grid, BLOCK_SIZE>>> (V[numCurrentV].values_gpu, H.values_cpu[matrixIndice], V[i].values_gpu, V[numCurrentV].numElems);
            if(!CUDA_Success(cudaThreadSynchronize())){
                cout << "Error en linea 1394" << endl;
            }
        }
        size_t matrixIndice = ((j+1) * (maxIter+1)) + j;
        copyMem(V[numCurrentV].values_gpu, V[numCurrentV].values_cpu, V[numCurrentV].numElems, 1);
        //copyMem(H.values_gpu, H.values_cpu, H.numElems, 1);
        H.values_cpu[matrixIndice] = calculateNorm(V[numCurrentV]);
        if(H.values_cpu[matrixIndice] < threshold)
            break;
        //copyMem(H.values_cpu, H.values_gpu, H.numElems, 0);

        divideVectorConstantValue <<< grid, BLOCK_SIZE >>> (V[numCurrentV].values_gpu, H.values_cpu[matrixIndice], V[numCurrentV].numElems);
        if(!CUDA_Success(cudaThreadSynchronize())){
            cout << "Error en linea 1407" << endl;
        }
        CUDA_CHECK(cudaFree(auxVec.values_gpu));

        numCurrentV++;

    }
    for(size_t i = 0; i < numCurrentV; i++){
        delete[] V[i].values_cpu;
        CUDA_CHECK(cudaFree(V[i].values_gpu));
    }
    delete[] V;

/*    for(size_t q = 0; q < maxIter+1; q++){
        for(size_t w = 0; w < maxIter+1; w++){
            cout << H.values_cpu[(q*(maxIter+1)) + w] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

    double max = eigVals(H, maxIter+1);

    delete[] H.values_cpu;

    return max;
}

/**     Memory management functions     **/
/*template <class Q>
void mallocAndCopyMem(Q*& CPU, Q*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(Q)));
    CUDA_CHECK(cudaMemcpy(GPU, CPU, size*sizeof(Q), cudaMemcpyHostToDevice));
}*/

void mallocAndCopyMem(double*& CPU, double*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(double)));
    CUDA_CHECK(cudaMemcpy(GPU, CPU, size*sizeof(double), cudaMemcpyHostToDevice));
}

void mallocAndCopyMem(size_t*& CPU, size_t*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(size_t)));
    CUDA_CHECK(cudaMemcpy(GPU, CPU, size*sizeof(size_t), cudaMemcpyHostToDevice));
}

/*template <class Q>
void malloc_(Q*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(Q)));
}*/

void malloc_(double*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(double)));
}

void malloc_(size_t*& GPU, size_t size){
    CUDA_CHECK(cudaMalloc((void**) &GPU, size*sizeof(size_t)));
}

template <class Q>
void copyMem(Q*& source, Q*& destiny, size_t size, unsigned short way){
    switch(way){
        case 0:
            CUDA_CHECK(cudaMemcpy(destiny, source, size*sizeof(Q), cudaMemcpyHostToDevice));
            break;
        case 1:
            CUDA_CHECK(cudaMemcpy(destiny, source, size*sizeof(Q), cudaMemcpyDeviceToHost));
            break;
        case 2:
            CUDA_CHECK(cudaMemcpy(destiny, source, size*sizeof(Q), cudaMemcpyDeviceToDevice));
            break;
    }
}

void deletingStuff(size_t* stuff){
    CUDA_CHECK(cudaFree(stuff));
}

void deletingStuff(double* stuff){
    CUDA_CHECK(cudaFree(stuff));
}

void GPU_VectorMultiply(double* sourceVec, double* destinyVec, size_t N){
	dim3 Grid = Build_Grid(N, BLOCK_SIZE);
	GPU_VectorVectorMultiplyElementWise_Kernel <<<Grid, BLOCK_SIZE>>> (N, sourceVec, destinyVec, destinyVec);
	GPUSparse::CUDA_Success(cudaThreadSynchronize());
}





}

}

//
// Compilation command
// make
