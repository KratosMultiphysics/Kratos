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

// Includes, project

#include "gpu_sparse.h"
#include "gpu_sparse_utils.h"
#include "gpu_sparse_kernels.h"

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
		CPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		CPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		CPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)

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
		CPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		CPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		CPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)

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
		GPU_Values = reinterpret_cast <double *> (CSR_Data);  // Size: NNZ * sizeof(double)
		GPU_Columns = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * sizeof(double));  // Size: NNZ * sizeof(size_t)
		GPU_RowIndices = reinterpret_cast <size_t *> (reinterpret_cast <size_t> (CSR_Data) + NNZ * (sizeof(double) + sizeof(size_t)));  // Size: (Size1 + 1) * sizeof(size_t)

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

	// Free CSR data; as it is allocated in one chunk of memory, we need only to free the begining address
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
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_A_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_B_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPU_VectorScaleAndAdd_1_E_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else if (A == -1.00)
	{
		if (B == 1.00)
			GPU_VectorScaleAndAdd_1_C_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPU_VectorScaleAndAdd_1_D_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPU_VectorScaleAndAdd_1_F_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else
	{
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

}

}

//
// Compilation command
// make
