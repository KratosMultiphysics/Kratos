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
#include <stdio.h>

// Includes, project

#include "gpu_sparse.h"
#include "gpu_sparse_utils.h"
#include "gpu_sparse_kernels.h"

//Includes, preconditioner
#include <lapackd.h>
#include <laslv.h>
#include <gmd.h>
#include <gmf.h>


using namespace std;

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
	haveDenseRepresentation = false;

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

GPUCSRMatrix::GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2): NNZ(_NNZ), Size1(_Size1), Size2(_Size2), CPU_Columns(0), CPU_RowIndices(0), CPU_Values(0), GPU_Columns(0), GPU_RowIndices(0), GPU_Values(0), Allocated(false), haveDenseRepresentation(false)
{
	// Nothing to do!
}

GPUCSRMatrix::~GPUCSRMatrix()
{
	// Free CSR data; as it is allocated in one chunk of memory, we need only to free the begining address
	cudaFreeHost(CPU_Values);

	if (Allocated)
		GPU_Free();
	if (haveDenseRepresentation)
		delete[] matAuxValues;
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

bool GPUCSRMatrix::GenerateDenseRepresentation(bool FortranRep){
	haveDenseRepresentation = true;
	numValuesDenseRep = Size2 * Size1;
	matAuxValues = new double[Size2 * Size1];
	size_t pointer = 0;
	size_t currentIndice = 0;
	for(size_t i = 0; i < Size1; i++){
		size_t nonZeros = CPU_RowIndices[i+1] - CPU_RowIndices[i];
		
		for(size_t j = 0; j < Size2; j++){

		    if(nonZeros > 0  && CPU_Columns[currentIndice] == j){
			matAuxValues[pointer] = (double)CPU_Values[currentIndice];
			currentIndice++;
			nonZeros--;
		    }else{
			matAuxValues[pointer] = 0.0;
		    }
		    pointer++;
		}
		
	}
	return true;
}

// Operations defined on GPUCSRMatrix and GPUVector

//
// CPUGPUCSRMatrixVectorMultiply
// Matrix-Vector multiply on CPU

bool CPUGPUCSRMatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y)
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
// GPUGPUCSRMatrixVectorMultiply
// Matrix-Vector multiply on GPU

bool GPUGPUCSRMatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y)
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
		GPUGPUCSRMatrixVectorMultiply_CSRGPUVectorized_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values, Y.GPU_Values);

		if (!GPUSparse::CUDA_Success(cudaGetLastError()))
			return false;
	}

	else

	{
		dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
		GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values, Y.GPU_Values);

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
// CPUGPUCSRMatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on CPU

bool CPUGPUCSRMatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X)
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
// GPUGPUCSRMatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on GPU

bool GPUGPUCSRMatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X)
{
	// Primary checks
	if (A.Size1 != A.Size2 || A.Size2 != X.Size || !A.Allocated || !X.Allocated)
		return false;

	dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPUGPUCSRMatrixGetDiagonals_Kernel <<<Grid, BLOCK_SIZE>>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values, X.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPUGPUCSRMatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on CPU

bool CPUGPUCSRMatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A)
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
// GPUGPUCSRMatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on GPU

bool GPUGPUCSRMatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A)
{
	// Primary checks
	if (X.Size != A.Size1 || !X.Allocated || !A.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);
	GPUGPUCSRMatrixMatrixDiagonalMultiply_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values, A.GPU_Columns, A.GPU_RowIndices, A.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPUGPUVectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on CPU

bool CPUGPUVectorPrepareDiagonalPreconditionerValues(GPUVector &X)
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
// GPUGPUVectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on GPU

bool GPUGPUVectorPrepareDiagonalPreconditionerValues(GPUVector &X)
{
	// Primary check
	if (!X.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);
	GPUGPUVectorPrepareDiagonalPreconditionerValues_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values);
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
// CPUGPUVectorVectorMultiply
// Vector-Vector multiply on CPU

bool CPUGPUVectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result)
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
// GPUGPUVectorVectorMultiply
// Vector-Vector multiply on GPU

bool GPUGPUVectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result)
{
	// Primary check
	if (X.Size != Y.Size || !X.Allocated || !Y.Allocated)
		return false;

	Result = cublasDdot(X.Size, X.GPU_Values, 1, Y.GPU_Values, 1);

	return CUBLAS_Success(cublasGetError());
}

//
// CPUGPUVectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on CPU

bool CPUGPUVectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y,  GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Z.CPU_Values[i] = X.CPU_Values[i] * Y.CPU_Values[i];

	return true;
}

//
// GPUGPUVectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on GPU

bool GPUGPUVectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size || !X.Allocated || !Y.Allocated || !Z.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	GPUGPUVectorVectorMultiplyElementWise_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, X.GPU_Values, Y.GPU_Values, Z.GPU_Values);
	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

//
// CPUGPUVectorNorm2
// Vector norm 2 on CPU

bool CPUGPUVectorNorm2(GPUVector &X, double &Result)
{
	Result = static_cast <double> (0);

	for (size_t i = 0; i < X.Size; i++)
		Result += X.CPU_Values[i] * X.CPU_Values[i];

	Result = sqrt(Result);

	return true;
}

//
// GPUGPUVectorNorm2
// Vector norm 2 on GPU

bool GPUGPUVectorNorm2(GPUVector &X, double &Result)
{
	// Primary check
	if (!X.Allocated)
		return false;

	Result = cublasDnrm2(X.Size, X.GPU_Values, 1);

	return CUBLAS_Success(cublasGetError());
}

//
// CPUGPUVectorScaleAndAdd
// Vector scale-and-add on CPU

// Variant 1: Z = A * X + B * Y

bool CPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Z.CPU_Values[i] = A * X.CPU_Values[i] + B * Y.CPU_Values[i];

	return true;
}

// Variant 2: Y = A * X + B * Y

bool CPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y)
{
	// Primary check
	if (X.Size != Y.Size)
		return false;

	for (size_t i = 0; i < X.Size; i++)
		Y.CPU_Values[i] = A * X.CPU_Values[i] + B * Y.CPU_Values[i];

	return true;
}

//
// GPUGPUVectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool GPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
{
	// Primary check
	if (X.Size != Y.Size || Y.Size != Z.Size || !X.Allocated || !Y.Allocated || !Z.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	if (A == 1.00)
	{
	//printf("A = 1.0\n");
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_1_A_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_B_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_1_E_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else if (A == -1.00)
	{
	//printf("A = -1.0\n");
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_1_C_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_D_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_1_F_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	else
	{
	//printf("B = 1.0\n");
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_1_G_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_H_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_1_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}

	if (!GPUSparse::CUDA_Success(cudaGetLastError()))
		return false;

	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
}

// Variant 2: Y = A * X + B * Y

bool GPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y)
{
	// Primary check
	if (X.Size != Y.Size || !X.Allocated || !Y.Allocated)
		return false;

	dim3 Grid = Build_Grid(X.Size, BLOCK_SIZE);

	if (A == 1.00)
	{
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_2_A_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_2_B_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_2_E_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
	}

	else if (A == -1.00)
	{
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_2_C_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_2_D_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_2_F_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
	}

	else
	{
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_2_G_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_2_H_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);

		else
			GPUGPUVectorScaleAndAdd_2_Kernel <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values);
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

bool GPUGPUVectorScaleAndAdd_addingVersion(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z)
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
			GPUGPUVectorScaleAndAdd_1_A_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);

			
		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_B_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
				
		else
			GPUGPUVectorScaleAndAdd_1_E_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	
	else if (A == -1.00)
	{
	//printf("A = -1.0\n");
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_1_C_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_D_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else
			GPUGPUVectorScaleAndAdd_1_F_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	
	else
	{
	//printf("B = 1.0\n");
		if (B == 1.00)
			GPUGPUVectorScaleAndAdd_1_G_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else if (B == -1.00)
			GPUGPUVectorScaleAndAdd_1_H_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
			
		else
			GPUGPUVectorScaleAndAdd_1_Kernel_addingVersion <<<Grid, BLOCK_SIZE>>> (X.Size, A, X.GPU_Values, B, Y.GPU_Values, Z.GPU_Values);
	}
	return GPUSparse::CUDA_Success(cudaThreadSynchronize());
	
}

/** temp variables for LU decomposition **/
LaVectorLongInt ipiv;
typedef struct{
    size_t Size;
    double* CPU_Values;
    double* GPU_Values;
}_Vector;

/** functions from scipy for mat-mat calculation **/
template <class I>
void csr_matmat_pass1(const I n_row,
                      const I n_col,
                      const I Ap[],
                      const I Aj[],
                      const I Bp[],
                      const I Bj[],
                            I Cp[]){
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
void multilevel(GPUCSRMatrix**& A, GPUCSRMatrix**& P, GPUCSRMatrix**& R, GPUCSRMatrix**& G, GPUVector& b, GPUVector& u,
			unsigned short lvl, unsigned short maxLevels, size_t* preSweeps, size_t* postSweeps, bool assumeZeros)
{
    bool vectorized = (A[lvl]->NNZ / A[lvl]->Size2) > 10;

    //calculateInstantVector(u, b, A[lvl], G[lvl]);
    if(lvl < maxLevels){

	double* vecR = new double[R[lvl]->Size1];
    	GPUVector r(R[lvl]->Size1, vecR);
	r.GPU_Allocate();
	
	//clock_t t1 = clock();
        if(assumeZeros) //we receive from the upper level a zero start vector
        {
            if(preSweeps[lvl] != 0)
            {
                //first iteration (does not require computation of residual
                dim3 Grid = Build_Grid(A[lvl]->Size1, BLOCK_SIZE);
                GPUGPUCSRMatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BLOCK_SIZE >>>(G[lvl]->Size1, G[lvl]->GPU_Columns,
		    G[lvl]->GPU_RowIndices, G[lvl]->GPU_Values, b.GPU_Values, u.GPU_Values);
                if(!CUDA_Success(cudaThreadSynchronize())){
                    printf("Error en linea 130\n");
                }
                //from the second sweel on we need to recompute the residual
                for(size_t i = 1; i < preSweeps[lvl]; i++){
		    if(!vectorized)
			calculateInstantVector(u, b, *A[lvl], *G[lvl]);
		    else
			calculateInstantVectorGPUVectorized(u, b, *A[lvl], *G[lvl]);
		}

                if(!vectorized)
                    generateResidual(*R[lvl], b, *A[lvl], u, r);
                else
                    generateResidualGPUVectorized(*R[lvl], b, *A[lvl], u, r);


            }
            else //preSweeps[0] == 0 case
            {
                //inefficient! -- in this case we do not need to recompute the residual
                if(!vectorized)
                    generateResidual(*R[lvl], b, *A[lvl], u, r);
                else
                    generateResidualGPUVectorized(*R[lvl], b, *A[lvl], u, r);

            }
        }
        else
        {
          //from the second sweel on we need to recompute the residual
                for(size_t i = 0; i < preSweeps[lvl]; i++){
		    if(!vectorized)
			calculateInstantVector(u, b, *A[lvl], *G[lvl]);
		    else
			calculateInstantVectorGPUVectorized(u, b, *A[lvl], *G[lvl]);
		}


                if(!vectorized)
                    generateResidual(*R[lvl], b, *A[lvl], u, r);
                else
                    generateResidualGPUVectorized(*R[lvl], b, *A[lvl], u, r);
        }
	double *vecV = new double[r.Size];
        GPUVector v(r.Size, vecV);
        v.GPU_Allocate();
        dim3 Grid = Build_Grid(v.Size, BLOCK_SIZE);
        fillWithZeros <<< Grid, BLOCK_SIZE >>>(v.GPU_Values, v.Size);
        if(!CUDA_Success(cudaThreadSynchronize())){
                printf("Error en linea 160\n");
        }

        multilevel(A, P, R, G, r, v, lvl+1, maxLevels, preSweeps, postSweeps, assumeZeros);

        GPUVector pv(P[lvl]->Size1);
	pv.GPU_Allocate();
        Grid = Build_Grid(pv.Size, BLOCK_SIZE);
        // here product matrix P with vector v
        GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(P[lvl]->Size1,
                P[lvl]->GPU_Columns, P[lvl]->GPU_RowIndices, P[lvl]->GPU_Values, v.GPU_Values, pv.GPU_Values);
        if(!CUDA_Success(cudaThreadSynchronize())){
            printf("Error en linea 173\n");
        }
        // here the addition of pv to u
        Grid = Build_Grid(u.Size, BLOCK_SIZE);
        sumVectorVector <<< Grid, BLOCK_SIZE >>> (pv.GPU_Values, u.GPU_Values, u.Size);
        if(!CUDA_Success(cudaThreadSynchronize())){
            printf("Error in the line 180");
        }
        //delete de pv, v i r
        delete[] vecV;
        delete[] vecR;

        //double norm2 = checkResidual(u, b, A[lvl]);

	for(size_t i = 0; i < postSweeps[lvl]; i++){
	    if(!vectorized)
		calculateInstantVector(u, b, *A[lvl], *G[lvl]);
	    else
		calculateInstantVectorGPUVectorized(u, b, *A[lvl], *G[lvl]);
	}

    }else{
	//clock_t t1 = clock();
        //here lapack direct solver

	u.Copy(GPU_CPU);
	b.Copy(GPU_CPU);


        LaGenMatDouble _A(A[lvl]->matAuxValues, A[lvl]->Size1, A[lvl]->Size2);
        LaGenMatDouble _b(b.CPU_Values, b.Size, 1);
        LaGenMatDouble _x(u.CPU_Values, u.Size, 1);
//	LaLinearSolve( _A, _x, _b );

    	_x.inject(_b);            // will throw exception if not conformant

	integer info = 0;
	int M = _A.size(0);
	integer Ml = M;
	integer lda = _A.inc(0) * _A.gdim(0);

	integer K = _x.size(1);
	integer ldx = _x.inc(0) * _x.gdim(0);
	F77NAME(dgetrs) ("No transpose", &Ml, &K, &_A(0,0), &lda, &ipiv(0), &_x(0,0), &ldx, &info);

	//int res = clapack_dgetrs(CblasRowMajor, CblasNoTrans, &Ml, &K, &_A(0,0), &lda, &ipiv(0), &_x(0,0), &ldx);

	//copyMem(u.CPU_Values, u.GPU_Values, u.Size, 0);
	u.Copy(CPU_GPU);

	//clock_t t2 = clock();

   }

}
/** This function is a wrapper for u += G ( A, b, u) **/
void calculateInstantVector(GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const GPUCSRMatrix& G)
{
		/** Au **/
	GPUVector auxAU(A.Size1);
	auxAU.GPU_Allocate();
	dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.Size1, A.GPU_Columns,
	    A.GPU_RowIndices, A.GPU_Values, u.GPU_Values, auxAU.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 238");
	}

		/** b - AU **/
	GPUVector auxABU(A.Size1);
	auxABU.GPU_Allocate();
	Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	subVectorVector <<< Grid, BLOCK_SIZE >>>(b.GPU_Values, auxAU.GPU_Values, auxABU.GPU_Values, auxABU.Size);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 249");
	}
		/** u += G ( b - Au ) **/
	Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPUGPUCSRMatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BLOCK_SIZE >>>(G.Size1, G.GPU_Columns,
	    G.GPU_RowIndices, G.GPU_Values, auxABU.GPU_Values, u.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 255");
	}
	//deleting structures

}
void calculateInstantVectorGPUVectorized(GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const GPUCSRMatrix& G)
{
		/** Au **/
	GPUVector auxAU(A.Size1);
	auxAU.GPU_Allocate();
	dim3 Grid = Build_Grid(A.Size1 * HWS, BS);
	GPUGPUCSRMatrixVectorMultiply_CSRGPUVectorized_Kernel <<< Grid, BS >>>(A.Size1, A.GPU_Columns,
	    A.GPU_RowIndices, A.GPU_Values, u.GPU_Values, auxAU.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 272");
	}
		/** b - AU **/
	GPUVector auxABU(A.Size1);
	auxABU.GPU_Allocate();
	Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	subVectorVector <<< Grid, BLOCK_SIZE >>>(b.GPU_Values, auxAU.GPU_Values, auxABU.GPU_Values, auxABU.Size);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 282");
	}
		/** u += G ( b - Au ) **/
	Grid = Build_Grid(A.Size1, BS);
	GPUGPUCSRMatrixVectorMultiply_CSR_Kernel_addingVersion <<< Grid, BS >>>(G.Size1, G.GPU_Columns,
	    G.GPU_RowIndices, G.GPU_Values, auxABU.GPU_Values, u.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
		printf("Error in the line 288");
	}
	//deleting structures

}

void generateResidual(const GPUCSRMatrix& R, const GPUVector& b, const GPUCSRMatrix& A, const GPUVector& u, GPUVector& r){
		/** Au **/
	GPUVector auxAU(A.Size1);
	auxAU.GPU_Allocate();
	dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.Size1, A.GPU_Columns,
	    A.GPU_RowIndices, A.GPU_Values, u.GPU_Values, auxAU.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
		printf("Error in the line 305");
	}
		/** b - AU **/
	GPUVector auxABU(A.Size1);
	auxABU.GPU_Allocate();
	Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	subVectorVector <<< Grid, BLOCK_SIZE >>>(b.GPU_Values, auxAU.GPU_Values,
	    auxABU.GPU_Values, auxABU.Size);
	if(!CUDA_Success(cudaThreadSynchronize())){
		printf("Error in the line 315");
	}
		/** r = R ( b - Au ) **/
	Grid = Build_Grid(R.Size1, BLOCK_SIZE);
	GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(R.Size1, R.GPU_Columns,
	    R.GPU_RowIndices, R.GPU_Values, auxABU.GPU_Values, r.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
		printf("Error in the line 325");
	}
    //deleting structures

}
void generateResidualGPUVectorized(const GPUCSRMatrix& R, const GPUVector& b, const GPUCSRMatrix& A, const GPUVector& u, GPUVector& r){
		/** Au **/
		GPUVector auxAU(A.Size1);
	auxAU.GPU_Allocate();
	dim3 Grid = Build_Grid(A.Size1 * HWS, BS);
	GPUGPUCSRMatrixVectorMultiply_CSRGPUVectorized_Kernel <<< Grid, BS >>>(A.Size1, A.GPU_Columns,
	    A.GPU_RowIndices, A.GPU_Values, u.GPU_Values, auxAU.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 340");
	}
		/** b - AU **/
	GPUVector auxABU(A.Size1);
	auxABU.GPU_Allocate();
	Grid = Build_Grid(A.Size1, BLOCK_SIZE);
	subVectorVector <<< Grid, BLOCK_SIZE >>>(b.GPU_Values, auxAU.GPU_Values,
	    auxABU.GPU_Values, auxABU.Size);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 350");
	}
		/** r = R ( b - Au ) **/
	Grid = Build_Grid(R.Size1 *  HWS, BS);
	GPUGPUCSRMatrixVectorMultiply_CSRGPUVectorized_Kernel <<< Grid, BS >>>(R.Size1, R.GPU_Columns,
	    R.GPU_RowIndices, R.GPU_Values, auxABU.GPU_Values, r.GPU_Values);
	if(!CUDA_Success(cudaThreadSynchronize())){
	printf("Error in the line 360");
	}
	//deleting structures

}

/** This function will return the hierarchy reconstructed of A and b;
    Additionally return the number of real hierarchy levels.
    All matrices and vectors are returned as GPU structures**/
size_t generateHierarchy(GPUCSRMatrix**& Matrices, GPUCSRMatrix**& Pmat, GPUCSRMatrix**& Qmat,
        GPUCSRMatrix**& Gmat, double W, size_t numLevelsRoh, size_t max_levels, size_t min_system_size)
{
    size_t i = 0;
    for(i = 0; i < max_levels; i++){
        /** This condition controls MAX_SYSTEM_SIZE for the last matrix in hierarchy **/
        if(Matrices[i]->Size1 < min_system_size || i == max_levels-1){
	    Matrices[i]->GenerateDenseRepresentation();
            LaGenMatDouble A(Matrices[i]->matAuxValues, Matrices[i]->Size1, Matrices[i]->Size2);

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

		
/*	    if( i > 0 ){ //A matrix is not the first one in hierarchy
		    if(Matrices[i].CPU_RowIndices != NULL) delete[] Matrices[i].CPU_Columns;
		    if(Matrices[i].CPU_Values != NULL) delete[] Matrices[i].CPU_Values;
		    if(Matrices[i].CPU_Columns != NULL) delete[] Matrices[i].CPU_RowIndices;
	    }*/
            Matrices[i]->numValuesDenseRep = A.inc(0) * A.gdim(0);

            break;
        }

	double *values, *values_2;
	size_t *ptr, *indices, *ptr_2, *indices_2;
	size_t NNZ, NNZ_2, Size1, Size2;
	double *vecDiag;
        createDiagonal_vCPU(*Matrices[i], NNZ, Size1, Size2, indices, ptr, values, vecDiag);
        GPUVector diag(Matrices[i]->Size1, vecDiag);
	Gmat[i] = new GPUCSRMatrix(NNZ, Size1, Size2, indices, ptr, values, false);
	Gmat[i]->GPU_Allocate();
	Gmat[i]->Copy(CPU_GPU, false);
	delete[] values;
	delete[] ptr;
	delete[] indices;
	
        /** Generating P and Q for the current A level **/
	
        generateP_vCPU(*Matrices[i], diag, W, numLevelsRoh, NNZ, Size1, Size2, indices, ptr, values);
	Pmat[i] = new GPUCSRMatrix(NNZ, Size1, Size2, indices, ptr, values, false);
	Pmat[i]->GPU_Allocate();
	Pmat[i]->Copy(CPU_GPU, false);
	delete[] values;
	delete[] ptr;
	delete[] indices;
        generateQ(*Pmat[i], NNZ, Size1, Size2, indices, ptr, values);
        Qmat[i] = new GPUCSRMatrix(NNZ, Size1, Size2, indices, ptr, values, false);
	Qmat[i]->GPU_Allocate();
	Qmat[i]->Copy(CPU_GPU, false);
	delete[] values;
	delete[] ptr;
	delete[] indices;
        //CUDA_CHECK(cudaFree(diag.GPU_Values));
	delete[] diag.CPU_Values;

        /** Allocating result matrix, and partialResult matrix **/

	//R * A = AUX
	ptr = new size_t[Qmat[i]->Size1+1];

        csr_matmat_pass1(Qmat[i]->Size1,
                      Matrices[i]->Size2,
                      Qmat[i]->CPU_RowIndices,
                      Qmat[i]->CPU_Columns,
                      Matrices[i]->CPU_RowIndices,
                      Matrices[i]->CPU_Columns,
                            ptr);
        NNZ = ptr[Qmat[i]->Size1];
        indices = new size_t[NNZ];
        values = new double[NNZ];
        csr_matmat_pass2(Qmat[i]->Size1,
      	              Matrices[i]->Size2,
      	              Qmat[i]->CPU_RowIndices,
                      Qmat[i]->CPU_Columns,
      	              Qmat[i]->CPU_Values,
      	              Matrices[i]->CPU_RowIndices,
                      Matrices[i]->CPU_Columns,
      	              Matrices[i]->CPU_Values,
      	                    ptr,
      	                    indices,
                            values);
        NNZ = ptr[Qmat[i]->Size1];
	//newA = AUX * P
        ptr_2 = new size_t[Qmat[i]->Size1+1];
        csr_matmat_pass1(Qmat[i]->Size1,
                      Pmat[i]->Size2,
                      ptr,
                      indices,
                      Pmat[i]->CPU_RowIndices,
                      Pmat[i]->CPU_Columns,
                            ptr_2);
        NNZ_2 = ptr_2[Qmat[i]->Size1];
        indices_2 = new size_t[NNZ_2];
        values_2 = new double[NNZ_2];
        csr_matmat_pass2(Qmat[i]->Size1,
                      Pmat[i]->Size2,
                      ptr,
                      indices,
      	              values,
      	              Pmat[i]->CPU_RowIndices,
                      Pmat[i]->CPU_Columns,
      	              Pmat[i]->CPU_Values,
      	                    ptr_2,
      	                    indices_2,
      	                    values_2);
        NNZ_2 = ptr_2[Qmat[i]->Size1];

	/** Store new matrix in the next lvl of hierarchy, i+1 **/
	Matrices[i+1] = new GPUCSRMatrix(NNZ_2, Qmat[i]->Size1, Pmat[i]->Size2, indices_2, ptr_2, values_2, false);

        sortMatrix(Matrices[i+1]->CPU_RowIndices, Matrices[i+1]->CPU_Columns, Matrices[i+1]->CPU_Values, Matrices[i+1]->Size1, Matrices[i+1]->NNZ, true);
        if(Matrices[i+1]->Size1 >= min_system_size || i+1 == max_levels-1 ){
		Matrices[i+1]->GPU_Allocate();
		Matrices[i+1]->Copy(CPU_GPU, false);
        }

        /** Free memory from useless structures **/
	delete[] values;
	delete[] ptr;
	delete[] indices;
	delete[] values_2;
	delete[] ptr_2;
	delete[] indices_2;

    }
    return i;
}

void generateP_vCPU(const GPUCSRMatrix& A, const GPUVector& diag, double W, size_t numLevelsRoh, size_t &pNNZ, size_t &pSize1, size_t &pSize2, size_t *&pindices, size_t *&pptr, double *&pvalues){
    size_t NNZ_ptent, Size1_ptent, Size2_ptent, *ptr_ptent, *indices_ptent;
    double *values_ptent;
    createPTent(A, NNZ_ptent, Size1_ptent, Size2_ptent, ptr_ptent, indices_ptent, values_ptent);

    /** Create P from Ptent **/
    //wDA

    size_t NNZ_prod, *ptr_prod, *indices_prod;
    double *values_prod;

    indices_prod = new size_t[A.NNZ];
    values_prod = new double[A.NNZ];
    ptr_prod = new size_t[A.Size1+1];

    ptr_prod[0] = 0;
    size_t currentIndice = 0;
    for(size_t i = 0; i < A.Size1; i++){
        ptr_prod[i+1] = ptr_prod[i];
        if(diag.CPU_Values[i] != 0.0){
            for(size_t r = A.CPU_RowIndices[i]; r < A.CPU_RowIndices[i+1]; r++){
                ptr_prod[i+1]++;
                indices_prod[currentIndice] = A.CPU_Columns[r];
                values_prod[currentIndice] = A.CPU_Values[r] * diag.CPU_Values[i];
                currentIndice++;
            }
        }
    }

    NNZ_prod = ptr_prod[A.Size1];
    GPUCSRMatrix prodMat(NNZ_prod, A.Size1, A.Size2, indices_prod, ptr_prod, values_prod);
    delete[] indices_prod;
    delete[] values_prod;
    delete[] ptr_prod;
    prodMat.GPU_Allocate();
    prodMat.Copy(CPU_GPU, false);

    double roh_ = roh(prodMat, numLevelsRoh);
    double W_ = W / roh_;

    for(size_t i = 0; i < prodMat.NNZ; i++){
        prodMat.CPU_Values[i] *= W_;
    }

    //I - wDA

    size_t subNNZ = prodMat.Size1 + prodMat.NNZ;
    size_t* ptr_sub = new size_t[prodMat.Size1+1];
    size_t* indices_sub = new size_t[subNNZ];
    double* values_sub = new double[subNNZ];

    subIdentityMatrix_cpu(prodMat, ptr_sub, indices_sub, values_sub);

    //(I - wDA) * PTent
    pSize1 = prodMat.Size1;
    pSize2 = Size2_ptent;
    pptr = new size_t[pSize1+1];
    
    csr_matmat_pass1(prodMat.Size1,
                      Size2_ptent,
                      ptr_sub,
                      indices_sub,
                      ptr_ptent,
                      indices_ptent,
                            pptr);
    pNNZ = pptr[pSize1];
    pindices = new size_t[pNNZ];
    pvalues = new double[pNNZ];
    csr_matmat_pass2(prodMat.Size1,
                      Size2_ptent,
                      ptr_sub,
                      indices_sub,
      	              values_sub,
      	              ptr_ptent,
                      indices_ptent,
      	              values_ptent,
      	                    pptr,
                            pindices,
                            pvalues);
    pNNZ = pptr[pSize1];

    sortMatrix(pptr, pindices, pvalues, pSize1, pNNZ, true);

    /** Free resources **/
    delete[] values_ptent;
    delete[] indices_ptent;
    delete[] ptr_ptent;

    delete[] indices_sub;
    delete[] values_sub;
    delete[] ptr_sub;

}

/** This is a simple function that transposes P assuming symmetric matrix
 *  Soon it will be need to implementate the non-symmetric construction of Q*/
void generateQ(const GPUCSRMatrix &P, size_t &NNZ, size_t &Size1, size_t &Size2, size_t *&indices, size_t *&ptr, double *&values){
    NNZ = P.NNZ;
    Size2 = P.Size1;
    Size1 = P.Size2;
    indices = new size_t[NNZ];
    values = new double[NNZ];
    ptr = new size_t[Size1+1];
    //ens es suficient amb utilitzar la funcio csr_tocsc i prendre-ho com csr
    csr_tocsc(P.Size1,
	           P.Size2,
	           P.CPU_RowIndices,
	           P.CPU_Columns,
	           P.CPU_Values,
	                 ptr,
	                 indices,
	                 values);
}

/** This function will create G from diagonal of A **/
void createDiagonal(const GPUCSRMatrix &A, size_t &Size, double *&values){
    Size = A.Size1;
    values = new double[Size];

    for(size_t i = 0; i < A.Size1; i++){
        values[i] = 0.0;
        for(size_t r = A.CPU_RowIndices[i]; r < A.CPU_RowIndices[i+1]; r++){
            if(A.CPU_Columns[r] == i){
		if(fabs(A.CPU_Values[r]) > 1e-30)
                	values[i] =1.0/A.CPU_Values[r];
		else{
			values[i]  =1.0;
		}
	
                break;
            }
        }
    }

}

void createDiagonal_vCPU(const GPUCSRMatrix& A, size_t &NNZ, size_t &Size1, size_t &Size2, size_t *&indices, size_t *&ptr, double *&values, double *&vecDiag){

    vecDiag = new double[A.Size1];

    Size1 = A.Size1;
    Size2 = A.Size2;

    ptr = new size_t[Size1+1];
    indices = new size_t[Size1];
    values = new double[Size1];

    size_t currentIndice = 0;
    ptr[0] = 0;
    for(size_t i = 0; i < A.Size1; i++){
        ptr[i+1] = ptr[i];
        vecDiag[i] = 0.0;
        for(size_t r = A.CPU_RowIndices[i]; r < A.CPU_RowIndices[i+1]; r++){
            if(A.CPU_Columns[r] == i){
                indices[currentIndice] = A.CPU_Columns[r];
                vecDiag[i] = values[currentIndice] = 1.0/A.CPU_Values[r];
                ptr[i+1]++;
                currentIndice++;
                break;
            }
        }
    }
    NNZ = ptr[Size1];
}
/** Aggregation of A to generate the colored graph **/
template <class I>
I standardAggregation(const I n_row,
                       const I Ap[],
                       const I Aj[],
                             I  x[])
{
    // Bj[n] == -1 means i-th node has not been aggregated
    //std::fill(x, x + n_row, 0);
    for(long q = 0; q < n_row; q++)
		x[q] = 0;

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
void createPTent(const GPUCSRMatrix &A, size_t &NNZ_ptent, size_t &Size1_ptent, size_t &Size2_ptent, size_t *&ptr_ptent, size_t *&indices_ptent, double *&values_ptent){
    long *x;
    long size = A.Size1;
    x = new long[size];
    size_t maxColumn = standardAggregation(size, (long*)A.CPU_RowIndices, (long*)A.CPU_Columns, x);

    ptr_ptent = new size_t[size+1];
    indices_ptent = new size_t[size];
    values_ptent = new double[size];

    long* auxValues = new long[size];
    //ini auxValues to 0
    for(long i = 0; i < size; i++){
        auxValues[i] = 0;
    }
    //calculate each value group num
    for(long i = 0; i < size; i++){
        auxValues[x[i]]++;
    }

    //assign indices and right values
    ptr_ptent[0] = 0;
    for(long i = 0; i < size; i++){
        size_t j = x[i];
        ptr_ptent[i+1] = ptr_ptent[i] + 1;
        indices_ptent[i] = j;
        values_ptent[i] = 1.0/sqrt((double)auxValues[j]);
    }

    Size2_ptent = maxColumn;
    Size1_ptent = size;
    NNZ_ptent = size;

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
    //std::fill(Bp, Bp + n_col, 0);
	for(size_t q = 0; q < n_col; q++)
		Bp[q] = 0;

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

double calculateNorm(_Vector &b){
    double finalNum = 0.0;
    for(size_t i = 0; i < b.Size; i++){
        finalNum += pow(b.CPU_Values[i], 2);
    }
    return sqrt(finalNum);
}

double calculateNorm_GPU(GPUVector &b){
    return cublasDnrm2(b.Size, b.GPU_Values, 1);
}

double checkResidual(const GPUVector &u, const GPUVector &b, const GPUCSRMatrix &A){
    /** Au **/
    GPUVector auxAU(A.Size1);
    auxAU.GPU_Allocate();
    dim3 Grid = Build_Grid(A.Size1, BLOCK_SIZE);
    GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< Grid, BLOCK_SIZE >>>(A.Size1, A.GPU_Columns,
            A.GPU_RowIndices, A.GPU_Values, u.GPU_Values, auxAU.GPU_Values);
    if(!CUDA_Success(cudaThreadSynchronize())){
        printf("Error in the line 1208");
    }

    /** b - AU **/
    GPUVector auxABU(auxAU.Size);
    auxABU.GPU_Allocate();
    Grid = Build_Grid(A.Size1, BLOCK_SIZE);
    subVectorVector <<< Grid, BLOCK_SIZE >>>(b.GPU_Values, auxAU.GPU_Values,
            auxABU.GPU_Values, auxABU.Size);
    if(!CUDA_Success(cudaThreadSynchronize())){
        printf("Error in the line 1219");
    }
    double finalNum = calculateNorm_GPU(auxABU);

    return finalNum;
}

bool checkConvergence(const GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const double lastResidual, const double threshold){
    double newResidual = checkResidual(u, b, A);

    if(newResidual/lastResidual < threshold){
        return true;
    }
    return false;

}

void computeDenseMatrix(const GPUCSRMatrix& A, double *& vec){
    vec = new double[A.Size2 * A.Size1];
    size_t pointer = 0;
    size_t currentIndice = 0;
    for(size_t i = 0; i < A.Size1; i++){
        size_t nonZeros = A.CPU_RowIndices[i+1] - A.CPU_RowIndices[i];
        size_t columnPointer = 0;
        for(size_t j = 0; j < A.Size2; j++){
            if(nonZeros > 0  && A.CPU_Columns[currentIndice] == columnPointer){
                vec[pointer] = (double)A.CPU_Values[currentIndice];
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

void subIdentityMatrix_cpu(const GPUCSRMatrix& A, size_t* CPU_RowIndices, size_t* CPU_Columns, double* CPU_Values){
    size_t currentIndex = 0;
    CPU_RowIndices[0] = 0;
    for(size_t i = 0; i < A.Size1; i++){
        bool haveDiagonal = false;
        size_t numElems = 0;
        long lastIndex = -1;
        for(size_t r = A.CPU_RowIndices[i]; r < A.CPU_RowIndices[i+1]; r++){
            lastIndex = A.CPU_Columns[r];
            if(A.CPU_Columns[r] > i && !haveDiagonal){
                CPU_Columns[currentIndex] = i;
                CPU_Values[currentIndex] = 1.0;
                currentIndex++;
                haveDiagonal = true;
                numElems++;
            }
            if(A.CPU_Columns[r] == i){
                haveDiagonal = true;
                if(A.CPU_Values[r] != 1.0){
                    CPU_Columns[currentIndex] = A.CPU_Columns[r];
                    CPU_Values[currentIndex] = 1.0 - A.CPU_Values[r];
                    currentIndex++;
                    numElems++;
                }
            }else{
                CPU_Columns[currentIndex] = A.CPU_Columns[r];
                CPU_Values[currentIndex] = -A.CPU_Values[r];
                currentIndex++;
                numElems++;
            }
        }
        if(lastIndex < (long)i){
            CPU_Columns[currentIndex] = i;
            CPU_Values[currentIndex] = 1.0;
            currentIndex++;
            numElems++;
        }
        CPU_RowIndices[i+1] = CPU_RowIndices[i] + numElems;
    }
}

double eigVals(_Vector& H, size_t finalIters){
    LaGenMatDouble A( H.CPU_Values, finalIters, finalIters);
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
    return max;
}

double roh(const GPUCSRMatrix& A, size_t iter){
    double threshold = 1e-10;

    size_t maxIter;
    maxIter = A.Size2 < iter ?  A.Size2 : iter;

    _Vector *V = new _Vector[maxIter+1];

    V[0].Size = A.Size2;
    V[0].CPU_Values = new double[V[0].Size];
    srand(0);
    for(size_t i = 0; i < V[0].Size; i++){
        V[0].CPU_Values[i] = (double)(((int)rand())%100000000)/100000000.0;

    }
    double v0Norm = calculateNorm(V[0]);
    for(size_t i = 0; i < V[0].Size; i++){
        V[0].CPU_Values[i] /= v0Norm;
    }
    mallocAndCopyMem(V[0].CPU_Values, V[0].GPU_Values, V[0].Size);
    //delete[] V[0].CPU_Values;

    _Vector H;
    H.Size = (maxIter+1) * (maxIter+1);
    H.CPU_Values = new double[H.Size];
    for(size_t q = 0; q < H.Size; q++){
        H.CPU_Values[q] = 0.0;
    }
    size_t numCurrentV = 1;

    size_t j;

    for(j = 0; j < maxIter; j++){
        V[numCurrentV].Size = A.Size1;
        malloc_(V[numCurrentV].GPU_Values, V[numCurrentV].Size);
        V[numCurrentV].CPU_Values = new double[V[numCurrentV].Size];
        dim3 grid = Build_Grid(V[numCurrentV].Size, BLOCK_SIZE);
        GPUGPUCSRMatrixVectorMultiply_CSR_Kernel <<< grid, BLOCK_SIZE >>> (A.Size1, A.GPU_Columns, A.GPU_RowIndices,
            A.GPU_Values, V[numCurrentV-1].GPU_Values, V[numCurrentV].GPU_Values);
        if(!CUDA_Success(cudaThreadSynchronize())){
            printf("Error in line 1379\n");
        }
        copyMem(V[numCurrentV].GPU_Values, V[numCurrentV].CPU_Values, V[numCurrentV].Size, 1);

        _Vector auxVec;
        auxVec.Size = V[numCurrentV].Size;
        malloc_(auxVec.GPU_Values, auxVec.Size);

        grid = Build_Grid(V[numCurrentV-1].Size, BLOCK_SIZE);
        for(size_t i = 0; i < numCurrentV; i++){
            size_t matrixIndice = (i*(maxIter+1))+j;
            double auxVal = H.CPU_Values[matrixIndice] = cublasDdot(V[i].Size, V[i].GPU_Values, 1, V[numCurrentV].GPU_Values, 1);
            subVectorConstantValue <<< grid, BLOCK_SIZE>>> (V[numCurrentV].GPU_Values, H.CPU_Values[matrixIndice], V[i].GPU_Values, V[numCurrentV].Size);
            if(!CUDA_Success(cudaThreadSynchronize())){
                printf("Error in line 1394\n");
            }
        }
        size_t matrixIndice = ((j+1) * (maxIter+1)) + j;
        copyMem(V[numCurrentV].GPU_Values, V[numCurrentV].CPU_Values, V[numCurrentV].Size, 1);
        //copyMem(H.GPU_Values, H.CPU_Values, H.Size, 1);
        H.CPU_Values[matrixIndice] = calculateNorm(V[numCurrentV]);
        if(H.CPU_Values[matrixIndice] < threshold)
            break;
        //copyMem(H.CPU_Values, H.GPU_Values, H.Size, 0);

        divideVectorConstantValue <<< grid, BLOCK_SIZE >>> (V[numCurrentV].GPU_Values, H.CPU_Values[matrixIndice], V[numCurrentV].Size);
        if(!CUDA_Success(cudaThreadSynchronize())){
            printf("Error in line 1407\n");
        }
        CUDA_CHECK(cudaFree(auxVec.GPU_Values));

        numCurrentV++;

    }
    for(size_t i = 0; i < numCurrentV; i++){
        delete[] V[i].CPU_Values;
        CUDA_CHECK(cudaFree(V[i].GPU_Values));
    }
    delete[] V;

    double max = eigVals(H, maxIter+1);

    delete[] H.CPU_Values;

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

void GPUGPUVectorMultiply(double* sourceVec, double* destinyVec, size_t N){
	dim3 Grid = Build_Grid(N, BLOCK_SIZE);
	GPUGPUVectorVectorMultiplyElementWise_Kernel <<<Grid, BLOCK_SIZE>>> (N, sourceVec, destinyVec, destinyVec);
	GPUSparse::CUDA_Success(cudaThreadSynchronize());
}





}

}

//
// Compilation command
// make
